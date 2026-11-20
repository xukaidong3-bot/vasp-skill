# VASP 编译选项、补丁与并行环境（00-build-env）

> 什么时候读这个文件：用户问「VASP 装好没有 / 用 vasp_std 还是 vasp_gam / 为什么并行不加速 / 报 `ifort: command not found` / 怎么开 HDF5、VTST、VASPsol / GPU 版参数」，或要给作业脚本定进程数、KPAR、NCORE 时。

## 1. 先判断走哪条路线

| 情形 | 做法 | 关键区别 |
|---|---|---|
| 集群已有 module（`vasp/6.6.1`） | 不编译，直接跳到 §2.4 空跑 + §2.5 选并行参数 | 先 `ldd $(which vasp_std) | grep "not found"` 确认为空 |
| Intel CPU 集群自编 | `cp arch/makefile.include.oneapi makefile.include`（ifx + Intel MPI + MKL） | oneAPI 2025 已移除 ifort，必须走 ifx；ifx 模板 6.4.3 起才有 |
| AMD EPYC / 个人电脑 / WSL2 | `gnu_omp`（GCC+OpenBLAS+ScaLAPACK+FFTW）或 `gnu_ompi_aocl`、`aocc_ompi_aocl` | 系统包要注释 `OPENBLAS_ROOT`/`SCALAPACK_ROOT`/`FFTW_ROOT`，改 `BLASPACK = -lopenblas`、`SCALAPACK = -lscalapack-openmpi` |
| NVIDIA GPU | `nvhpc_ompi_mkl_omp_acc`（Intel 宿主）或 `nvhpc_acc`；NVHPC ≥ 21.2，避开 22.1/22.2 | 保留模板默认 `-DACC_OFFLOAD -DNVCUDA -DUSENCCL -Dqd_emulate`；1 rank/GPU；`-gpu=cc80` |
| 要用 VTST（NEB/Dimer/IOPT）或 VASPsol | 打补丁后重编：`cp -r vtstcode6.6.1/* src/`；`cp VASPsol/src/solvation.F src/` + `CPP_OPTIONS += -Dsol_compat` | 补丁目录版本号必须与 VASP 完全相同，旧补丁配新源码是编译失败最常见原因 |

可执行程序判据：KPOINTS 只有 Γ 点（`Gamma 1 1 1`）→ `vasp_gam`；INCAR 有 `LSORBIT`/`LNONCOLLINEAR` → `vasp_ncl`；其余 `vasp_std`。用 `vasp_gam` 跑多 k 点会直接报错。

## 2. 标准流程

1. **选模板并备份**：`cp arch/makefile.include.oneapi makefile.include`；模板必须与源码同版本（6.5.0 把 `-D_OPENACC` 改成 `-DACC_OFFLOAD`）。
2. **只改 `## Customize as of this point!` 之后的内容**，最小改动如下（路径换成自己的）：
   ```makefile
   VASP_TARGET_CPU = -xHOST          # 只在"编译机就是计算节点且是 Intel"时保留
   CPP_OPTIONS += -DVASP_HDF5
   HDF5_ROOT   ?= $(HOME)/opt/hdf5-1.14.5-oneapi
   LLIBS       += -L$(HDF5_ROOT)/lib -lhdf5_fortran -Wl,-rpath,$(HDF5_ROOT)/lib
   INCS        += -I$(HDF5_ROOT)/include
   ```
   检查点：`ifx --version` 有输出、`echo $MKLROOT` 非空；HDF5 的 `hdf5.mod` 必须由**同一个 Fortran 编译器**编译。
3. **编译**：`make veryclean && make DEPS=1 -j 16 all 2>&1 | tee make.log`（`all` = std gam ncl）。
   检查点：`ls -lh bin/` 有三个程序；`grep -i -E " error |error:" make.log` 为空；`ldd bin/vasp_std` 无 `not found`。
4. **空跑拿并行规模**（不真算电子步，秒级）：
   ```bash
   cp INCAR INCAR.prod
   echo "ALGO = None" >> INCAR
   export OMP_NUM_THREADS=1; ulimit -s unlimited
   mpirun -np 64 vasp_std > dry.out 2>&1
   grep -E "NKPTS|NBANDS" OUTCAR | head -3
   mv INCAR.prod INCAR
   ```
   检查点：拿到不可约 k 点数 NKPTS 与 NBANDS，决定 KPAR 能整除谁。
5. **定并行参数并计时**：候选 3–4 组，`NELM = 6`、`NSW = 0` 各跑一遍，比 `LOOP:` 计时：
   ```bash
   grep "LOOP:" OUTCAR | awk '{s+=$NF; n++} END {print "mean real time per SCF step:", s/n, "s"}'
   ```
6. **验证安装**：`make test`（FAST 类，4 核约 1.5 h），`grep -c ERROR testsuite/testsuite.log` 应为 0。

## 3. 关键标签表

| 标签 / 变量 | 取值 | 何时用 / 注意 |
|---|---|---|
| `KPAR` | 默认 1；取能整除总进程数、并尽量整除 NKPTS 的值 | k 点并行，内存线性增长；让每个节点/NUMA 域各负责一组 |
| `NCORE` | 默认 1；≈ `sqrt(N_ranks/KPAR)`，常取一个 NUMA 域的核数（如 8） | 小体系/慢网络用 1；约 100 原子用 4，400 原子以上可 12–16 |
| `NSIM` | 默认 4；GPU 上 16–32 | CPU 上一般不改 |
| `NPAR` | 旧标签，`NPAR × NCORE = 可用进程数` | 不要与 NCORE 同时设置（同时出现时 NPAR 优先） |
| `OMP_NUM_THREADS` | 纯 MPI：`1`；GPU/混合：每 rank 线程数 | 必须显式 `export`，否则默认线程数导致严重超额订阅 |
| `OMP_STACKSIZE` | CPU `512m`；OpenMP offload GPU `2048m` | 配合 `ulimit -s unlimited` 防 SIGSEGV |
| `ALGO = None` | 空跑 | 只做 FFT 规划与 k 点约化，用来读 NKPTS/NBANDS |
| `LSOL`/`EB_K`/`TAU` | VASPsol 标签 | 只在打过 VASPsol 补丁并加 `-Dsol_compat` 的构建里有效 |
| `LCLIMB`/`IMAGES`/`ICHAIN` | VTST 标签 | 只在 `IMAGES`/`ICHAIN` 被设置时才调用 VTST，不影响其他计算 |

## 4. 易错点

1. **INCAR 里写了 `NCORE` 完全不生效** → 用了 OpenMP（`OMP_NUM_THREADS > 1`）或任何 GPU 版本时，VASP 自动把 NCORE 重置为 1。→ 纯 MPI 必须 `export OMP_NUM_THREADS=1`；GPU 版把这个变量当线程数用，不要再指望 NCORE。
2. **编译立刻失败 `ifort: command not found`** → oneAPI 2025 起不再提供 ifort，`mpiifort` 默认仍调 ifort。→ 用 6.4.3+ 的 `makefile.include.oneapi`，把 FC/FCL 改成 `mpiifort -fc=ifx` 或 `mpiifx`。
3. **链接报 `undefined reference to ...h5...`，或编译 HDF5 相关文件报无法读取 `hdf5.mod`** → `LLIBS` 缺 `-lhdf5_fortran`，或 HDF5 用了另一版本编译器编译。→ 检查 `$HDF5_ROOT/lib` 下的库名，并用编译 VASP 的同一 Fortran 编译器重编 HDF5。
4. **运行报 `error while loading shared libraries: libxxx.so: cannot open shared object file`** → 共享库不在 `LD_LIBRARY_PATH`。→ `ldd bin/vasp_std && grep "not found"` 定位，补 `LD_LIBRARY_PATH` 或链接时加 `-Wl,-rpath,...`。
5. **VTST 编译报 `CHAIN_FORCE` 参数不匹配 / 找不到 PyAMFF 对象** → `src/main.F` 没插 `TSIF,`，或 `src/.objects` 里 `dynmat.o neb.o dimer.o lanczos.o instanton.o` 没排在 `chain.o` 前，或 `src/makefile` 的 `LIB` 末尾缺 `pyamff_fortran`。→ 改完 `make veryclean` 重编；`nm bin/vasp_std | grep -i -c dimer` > 0 才算编进去。
6. **6.6.0 做 ISPIN=2 的 SCAN/r2SCAN 优化时结构"越优化越怪"** → known issue 94，自旋极化 meta-GGA 受力错误（能量正确）。→ 升级 6.6.1，或按官方说明给 6.6.0 打 `pot.F.patch` 后重编。
7. **纯 MPI 运行比预期慢很多** → `_omp` 模板未设 `OMP_NUM_THREADS=1`，每个 rank 按核数开线程造成超额订阅。→ 作业脚本里无条件 `export OMP_NUM_THREADS=1`。
8. **GPU 版比 CPU 还慢** → 每块 GPU 跑了多个 rank、体系太小、`NSIM` 太小。→ 1 rank/GPU，`KPAR` = GPU 数，`NSIM = 16–32`；小体系标准 DFT 不适合 GPU。

## 5. 算完了怎么验收

- `OUTCAR` 末尾有 `General timing and accounting informations for this job`：出现 = 正常结束；没有 = 被杀/崩溃，别读结果。
- 同一算例 1 核与 4 核：`grep "energy(sigma->0)" OUTCAR | tail -1` 之差应在 10⁻⁵ eV 量级（与 EDIFF 同量级）。差到 meV 以上说明 NBANDS 被并行度改动或 SCF 没收敛，固定 `NBANDS` 与并行设置重跑。
- `grep "NKPTS" OUTCAR`：与空跑一致，说明并行分组没把网格弄掉；改进程数时 VASP 可能自动增大 NBANDS 并打印 `The number of bands has been changed`，做 LOBSTER/LOPTICS/MLFF 必须固定并行设置并手写 `NBANDS`。
- `grep -c ERROR testsuite/testsuite.log` = 0；`mpirun --report-bindings` 中 rank 落在同一 NUMA 域内才正常，GPU 作业再用 `nvidia-smi` 确认每卡一个 rank。

## 6. 本地完整输入

- `$VASP_TUTORIALS_ROOT/tutorials/01-setup-compile.md`
- `$VASP_TUTORIALS_ROOT/examples/01-setup-compile/makefile.include.oneapi.changes`
- `$VASP_TUTORIALS_ROOT/examples/01-setup-compile/vtst_patch_steps.md`
- `$VASP_TUTORIALS_ROOT/examples/01-setup-compile/job_cpu.slurm`
- `$VASP_TUTORIALS_ROOT/examples/01-setup-compile/test_Si/INCAR`

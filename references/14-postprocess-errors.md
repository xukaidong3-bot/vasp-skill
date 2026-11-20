# 后处理生态与报错排查

> 什么时候读这个文件：要选后处理工具（VASPKIT / pymatgen / ASE / py4vasp）、批量提交与作业数组、加自动纠错（custodian / atomate2），或作业失败要按报错原文定位原因时。

## 1. 先判断走哪条路线

| 情形 | 做法 | 关键区别 |
|---|---|---|
| VASP ≥ 6.4 且编译了 `-DVASP_HDF5` | py4vasp 读 `vaspout.h5`，pymatgen/VASPKIT 交叉验证 | 输出最全（KPOINTS_OPT、部分电荷、声子、MLFF） |
| 无 HDF5 或 VASP 5.4.4 | pymatgen 读 vasprun.xml + VASPKIT 读 DOSCAR/PROCAR/EIGENVAL | 依赖文件完整，中断作业读不了 |
| 只拿带隙/DOS 数字 | VASPKIT 911 / 11x，或 `eigenvalue_band_properties` | 菜单式，零脚本 |
| 计算数 < 50 | Shell 循环 + Slurm 作业数组 + 检查脚本 | 手工可控 |
| 50–1000 / >1000 或多步依赖 | 作业数组 + custodian；atomate2 + jobflow-remote/FireWorks/AiiDA | 自动纠错会改参数，必须查 `custodian.json` |

## 2. 标准流程

1. **装齐并核对版本**：VASPKIT 1.5.1、pymatgen 2026.10.24（含 `pymatgen-core`）、ASE 3.29.0（`from ase.mep import NEB` 需 ≥ 3.23）、py4vasp 0.11.3（需 `-DVASP_HDF5`）。检查：`python -c "import pymatgen,ase,py4vasp"`。
2. **VASPKIT 配置**：`cp -f how_to_set_environment_variable ~/.vaspkit`，改 `PBE_PATH`、`RECOMMENDED_POTCAR .TRUE.`、`PLOT_MATPLOTLIB .TRUE.`。检查：`vaspkit -task 102 -kpr 0.04` 生成 KPOINTS。
3. **单算例解析**（三种入口）：

```python
from pymatgen.io.vasp import Vasprun        # parse_potcar_file=False 跳过 POTCAR 校验
v = Vasprun("vasprun.xml"); print(v.eigenvalue_band_properties, v.converged_electronic)
from ase.io import read; at = read("OUTCAR", index=":")   # XDATCAR 用 index="::10"
import py4vasp; py4vasp.Calculation.from_path(".").bandgap.print()
```

检查：`converged_electronic` 为 True 且 OUTCAR 有 `General timing`，否则结果不可信。
4. **批处理提数**：`for d in */; do (cd "$d" && echo -e "911\n" | vaspkit > gap.log); done`。**不要**按行号截 VASPKIT 屏幕输出；读它写出的 `BAND_GAP`、`PDOS_*.dat`。检查：一个算例一个目录（功能会写固定文件名）。
5. **Slurm 作业数组**（`sbatch --array=1-${N}%10 array_job.slurm`，`%10` 限并发）：

```bash
DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" dirs.txt); cd "$DIR" || exit 1
if [ -f OUTCAR ] && grep -q "General timing" OUTCAR; then exit 0; fi   # 重投安全
mpirun -np ${SLURM_NTASKS:-32} vasp_std > vasp.out 2>&1
```

脚本里要有 `ulimit -s unlimited`。检查：`grep -c "General timing" */OUTCAR` 全为 1。
6. **自动纠错**：`python custodian_run.py --cmd "mpirun -np 32 vasp_std" --mode double --max-errors 5`。检查：每批结束读 `custodian.json`。
7. **工作流**：`RelaxBandStructureMaker().make(struct)` + `run_locally(flow, create_folders=True)`；配 `~/.atomate2.yaml` 的 `VASP_CMD`（或 `ATOMATE2_VASP_CMD`）、`pmg config --add PMG_VASP_PSP_DIR`。检查：`job_<uuid>` 目录数 = 预期任务数。

## 3. 关键标签表

| 标签 | 取值 | 何时用 / 注意 |
|---|---|---|
| `ALGO` | `Normal` / `All` | `All` 治 ZHEGV、ZBRENT、EDDRMM；6.5.1 起可加 `ISEARCH = 1` |
| `ISYM` / `SYMPREC` | 0 / 1E-6 | 对称性内部错误：先改 Γ 中心网格，再减 `SYMPREC`（≤ 1E-4 内调），最后 `ISYM = 0` |
| `NCORE` / `KPAR` / `NPAR` | 进程数的因子 | ELF 需 `KPAR = 1` + `NPAR = 1`；DFPT 删 `NCORE`/`NPAR`；`M_divide` 查整除 |
| `LMAXMIX` | 4（d）/ 6（f） | DFT+U 混合不收敛时加 |
| `LSCALAPACK` | `.FALSE.` | `PSSYEVX`/`PDSYEVX` 失败时关 |
| 环境变量 | `PMG_VASP_PSP_DIR`、`ASE_VASP_COMMAND`、`VASP_PP_PATH`、`ATOMATE2_VASP_CMD`、`~/.vaspkit` | 各工具找 POTCAR/可执行文件方式不同，别只设一处 |

## 4. 易错点：报错原文 → 原因 → 修复

先定位：`grep -c "General timing" OUTCAR`（0 = 中途退出）、`tail -20 vasp.out`、`grep -n -i -E "error|very bad news" vasp.out`、`scripts/check_vasp.py -h`。

| 报错原文 / 现象 | 原因 | 修复 |
|---|---|---|
| `Error reading item '<tag>' from file INCAR` | 值格式错、Tab、CRLF、MAGMOM 个数不对 | `dos2unix INCAR`；空格代替 Tab；MAGMOM 写满 NIONS（非共线 3×NIONS） |
| 标签被静默忽略（OUTCAR 回显没有） | 拼写错、行首不可见字符、写在注释后 | `grep -i 标签名 OUTCAR` 核对 |
| `ERROR: charge density could not be read from file CHGCAR for ICHARG>10` | 缺 CHGCAR 或损坏 | 从自洽目录复制完整 CHGCAR |
| `WARNING: dimensions on CHGCAR file are different` | 网格与当前 ENCUT/PREC/晶胞不符 | 两步设置一致，或删旧 CHGCAR 重算 |
| `while reading WAVECAR` / `while reading plane` | WAVECAR 与当前 ENCUT/NBANDS/k 点/晶胞不匹配或损坏 | 删 WAVECAR 重头算；或保持与生成时一致 |
| `ERROR: the triple product of the basis vectors is negative` | 晶格矢量构成左手系 | 交换任意两个晶格矢量（同时交换坐标分量），或用 ASE/pymatgen 重写 POSCAR |
| `The distance between some ions is very small` | 原子重叠：建模错、Direct/Cartesian 单位写反 | VESTA 查结构；确认 POSCAR 第 8 行坐标类型 |
| 不报错但结果离谱（元素与 POTCAR 不匹配） | POTCAR 拼接顺序与 POSCAR 元素行不一致 | `grep TITEL POTCAR` 与 POSCAR 元素行逐一对照 |
| `VERY BAD NEWS! internal error in subroutine IBZKPT: Reciprocal lattice and k-lattice belong to different class of lattices` | 六角/三方体系用了偶数 Monkhorst–Pack 网格 | 改 Γ 中心网格 |
| `Tetrahedron method fails (number of k-points < 4)`、`Routine TETIRR needs special values` | `ISMEAR = -5` 但 k 点 < 4 或非 Γ 中心 | 加 k 点并 Γ 中心；分子/slab 改 `ISMEAR = 0`、`SIGMA = 0.05` |
| `internal error in subroutine PRICEL`、`Found some non-integer element in rotation matrix`、`SGRCON` | 结构对称性数值噪声；MP 网格破坏对称性 | spglib/vaspkit 602/603 整理结构；减 `SYMPREC`；`ISYM = 0` |
| `... rotation matrix was not found (increase SYMPREC)`、`POSMAP`/`RHOSYG`、`PRICELV: current lattice and primitive lattice are incommensurate` | 对称性容差不当；弛豫中对称性改变 | 先改 Γ 中心网格，再适度增大 `SYMPREC`（≤ 1E-4），最后 `ISYM = 0`；用 CONTCAR 重启 |
| `Inconsistent Bravais lattice types found for crystalline and reciprocal lattice` | 对称性容差问题 | `SYMPREC = 1E-6` 或更小；`ISYM = 0` |
| `WARNING: Sub-Space-Matrix is not hermitian in DAV` | 原子重叠、步长过大、`LREAL` 误差、并行过度 | 查结构；`LREAL = .FALSE.`；`PREC = Accurate`；`ALGO = All`；减核数 |
| `Error EDDDAV: Call to ZHEGV failed` | 矩阵不正定（原因同上） | `ALGO = All`（6.5.1 起可加 `ISEARCH = 1`）；删 CHGCAR/WAVECAR；减 POTIM |
| `WARNING in EDDRMM: call to ZHEGV failed` | RMM-DIIS（`ALGO = Fast/VeryFast`）不稳定 | `ALGO = Normal`；减 POTIM；从头开始 |
| `EDWAV: internal error, the gradient is not orthogonal` | 编译优化等级过高或 `ALGO` 组合不当 | 换 `ALGO = Normal`；降低编译优化等级；减核数 |
| `ERROR RSPHER`、`REAL_OPTLAY: internal error` | 实空间投影（`LREAL`）在部分体系/网格失败 | `LREAL = .FALSE.` |
| `ERROR in subspace rotation PSSYEVX` / `PDSYEVX` | ScaLAPACK 子空间旋转失败 | `ALGO = Normal`；调 NCORE；`LSCALAPACK = .FALSE.` |
| `BRMIX: very serious problems the old and the new charge density differ` | 读入不匹配 CHGCAR、对称性不一致、混合不稳 | 删 CHGCAR/WAVECAR；`ISYM = 0` 或调 `SYMPREC`；`IMIX = 1` |
| 到 NELM 仍不收敛（dE 振荡） | 金属/磁性混合难、初始磁矩离谱、`SIGMA` 太小 | 增 `NELM`；`ALGO = All`；调 `AMIX/BMIX`（磁性用 `AMIX_MAG/BMIX_MAG`）；DFT+U 加 `LMAXMIX` |
| `One of the lattice vectors is very long (>50 A), but AMIN is rather large` | 长真空轴 slab | `AMIN = 0.01` |
| `WARNING: Your highest band is occupied at some k-points!` | NBANDS 太少（金属或高温 smearing） | 增大 `NBANDS`（约 1.2–1.5 倍） |
| `Spin polarized Harris functional dynamics is a good joke` | `ICHARG ≥ 10`、`ISPIN = 2` 且 `NSW > 0` | 非自洽计算设 `NSW = 0` |
| `ZBRENT: fatal error in bracketing` | CG 线搜索失败：受力噪声大（SCF 精度不足） | `cp CONTCAR POSCAR` 续算；`EDIFF` 收到 1E-6、`NELMIN = 6–8`；减 `POTIM`；`IBRION = 1` 或 VTST FIRE（`IBRION = 3`、`POTIM = 0`、`IOPT = 7`） |
| `BRIONS problems: POTIM should be increased` | `IBRION = 1` 步长不合适 | 略增 POTIM，或改 `IBRION = 2` |
| `LAPACK: Routine ZPOTRF failed` | 原子过近（首步）或步长过大；小体系 NCORE 过大 | 查结构；减 POTIM；`ISYM = 0`；小体系减小 NCORE |
| 结构"乱飞"、能量突然升高 | POTIM 太大；MD 中含 H 体系 POTIM 过大 | 减 POTIM；MD 含 H 用 0.5–1 fs |
| NEB 图像间距异常、不收敛 | 初末态原子顺序不一致；进程数不能被 IMAGES 整除；过早开 CI | `dist.pl` 检查；进程数取 IMAGES 整数倍；先普通 NEB 再 `LCLIMB = .TRUE.` |
| 优化后体积系统性偏小/偏大 | Pulay 应力，ENCUT 太低 | `ENCUT ≥ 1.3×ENMAX`；用 CONTCAR 反复优化到体积不变 |
| `ELF: KPAR>1 not implemented` | ELF 不支持 k 点并行 | `KPAR = 1` 并设 `NPAR = 1` |
| `WARNING: ELF not implemented for non collinear case` | 非共线 | 用共线计算求 ELF |
| `PEAD routines do not work for NCORE` / `remove the tag NPAR from the INCAR file` | `LEPSILON`/`LCALCEPS`/DFPT 不支持能带并行 | 删除 `NCORE`、`NPAR` |
| DFPT（`IBRION = 7/8`）中 vdW 修正无效 | DFPT 不支持 `IVDW` | 改用有限位移（`IBRION = 5/6` 或 phonopy） |
| MLFF 续训报找不到 ML_AB / run 模式找不到 ML_FF | 没把 ML_ABN 拷成 ML_AB；没有 ML_FF | 续训前 `cp ML_ABN ML_AB`；`cp ML_FFN ML_FF` 后 `ML_MODE = run` |
| SOC 用 `vasp_std` 报错或 MAGMOM 个数错 | 非共线必须 `vasp_ncl`，MAGMOM 每原子 3 分量 | 用 `vasp_ncl`；MAGMOM 写 3×NIONS（6.5.1 前同时设 `ISPIN=2` 与 `LNONCOLLINEAR` 会报错） |
| 6.6.0 自旋极化 meta-GGA 受力错 | Known issue 94 | 升 6.6.1 或打官方 `pot.F.patch` 重编译 |
| `forrtl: severe (174): SIGSEGV, segmentation fault occurred` | 栈/内存不足、MPI 与编译器不匹配 | `ulimit -s unlimited`；减每节点进程；查 MPI 版本 |
| 作业被杀（OOM、`Killed`、`exit code 137`） | 内存不足 | 增 NCORE、减 KPAR；少放进程；降 NBANDS/ENCUT |
| `M_divide: can not subdivide N nodes by M` | NCORE/NPAR/KPAR 与进程数不整除 | 取进程数的因子 |
| GPU 版本很慢 | 每块 GPU 多个 MPI 进程；设了 NCORE>1 | 每块 GPU 1 个 rank 用 OpenMP 线程；`KPAR` = GPU 数；适当增大 `NSIM` |

## 5. 算完了怎么验收

- 正常结束：`grep -c "General timing" OUTCAR` = 1；`grep -l "reached required accuracy" */OUTCAR` 覆盖所有弛豫目录；`grep " F= " OSZICAR | tail -5` 看能量是否平稳。
- 收敛：`grep -E "^ *(DAV|RMM|CG|SDA) *:" OSZICAR | awk '{print $2}' | sort -n | tail -1` 未触 NELM。
- 结构：`grep -A3 "neighbor table" OUTCAR` 最近邻距离不像 0.5 Å；续算前确认 CONTCAR 非空完整（被杀时可能只写一半）。
- 可忽略：`WARNING: random wavefunctions but no delay for mixing`、第一个离子步电子步偏多、`The number of bands has been changed from X to Y`（LOBSTER/LOPTICS/MLFF 续训例外）。
- 结果可信：pymatgen `converged_electronic` 为 True 才用 band 属性；XML 解析失败 = 作业中断；批处理一致性靠 `custodian.json` 复核（它会把 `ISMEAR = -5` 改成 0、关对称性、换网格类型）。
- 合规：POTCAR、源码、`makefile.include` 模板、WAVECAR/CHGCAR/vaspout.h5 不入库（用 `.gitignore_vasp`）；误提交 POTCAR 需 `git filter-repo` 改写历史并强制更新远程；大文件 `xz -9` 压缩；按 VASP 2026-06-16 AI 规范，不得把源码与赝势库喂给外部模型。

## 6. 本地完整输入

- `$VASP_TUTORIALS_ROOT/examples/18-postprocess-toolkit/`
- `$VASP_TUTORIALS_ROOT/examples/18-postprocess-toolkit/error_table.md`
- `$VASP_TUTORIALS_ROOT/examples/18-postprocess-toolkit/array_job.slurm`
- `$VASP_TUTORIALS_ROOT/examples/18-postprocess-toolkit/custodian_run.py`
- `$VASP_TUTORIALS_ROOT/examples/18-postprocess-toolkit/atomate2_relax_band.py`
- `$VASP_TUTORIALS_ROOT/tutorials/18-postprocess-troubleshooting.md`

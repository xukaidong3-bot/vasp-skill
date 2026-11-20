---
name: vasp
description: >-
  Use for any VASP calculation task — 写/查 INCAR·POSCAR·KPOINTS·POTCAR、ENCUT/k 点/ISMEAR 收敛测试、结构优化与 EOS、能带与 DOS、HSE06/r2SCAN/mBJ、磁性/DFT+U/SOC、vdW 与二维材料、slab 与表面能/功函数、吸附能与自由能台阶图、NEB/Dimer 过渡态、声子、弹性常数/介电/光学、AIMD/机器学习力场、Bader/差分电荷/ELF/COHP、Wannier/有效质量/缺陷，以及读 VASP 报错。English: INCAR, convergence, band structure, DFT+U, NEB, phonon, Bader, MLFF, VASP errors.
whenToUse: >-
  用户要"用 VASP 算某个物理量"、要你写/查一组 INCAR、或问某个 VASP 报错时加载本 skill——它给参数决策、标准流程、验收判据与踩坑清单。只是把写好的作业送到远程服务器跑、排队、看进度，用 remote-dft-server。
---

# VASP 计算：参数决策与标准流程

本 skill 负责**算什么、标签怎么写、结果怎么判断**。精确参数与真实可运行输入不在这里，在本地讲义仓库里——需要照抄 INCAR 时去读 `examples/<task>/`，不要凭记忆编。

## 0. 用之前必须知道的四件事

1. **准备外部讲义（可选）**：本 skill 的示例路径指向外部的 VASP 讲义仓库。将该仓库 clone 到本地后设置：
   `export VASP_TUTORIALS_ROOT=/path/to/DFT-Tutorials/VASP`
   也可自行整理本地速查笔记；本公开仓库不包含私人笔记。
2. **要真跑就走远程服务器**：提交、POTCAR 拼接、排队、取结果全用 `remote-dft-server` skill 的 `qe_queue` / `qe_submit`（`potcarElements` 顺序 = POSCAR 物种顺序，`nproc:32`）。本仓库不含 POTCAR。
3. **版本门槛**（讲义按 6.4–6.6 写，当前 6.6.1）：6.6.0 的自旋极化 meta-GGA（SCAN/r2SCAN）受力是错的（known issue 94）→ 用 6.6.1；`LATTICE_CONSTRAINTS` 配 `IBRION=1/2` 要 ≥ 6.4.3；`ML_MODE` 要 ≥ 6.4.0；`CHECKPOINT_FD` 要 ≥ 6.6 且编译开 HDF5。远程机版本不确定时先看任一 `OUTCAR` 头部的版本行。
4. **不套默认值**：下面的骨架只解决常规体系。f 电子、带电体系、金属、强关联、低维，先按对应 reference 的判据确认参数再动手。

## 1. 任务 → 读哪个文件（一次只读一个）

| 要算的东西 / 遇到的问题 | 读 | 里面有什么 |
|---|---|---|
| 编译、补丁（VTST/VASPsol）、`KPAR`/`NCORE`/`NSIM` | `references/00-build-env.md` | makefile 模板与预编译宏、打补丁顺序、并行参数与线程超额订阅 |
| 输入文件写法与坑、输出看什么、收敛测试、续算 | `references/01-inputs.md` | 四件套逐项陷阱、`ISTART`/`ICHARG`、ENCUT 与 k 点收敛脚本 |
| 优化晶胞/离子、Pulay 应力、二维只松面内、分子、EOS | `references/02-structure.md` | `IBRION`/`ISIF`/`EDIFFG` 选择、CONTCAR 重启、Birch–Murnaghan |
| 能带、DOS/PDOS、带隙、HSE06、r2SCAN、mBJ | `references/03-electronic.md` | 两步法与 `KPOINTS_OPT`、高对称路径、`LORBIT=11`、零权重 k 点法 |
| 自旋极化、FM/AFM、DFT+U、线性响应、J 与 T_C、SOC、MAE | `references/04-magnetism.md` | `MAGMOM` 写法、AFM 超胞、U 的来源、力定理算 MAE |
| 色散修正、层间结合、真空层、异质结 | `references/05-vdw-2d.md` | `IVDW` 取值与编译前提、`LUSE_VDW`、二维库仑截断 |
| slab 建模、层数收敛、表面能、偶极校正、功函数 | `references/06-surfaces.md` | 切表面、两种表面能公式、`LDIPOL`/`IDIPOL`、平面平均 |
| 吸附能、ZPE/熵、CHE 台阶图、d 带中心、VASPsol | `references/07-catalysis.md` | 位点与覆盖度、气相参考、自由能校正、过电位、隐式溶剂两步流程 |
| NEB / CI-NEB / Dimer、能垒、唯一虚频验证 | `references/08-transition-state.md` | VTST 两阶段、`LCLIMB`、插值、`nebresults`、Dimer 精修 |
| 声子谱、热力学、DFPT、NAC、虚频诊断 | `references/09-phonon.md` | phonopy 有限位移与 DFPT、`band.conf`/`mesh.conf`、QHA/phono3py |
| 弹性常数、介电、光学、GW/BSE | `references/10-elastic-optics.md` | `IBRION=6`、Voigt 重排与单位、Born 判据、VRH、`LOPTICS`/`LEPSILON`/`LCALCEPS` |
| AIMD、NVT/NpT、RDF/MSD、MLFF | `references/11-md-mlff.md` | 系综选择、盒子塌缩对策、轨迹分析、train→refit→run 与日志 |
| Bader、差分电荷、ELF、部分电荷/STM、COHP | `references/12-charge-analysis.md` | 全电子参考链路、三次一致性、`LPARD`、spilling 判据 |
| Wannier90、有效质量、反折叠、费米面、缺陷 | `references/13-advanced-electronic.md` | `win` 关键词、单位换算、easyunfold、IFermi、电荷修正 |
| 后处理工具选型、批处理、**报错速查** | `references/14-postprocess-errors.md` | VASPKIT/pymatgen/ASE/py4vasp 分工、Slurm 数组、报错大表、合规红线 |

## 2. 通用骨架（从这里改，别从零写）

```
SYSTEM = <体系名>
PREC   = Accurate        # 应力/声子/弹性必须 Accurate；粗测可用 Normal
ENCUT  = 450             # ≥ 1.3 × max(ENMAX)；变胞与应力计算尤其
EDIFF  = 1E-6            # 生产 1E-6；声子/弹性 1E-7~1E-8；杂化 1E-5 起步
NELM   = 100
ISMEAR = 0               # 半导体/分子 0；金属 1；DOS 与静态 -5（见 §3）
SIGMA  = 0.05            # 金属 0.10；分子 0.01
LREAL  = .FALSE.         # 小体系关闭；>100 原子可 Auto
LASPH  = .TRUE.
LWAVE  = .TRUE.          # 下一步要读 WAVECAR 时必需
LCHARG = .TRUE.          # 做 Bader/差分电荷时必需
```

KPOINTS（Γ 中心最安全，避开 fcc/六方偶数 MP 网格的 `IBZKPT` 报错）：

```
Automatic generation
0
Gamma
  9  9  9
  0  0  0
```

三条不可省的核对：
- **`ISMEAR` 必须显式写**。VASP 默认 `ISMEAR=1, SIGMA=0.2` 是给金属的，半导体忘了写会算错占据、能量、声子频率。
- **POTCAR 顺序 = POSCAR 第 6 行顺序 = 第 7 行原子数分组**。第 6 行只是注释，写错不报错但算错：写完 `grep TITEL POTCAR` 对 `sed -n 6p POSCAR`。
- **标签拼错会被静默忽略**（`ISMAER` → 实际用默认 1）；Tab 与 Windows CRLF 是 INCAR 头号杀手，`dos2unix` + `cat -A` 检查。

## 3. 三个关键参数怎么定

- **ENCUT**：`≥ 1.3 × max(ENMAX)`（`grep ENMAX POTCAR`）。弹性/声子/应力再 +30%；要出数字就给 15% 步长的收敛曲线。变胞优化时截断不足会因 Pulay 应力把晶胞系统性算小。
- **k 点**：金属/导体取密（`KSPACING` 0.15–0.2 Å⁻¹ 量级），半导体可疏；**六角/三方/四方一律 Γ 中心**；分子、大真空盒子用 `Gamma 1 1 1` + `vasp_gam`。
- **`ISMEAR`**：金属 `1` 且每原子 `|T*S| < 1 meV`（不够就减小 `SIGMA` 并加密 k 点）；半导体/绝缘体/分子 `0`；DOS 与静态用 `-5`（**要求 Γ 中心且不可约 k 点 ≥ 4**，否则 `Tetrahedron method fails`）；`-5` 的力与应力在金属里偏差 5–10%，不能用于金属结构优化。

## 4. 十条硬约束（违反会产生错误结果或直接报错）

1. **杂化泛函与 meta-GGA 不能用 `ICHARG = 11` 做能带**（哈密顿量还依赖规则网格上的轨道/动能密度）→ 零权重 k 点法或 `KPOINTS_OPT`（≥ 6.3）。
2. **`ISMEAR = -5` 要 Γ 中心且 ≥ 4 个不可约 k 点**；slab、分子、Γ 点计算别用。
3. **slab 不用 `ISIF = 3`**（真空方向无应力支撑，晶胞会塌缩）→ 一律 `ISIF = 2`，面内晶格取自同设置的体相优化。
4. **相互比较的能量必须同设置**：POTCAR / ENCUT / PREC / LASPH / LREAL / 泛函 / vdW 方法 / U 值全同，否则误差直接进能量差——"吸附能为正且很大"的头号原因。
5. **位移超胞做声子绝不能弛豫**：`IBRION = -1`、`NSW = 0`，否则位移被抹平，结果是整体大虚频。
6. **DFPT（`IBRION = 7/8`）不支持 `IVDW`**（要 vdW 就改有限位移），且不能出现 `NCORE`/`NPAR`（`PEAD routines do not work for NCORE`）。
7. **SOC/非共线用 `vasp_ncl`**，INCAR 里**不写 `ISPIN`**（6.5.0 中 `ISPIN=2` + `LSORBIT` 会因 MAGMOM 个数检查报错）；MAE 必须 `ISYM = -1`。
8. **`ICHARG = 11` 必须与自洽步同 ENCUT/PREC/POTCAR/ISPIN**，DFT+U 还要 `LMAXMIX = 4`（f：6），否则 `dimensions on CHGCAR file are different`。
9. **`IBRION = 0` 时 `POTIM` 没有默认值**，必须显式给；NEB 的总 MPI 进程数必须是 `IMAGES` 的整数倍。
10. **`LATTICE_CONSTRAINTS` 配 `IBRION = 1/2` 在 < 6.4.3 被静默忽略**（c 轴照样动、真空层被压），推荐 6.5.1+。

## 5. 算完的通用验收（别只看 `reached required accuracy`）

- 电子步收敛：`grep "aborting loop because EDIFF is reached" OUTCAR`；确认没撞 `NELM`（撞了就是没收敛）。
- 总能：金属取 `energy(sigma->0)`（OUTCAR 的 `energy  without entropy` 段），不要用 `free  energy   TOTEN`。
- 离子步：`grep "reached required accuracy" OUTCAR`，再自己看最后一步最大力。
- 变胞/应力：`grep "external pressure" OUTCAR`，|P| 在 1–2 kB 内才算真弛豫。
- 参数是否真被接受：`grep -i "warn\|not found\|unknown" OUTCAR`。
- 目录体检：`python $VASP_TUTORIALS_ROOT/scripts\check_vasp.py`；已知可忽略的输出 vs 真错误见 `14-postprocess-errors.md`。

## 6. 相关资源

- 讲义仓库（可运行输入）：`$VASP_TUTORIALS_ROOT/`
- 服务器提交与排队：`remote-dft-server` skill（`qe_*` 工具族）
- Quantum ESPRESSO: see the official Quantum ESPRESSO documentation.

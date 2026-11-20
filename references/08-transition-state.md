# 过渡态：NEB / CI-NEB / Dimer / 虚频验证

> 用户要反应路径与能垒、扩散/解离势垒，要给 NEB 或 CI-NEB 输入、读 nebresults 结果、Dimer 精修鞍点，或用频率证明"这到底是不是过渡态"时读这个文件。

## 1. 先判断走哪条路线

| 情形 | 做法 | 关键区别 |
|---|---|---|
| 初态、末态都已知 | VTST CI-NEB 两阶段：普通 NEB 粗收敛 → `LCLIMB=.TRUE.` 精修 | 给出整条 MEP；能垒取自爬坡图像/样条极值 |
| 已有 NEB，要更严的鞍点 | `neb2dim.pl` 起 Dimer（`ICHAIN=2`） | 单构型，每步 2 次力计算；鞍点精度受 `EDIFF` 限制 |
| 只知道初态/方向 | Dimer（`MODECAR` 给初始方向） | 可能收敛到无关鞍点，必须用频率验阶 |
| 手上没有 VTST 版 VASP | 内置 NEB（`IBRION=3`+`SMASS=2`）或多加密 `IMAGES`，或用 ASE NEB | 无爬坡图像；内置 NEB 不能处理变胞 |
| 要确认鞍点阶数 / 算指前因子 | `IBRION=5` 频率（或 VTST 动力学矩阵） | 一阶鞍点 = 有且仅有一个虚频 |

## 2. 标准流程

**步骤 1 端点优化**

```bash
PREC=Normal   ENCUT=350    EDIFF=1E-6      # Pt 的 ENMAX≈230，350 够
ISMEAR=1      SIGMA=0.20   LREAL=Auto
IBRION=2      ISIF=2       NSW=200   EDIFFG=-0.02
```
检查：
- `grep "reached required accuracy" OUTCAR`，两端都要有。
- `dist.pl IS/CONTCAR FS/CONTCAR`：IS→FS 总位移，用来定图像数。
- `python idpp_interpolate.py IS/CONTCAR FS/CONTCAR --check-only`：查原子数、元素顺序、晶胞、Selective dynamics，并列出位移最大的原子。位移最大的不止一个且 >3 Å = 原子顺序不一致（跨周期边界的 0.98 Å 应看作 0.02 Å）。
- 端点 `EDIFFG` 要比 NEB 更严。

**步骤 2 插值**

```bash
bash make_images.sh IS FS 3           # IDPP 默认；自动拷 IS/OUTCAR -> 00/、FS/OUTCAR -> 04/
bash make_images.sh IS FS 3 linear    # 只想看线性路径时用
nebmake.pl IS/CONTCAR FS/CONTCAR 3    # 无 VTST 时的等价命令
```
图像数按相邻间距 ≲0.8–1 Å 定（单原子跳 2.81 Å → 3 图像，间距 ~0.7 Å）；对称路径用奇数；中间有稳定极小就拆成两段。检查 `00/POSCAR` 元素行与 IS 相同、`00/OUTCAR` 存在。

**步骤 3 stage1：普通 NEB 预收敛**

```bash
IMAGES=3   SPRING=-5   ICHAIN=0   LCLIMB=.FALSE.
IBRION=3   POTIM=0     IOPT=3
ISIF=2     NSW=150     EDIFFG=-0.10
```
提交前必查：总 MPI 进程数必须是 `IMAGES` 的整数倍（3 图像 × 32 ranks = 96，job.slurm 里含 `NP % IMAGES` 检查）。`NCORE`/`KPAR` 在每个图像内部生效。

**步骤 4 监控并切到 stage2（CI-NEB）**

```bash
nebef.pl          # 图像号 | 最大受力 | 能量 | 相对 00 的能量
vfin.pl stage1    # 归档 stage1，并把各图像 CONTCAR 覆盖回 POSCAR
```
换 stage2 的 INCAR 关键行：
```bash
LCLIMB=.TRUE.   IOPT=1   NSW=300   EDIFFG=-0.03
```

**步骤 5 出能垒**

```bash
nebresults.pl     # 生成 neb.dat spline.dat exts.dat mep.eps movie.xyz
python neb_tools.py collect stage2 --out summary.dat   # 仓库辅助脚本
```
`exts.dat` 的极值即 `E_TS`；`Ea_f = E_TS - E_IS`，`Ea_r = E_TS - E_FS`，`ΔE = E_FS - E_IS`。ASE 等价：`NEBTools(images).get_barrier()`。

**步骤 6 Dimer 精修（可选）**

```bash
neb2dim.pl        # 由爬坡图像生成 saddle_dimer/（含 MODECAR）
```
```bash
ICHAIN=2   DdR=0.005   DRotMax=4   DFNMin=0.01   DFNMax=1.0
IBRION=3   POTIM=0     IOPT=2      EDIFF=1E-7    EDIFFG=-0.01
```
不加 `IMAGES`。判据见 `DIMCAR`：Force < |EDIFFG| 且 Curvature < 0。续算：`cp CENTCAR POSCAR; cp NEWMODECAR MODECAR`。

**步骤 7 频率验证**

```bash
python set_freq_selective.py 02/CONTCAR -i 37 -o POSCAR   # --zmin 13.0 亦可
```
```bash
EDIFF=1E-7   IBRION=5   POTIM=0.015   NFREE=2   NSW=1
```
只放开与反应相关的原子，且与 IS 放开的原子**完全相同**。指前因子 `ν* = ∏ν_IS / ∏_{j≠imag}ν_TS`，典型 1e12–1e13 s⁻¹。

## 3. 关键标签表

| 标签 | 取值 | 何时用 / 注意 |
|---|---|---|
| `IMAGES` | 3–7（相邻 ≲1 Å） | 只在 NEB；总 ranks 必须是它的整数倍 |
| `SPRING` | −5 | 所有教程 NEB 输入都用 −5，全程别改 |
| `ICHAIN` | 0 = NEB；2 = Dimer | 用 `ICHAIN≥1` 就必须配 IOPT（需 VTST 版） |
| `LCLIMB` | stage1 `.FALSE.` → stage2 `.TRUE.` | VTST 版默认就是 `.TRUE.`，粗收敛必须显式关掉 |
| `IBRION`/`POTIM` | 3 / 0（NEB、Dimer） | 不能用 `IBRION=2`：CG 线搜索依赖能量 |
| `IOPT` | 1=LBFGS 2=CG 3=QM 4=SD 7=FIRE 0=VASP 优化器 | 3 最稳健；Dimer 用 2 |
| `EDIFF` | NEB 1E-6；Dimer、频率 1E-7 | 力噪声会让爬坡图像抖 |
| `EDIFFG` | stage1 −0.10；stage2 −0.03；Dimer −0.01 | 负值 = 力判据 |
| `MAXMOVE`/`TIMESTEP`/`ILBFGSMEM` | 默认 0.2 Å / 0.1 / 20 | 能量振荡时减小 `MAXMOVE`/`TIMESTEP` |
| `DdR`/`DRotMax`/`DFNMin`/`DFNMax` | 0.005 / 4 / 0.01 / 1.0 | 仅 Dimer |
| `ISIF` | 2（固定晶胞） | 要变胞只能用 SS-NEB：`LNEBCELL=.TRUE.` + `ISIF=3` + `IOPT=3` |
| `SMASS` | 2 | 仅内置 NEB（配合 `IBRION=3` 阻尼动力学） |

## 4. 易错点

1. **现象** 报 `M_divide: can not subdivide ...` → **原因** 总进程数不是 `IMAGES` 的整数倍 → **处置** 改成 96 = 3×32 之类，或用 `nebef` 看实际起了几个图像后重排。
2. **现象** 图像乱飞、最大受力几百 eV/Å → **原因** 原子顺序不一致、跨周期边界算成 0.98 Å、只做了线性插值 → **处置** `idpp_interpolate.py --check-only` 查顺序后重新 IDPP 插值。
3. **现象** 设了 `IOPT` 但原子完全不动 → **原因** 用的不是 VTST 版，`POTIM=0` 时普通 VASP 不移动原子 → **处置** 换 VTST 编译，或改内置 NEB。
4. **现象** 以为在跑普通 NEB，最高图像却一直在爬 → **原因** VTST 版 `LCLIMB` 默认 `.TRUE.` → **处置** stage1 显式写 `LCLIMB=.FALSE.`。
5. **现象** 用 `IBRION=2` 跑 NEB/Dimer，报 `ZBRENT: fatal error in bracketing` → **原因** CG 线搜索要求能量对位置连续 → **处置** 一律 `IBRION=3` + `POTIM=0`，优化器交给 IOPT。
6. **现象** `nebresults.pl` 端点能量为 0 或找不到 00 的能量 → **原因** 没把端点 OUTCAR 拷进 `00/` 与 `0N/` → **处置** 补拷（`make_images.sh` 会自动拷，`nebmake.pl` 不会）后重跑分析。
7. **现象** 能垒曲线出现两个极大 → **原因** 路径中间存在稳定中间体 → **处置** 在极小处切成两段 NEB 分别算。
8. **现象** Dimer 提示需要 `IOPT`；或曲率一直为正、Torque 不降 → **原因** 非 VTST 版、初始 `MODECAR` 方向不对、力太噪 → **处置** `EDIFF=1e-8`、`DdR=0.01`、增大 `DRotMax`、换初始方向。
9. **现象** TS 频率出现 2 个以上虚频 → **原因** 高阶鞍点、结构还没到鞍点，或 TS 放开的原子比 IS 多 → **处置** 回到 CI-NEB/Dimer 用更严判据重精修，让 IS 与 TS 放开同一批原子；只要一个虚频，本征矢应沿反应方向。
10. **现象** CI-NEB 能量来回振荡不收敛 → **原因** 一开始就开爬坡、LBFGS 步子太大、SCF 不够严 → **处置** 先 stage1 粗收敛，减小 `TIMESTEP`/`MAXMOVE`，`EDIFF=1E-6`。

## 5. 算完了怎么验收

- `neb.dat`/`spline.dat` 横坐标是累计路径长度，`exts.dat` 给出每个极值（能量 + 位置）；最高点须在中间图像附近，落在端点上说明图像太少或端点没优化好。
- `nebef.pl` 最后一步的最大受力：收敛时应小于 |EDIFFG|（stage2 为 0.03 eV/Å）。没达标 = 还没收敛，能垒不可引用。
- `DIMCAR`：最后几行 Force < |EDIFFG| 且 Curvature 为负；表里的 Energy 列**不是**鞍点能量，鞍点能量取该 OUTCAR 最后一个 `energy(sigma->0)`。
- 虚频：OUTCAR 里 `f/i=` 只出现一次，且本征矢沿反应坐标（对照 IS→TS 的 CONTCAR 位移方向）。无虚频 = 停在极小；<~50 cm⁻¹ 的虚频是数值噪声。
- 一致性：`Ea_f - Ea_r` 应约等于 `E_FS - E_IS`；差得多说明端点用了不同电子参数。
- 指前因子 `ν*` 落在 1e12–1e13 s⁻¹ 量级；差一个数量级以上多半是虚频或模数取错。
- 失败信号：能垒 <0.05 eV 或为负、图像 CONTCAR 里原子跑进真空层、`movie.xyz` 里原子集体平移 → 插值或原子顺序错了，重做步骤 1–2。

## 6. 本地完整输入

$VASP_TUTORIALS_ROOT/examples/11-pt-adatom-neb/
$VASP_TUTORIALS_ROOT/examples/11-li-si-neb/
$VASP_TUTORIALS_ROOT/examples/11-dimer-from-neb/
$VASP_TUTORIALS_ROOT/examples/11-ts-freq-check/
$VASP_TUTORIALS_ROOT/tutorials/11-neb-dimer.md
$VASP_TUTORIALS_ROOT/scripts/neb_tools.py

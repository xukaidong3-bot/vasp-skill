# 能带与 DOS/PDOS：两步法、KPOINTS_OPT 与超越 GGA 的带隙

> 什么时候读这个文件：用户要能带图、带隙数值、DOS/PDOS、fat band，或嫌 PBE 带隙太小要上 HSE06 / r2SCAN / mBJ 时。

## 1. 先判断走哪条路线

| 情形 | 做法 | 关键区别 |
|---|---|---|
| PBE/LDA/PBE+U 能带 | 自洽 → `ICHARG = 11` 两步法 | 全版本通用；两步参数须一致 |
| 同上，VASP ≥ 6.3 | `KPOINTS` 网格 + `KPOINTS_OPT` 路径 | 单次作业；必须 `vasp_std`，不能用 `vasp_gam` |
| DOS / PDOS | `ICHARG = 11` + 密网格 + `ISMEAR = -5` + `LORBIT = 11` | 网格比 SCF 密 1.5–2 倍；`-5` 不能配 Line-mode |
| HSE06 / r2SCAN 能带 | 零权重 k 点 或 `KPOINTS_OPT` | 禁 `ICHARG = 11`：Fock 算符与 τ 依赖规则网格占据态 |
| 只要一个可信带隙 | `METAGGA = MBJ` + `KPOINTS_OPT` | 只有交换势、无能量泛函：不能优化、不能比总能 |

## 2. 标准流程

**1) 自洽（所有路线起点）**

```text
PREC    = Accurate
ENCUT   = 400        # 每步都相同
EDIFF   = 1E-7
ISMEAR  = 0
SIGMA   = 0.01
LREAL   = .FALSE.
NSW     = 0
IBRION  = -1
LCHARG  = .TRUE.     # CHGCAR 供 ICHARG=11
LWAVE   = .FALSE.
# BANDGAP = KPOINT   # 带隙报告详略，>= 6.4.3
```

KPOINTS `Gamma` `8 8 8`。检查：`CHGCAR` 非空；`grep "E-fermi" OUTCAR | tail -1` 记下 E_F；`grep -A12 "Band structure" OUTCAR` 读网格带隙。

**2) 能带（两步法）** Line-mode KPOINTS（Si 路径 `L-G-X-U|K-G`，每段 40 点含端点，`Line-mode` + `Reciprocal`）：

```text
Si band path L-G-X-U|K-G
40
Line-mode
Reciprocal
  0.500 0.500 0.500  L
  0.000 0.000 0.000  GAMMA

  0.000 0.000 0.000  GAMMA
  0.500 0.000 0.500  X

  0.500 0.000 0.500  X
  0.625 0.250 0.625  U

  0.375 0.375 0.750  K
  0.000 0.000 0.000  GAMMA
```

INCAR 复制上一步，只改：

```text
ICHARG = 11
NBANDS = 16     # 占据带 + 空带；Si 原胞 8 电子
LORBIT = 11     # fat band
LCHARG = .FALSE.
```

检查 `grep "NKPTS" OUTCAR | head -1` 应为 160（4×40）；一致性判据见 §4 第 2 条。

**3) DOS/PDOS** KPOINTS 换 `Gamma` `16 16 16`，INCAR 改：

```text
ISMEAR = -5     # 四面体法：要求 Γ 中心且不可约 k >= 4
NEDOS  = 3001
EMIN   = -15    # 绝对能量，不是相对 E_F
EMAX   = 10
```

**4) 一步法（≥ 6.3）** `KPOINTS` 仍 8×8×8，另建 `KPOINTS_OPT`（同一路径）；INCAR 同自洽 + `NBANDS`、`LORBIT = 11`，`vasp_std` 跑。能带读 `vaspout.h5`（`kpoints_opt`）或 `vasprun.xml` 的 `<eigenvalues_kpoints_opt>`，**不读 EIGENVAL**；fat band 取 `PROCAR_OPT`。

**5) HSE06 自洽** 先跑 PBE（同网格、`NBANDS = 16`、`LWAVE = .TRUE.`）拿 WAVECAR，再：

```text
ISTART   = 1
LHFCALC  = .TRUE.
GGA      = PE
HFSCREEN = 0.2      # 0.3 → HSE03；不写 → PBE0
AEXX     = 0.25
ALGO     = Damped   # 绝不用 Fast：不被正确支持且不报警告
TIME     = 0.4      # 振荡时降 0.2–0.3
PRECFOCK = Normal
HFRCUT   = -1       # 有带隙体系与能带必须
# NKRED  = 2        # 金属/掺杂禁用，每个方向须被 2 整除
```

KPOINTS 与 PBE 相同的 `6 6 6`。

**6) HSE 能带** 二选一：
- 零权重：`make_zero_weight_kpoints.py --ibzkpt ../0-pbe-scf/IBZKPT --kpath KPATH --out KPOINTS`；先 `INCAR.pbe`（同一 KPOINTS）得 WAVECAR，再 `INCAR.hse`。等价 `vaspkit -task 251`、`sumo-kgen --hybrid`；作图须剔除非零权重点（`vaspkit -task 252`、`sumo-bandplot`）。
- `KPOINTS_OPT`：同第 4 步，且 `HFRCUT = -1` 必写。

**7) r2SCAN / mBJ**

```text
METAGGA = R2SCAN    # >= 6.2.0
LASPH   = .TRUE.
LMIXTAU = .TRUE.    # 让 τ 参与混合，收敛更稳
ALGO    = Normal    # 困难时 All
```

```text
METAGGA = MBJ
LASPH   = .TRUE.
ALGO    = All
TIME    = 0.1
NELM    = 200
# CMBJ  = 1.20      # slab/界面/续算必须固定：grep CMBJ OUTCAR | tail -1
```

两者都先 `grep kinetic POTCAR` 确认势含芯电子动能密度。

## 3. 关键标签表

| 标签 | 取值 | 何时用 / 注意 |
|---|---|---|
| `ICHARG` | `11` | 仅半局域泛函可用 |
| `NBANDS` | 占据带 + 若干空带（Si：16） | 太少 → 导带与 DOS 高能处被截断 |
| `LORBIT` | `11` | PAW 投影，无需 `RWIGS`；给 PROCAR 与 lm 分解 DOSCAR |
| `NEDOS`/`EMIN`/`EMAX` | `3001` / `-15` / `10` | 默认 301 峰太粗；EMIN/EMAX 是绝对能量，先看 E_F 再定 |
| `ISMEAR` | DOS 用 `-5`，能带步用 `0` | `-5` 要求 Γ 中心、不可约 k ≥ 4、非 Line-mode |
| `EFERMI` | `MIDGAP` | 半导体把 E_F 放带隙中间（≥ 6.4） |
| `LKPOINTS_OPT`/`KPOINTS_OPT_NKBATCH`/`KSPACING_OPT` | `.FALSE.` / 10–20 / 0.5 | 分别：禁用、控内存（≥6.3）、自动更密网格（≥6.6） |
| `HFSCREEN`/`AEXX` | `0.2` / `0.25` | 把 `AEXX` 调到 0.3–0.4 拟合宽带隙氧化物时，泛函已不是 HSE06 |

## 4. 易错点

1. **杂化/meta-GGA 用了 `ICHARG = 11`** → 结果错或跑不出能带。原因：Fock 算符要规则网格上的占据轨道、meta-GGA 要 τ，只读 CHGCAR 重建不出。处置：零权重 k 点或 `KPOINTS_OPT`。
2. **两步法参数不一致** → 路径上 Γ 点本征值与 SCF 的 Γ 差几十 meV。原因：`ENCUT`/`PREC`/`GGA`/`METAGGA`/`LASPH`/vdW/`ISPIN`/`MAGMOM`/`LDAU*`/`LMAXMIX` 有一处不同，或 CHGCAR 不属于该结构。处置：复制 SCF 的 INCAR，只改 `ICHARG`/`KPOINTS`/`NBANDS`/`LORBIT`。
3. **拿路径算出的 E-fermi 作图** → 费米面附近整条带平移（`ICHARG=11` 步的 E_F 由路径 k 点算出）。处置：用 SCF 或密网格 DOS 的 E_F，或 `EFERMI = <数值>` 固定。
4. **`ISMEAR = -5` 配 Line-mode 或小网格** → `Tetrahedron method fails (number of k-points < 4)` / `Routine TETIRR needs special values`。处置：DOS 用 Γ 中心密网格，能带步改 `ISMEAR = 0`。
5. **六角/三方体系用偶数 MP 网格** → `VERY BAD NEWS! internal error in subroutine IBZKPT: Reciprocal lattice and k-lattice belong to different class of lattices`。处置：一律 Γ 中心。
6. **高对称点坐标与晶胞约定不匹配** → 能带"看着怪"、带隙位置不对。处置：用生成路径工具输出的原胞重做 SCF（pymatgen `kpath.prim`、VASPKIT `PRIMCELL.vasp`），或改用 ASE `cell.bandpath()`。
7. **HSE 的两个隐性坑** → `ALGO = Fast` 不被正确支持且不给警告，改 `Damped`（`TIME = 0.4`，振荡降 0.2–0.3）/`All`/`Normal`；路径 k 点上默认 `HFRCUT = 0`（辅助函数法）会让 Γ 附近出现台阶，改 `HFRCUT = -1` 重算。
8. **零权重法从已收敛的 HSE WAVECAR 续算** → 几步就"收敛"但能带乱（判据只看总能）。处置：从同 k 集的 PBE WAVECAR 起步，或加大 `NELMIN`。
9. **`KPOINTS_OPT` 的版本与参数坑** → `vasp_gam` 跑会报错（6.4.3 起明确提示），换 `vasp_std`；`NKRED = 2` 后带隙变化 > 0.1 eV 说明降采样不适用（金属/部分占据），改回 1。
10. **meta-GGA/mBJ 的两个坑** → 第一步发散多为 POTCAR 缺芯电子动能密度或 τ 未参与混合：`grep kinetic POTCAR` 换势、`LMIXTAU = .TRUE.`、`ISTART = 1` 从 PBE WAVECAR 起步；mBJ 没有对应能量泛函（VASP 用 LDA 表达式占位），不能优化结构、不能比总能，结构用 PBE/r2SCAN 优化、mBJ 只取带隙。
11. **`PROCAR_OPT` 损坏 / 自旋极化 PDOS 缺下自旋通道** → 6.5.0–6.5.1 的 issue 45、6.6.0 的 issue 69。处置：用 6.6.1。

## 5. 算完了怎么验收

- **文件**：能带目录有 `EIGENVAL`+`PROCAR`（两步法）或 `PROCAR_OPT`+`vaspout.h5`（一步法）；DOS 目录有 `DOSCAR`。
- **两步法自检**：`grep -A7 "k-point     1 :" ../1-scf/OUTCAR | tail -6` 与 `grep -A7 "k-point    40 :" OUTCAR | tail -6` 本征值差应 < 1 meV。
- **带隙量级**：Si 原胞 PBE 约 0.5–0.6 eV、HSE06 约 1.1–1.3 eV、mBJ 约 1.1–1.2 eV、实验 1.17 eV。PBE 低估 40–60% 是系统误差、不是算错。Si 的 CBM 在 Γ–X 线约 85% 处，**不在** 8×8×8 网格点上——网格报告的带隙会偏大。SCF 网格 8×8×8 → 12×12×12 后路径能带变化应 < 10 meV，DOS 网格加密一级主峰位移 < 0.05 eV。
- **DOS**：`sed -n 6p DOSCAR` 的 `EMAX EMIN NEDOS E-fermi` 与 INCAR 一致；E_F 处积分 DOS 应等于 `NELECT`（Si 原胞 8）。PDOS 之和小于 TDOS 属正常（PAW 球间隙区未投影）。
- **判定算错的信号**：带隙为 0 或负；某条带半路断掉（`NBANDS` 不足）；DOS 呈锯齿（`NEDOS` 太小）；带隙随路径形状剧变（路径与晶胞不匹配）。

## 6. 本地完整输入

- `$VASP_TUTORIALS_ROOT/tutorials/05-band-dos.md`
- `$VASP_TUTORIALS_ROOT/tutorials/06-hybrid-metagga.md`
- `$VASP_TUTORIALS_ROOT/examples/05-band-dos-si/`（`1-scf`/`2-band-icharg11`/`3-dos`/`4-band-kpoints-opt`）
- `$VASP_TUTORIALS_ROOT/examples/06-hse06-band-si/`（`0-pbe-scf`/`1-hse-scf`/`2-hse-band-zero-weight`/`3-hse-band-kpoints-opt`）
- `$VASP_TUTORIALS_ROOT/examples/06-r2scan-relax-zno/`
- `$VASP_TUTORIALS_ROOT/examples/06-mbj-gap-si/`

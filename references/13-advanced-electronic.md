# 高级电子结构：Wannier90、有效质量、反折叠、费米面、拓扑与缺陷

> 什么时候读这个文件：要 Wannier 插值能带/紧束缚模型、有效质量 m*、超胞能带反折叠、金属费米面、Berry 曲率/陈数/Z2/Rashba，或带电点缺陷形成能时。普通 PBE 能带与 DOS 不在本文件。

## 1. 先判断走哪条路线

| 情形 | 做法 | 关键区别 |
|---|---|---|
| 插值能带、`wannier90_hr.dat`（喂 postw90/WannierTools） | `-DVASP2WANNIER90` 接口 + 独立 `wannier90.x` | 需 `ISYM = -1` 完整均匀网格 |
| 只要带边曲率 → m* | 带边 ±Δk 的 `ICHARG = 11` 非自洽，或 `KPOINTS_OPT`（≥6.3） | 不碰 Wannier；对 EDIFF、步长、泛函极敏感 |
| 超胞（缺陷/合金/表面）能带 | easyunfold（或 vaspkit 281–285） | 必须 `LWAVE = .TRUE.`；k 点用**原胞**路径 |
| 金属费米面、费米速度 | 稠密均匀非自洽网格 + IFermi | 必须用原胞（超胞会折叠）；Cu 用 25×25×25 |
| Berry 曲率、陈数、Z2、Rashba | `vasp_ncl` + `LSORBIT`；Wannier 化后交 postw90/WannierTools/Z2Pack | 需 SOC + 完整网格 |

带电缺陷也不需新流程：`NELECT` 扫电荷态 + FNV/eFNV（§2 第 9 步）。

## 2. 标准流程

1. **编译**：VASP 加 `CPP_OPTIONS += -DVASP2WANNIER90`、`LLIBS += -L$(WANNIER90_ROOT)/lib -lwannier`；Wannier90 3.1.0 必须**串行**编库（注释 `make.inc` 的 `COMMS = mpi`）再 `make lib`。出现 `wannier90.amn/.mmn/.eig` 才算接口编好。
2. **自洽**：`PREC = Accurate`、`ENCUT = 400`、`EDIFF = 1E-8`、`ISMEAR = 0`、`SIGMA = 0.05`、`LCHARG = .TRUE.`，Γ 中心 8×8×8；记下 `grep E-fermi OUTCAR | tail -1`（窗口同零点）。
3. **非自洽 + 写接口文件**（KPOINTS 换 6×6×6 完整网格）：

```fortran
ICHARG  = 11         ! 读自洽 CHGCAR；ENCUT/PREC 必须与第 2 步一致
ISYM    = -1         ! 完整均匀网格：6×6×6 = 216 个 k 点
NBANDS  = 16         ! 窗口内每 k 点态数 ≥ NUM_WANN
NUM_WANN = 8         ! Si 两原子各 4 个 sp3 函数
LWANNIER90 = .TRUE.  ! 写 wannier90.win/.amn/.mmn/.eig
LCHARG = .FALSE.
LWAVE  = .FALSE.
```

`LWANNIER90_RUN = .TRUE.` 走库模式一次算完，但默认不写 `.amn/.mmn`（保留需 `LWRITE_MMN_AMN = .TRUE.`）。
4. **`wannier90.win` 窗口**：绝对能量，与 VASP 本征值同零点；Si 例 `dis_froz_max = E_F+2`、`dis_win_max = E_F+11`（E_F≈6 eV 时即 8.0/17.0）。

```text
num_wann = 8
num_iter = 200
dis_num_iter = 500
dis_froz_max = 8.0
dis_win_max  = 17.0
begin projections
Si:sp3
end projections
guiding_centres = true
write_hr = true
bands_plot = true
begin kpoint_path
L 0.5 0.5 0.5  G 0 0 0
G 0 0 0  X 0.5 0 0.5
end kpoint_path
```

冻结上限取在「目标子空间完整、任一 k 点窗口内态数 ≤ `num_wann`」处；别把 win 塞进 `WANNIER90_WIN`（8192 字符上限，issue 64）。
5. **独立运行查质量**：`wannier90.x wannier90`。`grep "<-- DIS" wannier90.wout | tail -3`（Ω_I 不再变）、`grep "<-- CONV"`（spread 不再降）、`grep -A 12 "Final State"`（8 个中心在 Si–Si 键上）。
6. **有效质量**：`KPOINTS_fd` 用 Cartesian、单位 2π/a，Γ + [100]/[110]/[111] 各 5 点、步长 0.005·2π/a（a=5.653 Å 时 ≈0.0056 Å⁻¹）；INCAR：`EDIFF = 1E-8`、`SIGMA = 0.01`、`ICHARG = 11`、`NBANDS = 16`、`PREC = Accurate`。

```bash
python3 fit_effmass.py --eigenval EIGENVAL --poscar POSCAR  # 默认 CBM+VBM
python3 fit_effmass.py --npoints 3  # 内侧 3 点查非抛物性
```

拟合 `E = E0 + a·Δk²`，**m\*/m0 = 3.80998/a**（a 单位 eV·Å²，k 单位 Å⁻¹ 且**含 2π**）。HSE 不能 `ICHARG = 11`：复制 `IBZKPT` 为 KPOINTS、路径点权重给 0；或（≥6.3）用 `KPOINTS_OPT`（PBE0 加 `HFRCUT = -1`）读 `EIGENVAL_OPT`。
7. **超胞反折叠**：

```bash
sumo-kgen -p POSCAR_prim  # 原胞路径 -> KPOINTS_band
easyunfold generate POSCAR_prim POSCAR_sc KPOINTS_band --matrix "3 3 3"
# 超胞非自洽：ICHARG = 11  NBANDS = 144  LORBIT = 11  LWAVE = .TRUE.
easyunfold unfold calculate WAVECAR
easyunfold unfold plot --eref <E-fermi> --emin -4 --emax 4
```

`--matrix` 是行向量约定 `cell_super = M @ cell_primitive`。`vasp_gam`/`vasp_ncl` 的 WAVECAR 分别加 `--gamma`/`--ncl`；杂化泛函用 `--scf-kpoints IBZKPT --nk-per-split 60`。
8. **费米面**：`ICHARG = 11`、`ISMEAR = 1`、`SIGMA = 0.1`、`NBANDS = 12`（Cu 第 6 条带穿 E_F）、`LORBIT = 11`、`NEDOS = 3001`，Γ 中心 25×25×25。再 `ifermi info`、`ifermi info --property velocity`、`ifermi plot --interpolation-factor 10`、`ifermi plot --slice 1 1 0 0`。IFermi 读 vasprun.xml 插值；`--property spin` 需非共线 + `ISYM = -1` + 全 BZ。
9. **拓扑/自旋与带电缺陷**：`LSORBIT = .TRUE.`（必须 `vasp_ncl`）、`ICHARG = 11`、`ISYM = -1`、`GGA_COMPAT = .FALSE.`。

Rashba：沿 −M→Γ→M 取密点，读交叉点与极值的能量差 E_R、动量差 k_R，`alpha_R = 2E_R/k_R`；自旋纹理用 vaspkit 651+652/653 或 pyprocar。Berry 曲率/陈数：WAVECAR 交 VASPBERRY，或 Wannier 化后 `postw90.x`（`berry = true`、`berry_task = ahc`）；Z2 用 WannierTools（`hr.dat` + `wt.in`）或 Z2Pack。缺陷：`NELECT` 比中性少 q 个电子即 q=+1，`E_f = E_tot[X^q] − E_tot[bulk] − Σn_iμ_i + q(E_F+E_VBM) + E_corr`，E_corr 用 FNV/eFNV（介电张量 + LOCPOT），工具 doped、pydefect。

## 3. 关键标签表

| 标签 | 取值 | 何时用 / 注意 |
|---|---|---|
| `LWANNIER90` / `LWANNIER90_RUN` | `.TRUE.` | 前者写 `.win/.amn/.mmn/.eig`；后者库模式算完（默认不写 `.amn/.mmn`） |
| `NUM_WANN` / `NBANDS` | 8 / 16 | 窗口内每 k 点态数 ≥ `NUM_WANN` |
| `WANNIER90_WIN` | 多行串 | INCAR 内嵌 win 有 8192 字符上限，长 win 用独立文件 |
| `EDIFF` / `SIGMA` | 1E-8 / 0.01 | m* 靠 meV 级差；金属改 `ISMEAR = 1`、`SIGMA = 0.1` |
| `LWAVE` | `.TRUE.` | 反折叠唯一必需项（WAVECAR 数 GB，算完删） |
| `LSORBIT` + `GGA_COMPAT` | `.TRUE.` / `.FALSE.` | SOC 用 `vasp_ncl` + `ISYM = -1` |
| `NELECT` | 默认值 − q | q = 中性电子数 − NELECT |
| win `dis_froz_max`/`dis_win_max` | E_F+2 / E_F+11 | 绝对能量（VASP 本征值同零点） |

## 4. 易错点

1. **设了 `LWANNIER90 = .TRUE.` 却无 `.amn/.mmn/.eig`** → 未用 `-DVASP2WANNIER90` 编译（链接报 `wannier_setup` = 库为 MPI 版）→ 3.1.0 串行 `make lib`。
2. **`dis_windows: More states in the frozen window than target WFs`** → 窗口太高或 `num_wann` 太小 → 降 `dis_froz_max` 或增 WF。
3. **k 点 / `mp_grid` 不一致** → 对称约化或网格非均匀 → `ISYM = -1` + 均匀 Γ 中心网格。
5. **m* 差一个数量级** → 漏 2π 或用分数坐标 → 用含 2π 的 Å⁻¹ 与 `m*/m0 = 3.80998/a`，确认 k0 为极值点。
6. **PBE 下 GaAs 电子质量远小于实验 0.067 m0** → 低估带隙使 m* 偏小 → 定量用 HSE06/mBJ；空穴质量须 SOC 下算。
7. **HSE 用 `ICHARG = 11`** → 不支持非自洽定密度 → 零权重 k 点法或 `KPOINTS_OPT`。
8. **反折叠几乎无权重** → `--matrix` 与超胞不符或 KPOINTS 未换 `KPOINTS_easyunfold` → `easyunfold unfold status` 核对重做；读 WAVECAR 失败加 `--gamma`/`--ncl`。
9. **`KPOINTS_WAN` 在立方晶体上明显错** → issue 105（等距晶格矢量只留最后一个）、issue 106（`LSCDM` 不经库时窗口用两次）→ 用 `bands_plot` 或与 `KPOINTS_OPT` 对比。
10. **IFermi 无费米面/面锯齿** → E_F 在带隙或网格太粗 → `--mu` 移到带边、加密网格、加大 `--interpolation-factor`。
11. **缺陷形成能随超胞剧变** → 未做镜像电荷与势对齐修正 → FNV/eFNV（doped、pydefect），≥10 Å 超胞验证。

## 5. 算完了怎么验收

- `wannier90.wout`：`<-- DIS` 的 Ω_I 与 `<-- CONV` 的 spread 不再变；`Final State` 的 8 个中心在 Si–Si 键上；插值能带在窗口内与 DFT 重合到 meV 级、无毛刺。
- 有效质量：`m*_fit` 与 `m*_2pt` 接近 ⇒ 步长合适；`--npoints 3` 与 5 点差别小 ⇒ 抛物线成立；正=电子、负=空穴。
- 反折叠：先用**完美** 3×3×3 超胞验证——权重只有 0/1、EBS 与原胞能带重合；空位超胞在带隙出现弱权重平带。
- 费米面：`ifermi info` 报 1 个穿 E_F 的面（Cu 第 6 条带）；近球面 + L 方向 neck；11³/17³/25³ 面积不再变 ⇒ 收敛。
- 缺陷：修正后形成能随超胞变化 < 约 0.1 eV；转变能级在带隙内且排序合理。算错信号：spread 仍在降、m* 量级异常、`KPOINTS_WAN` 与 `KPOINTS_OPT` 形状不一致却照用。

## 6. 本地完整输入

- `$VASP_TUTORIALS_ROOT/examples/17-wannier90-si/`
- `$VASP_TUTORIALS_ROOT/examples/17-effmass-gaas/`
- `$VASP_TUTORIALS_ROOT/examples/17-unfold-si-vacancy/`
- `$VASP_TUTORIALS_ROOT/examples/17-fermi-surface-cu/`
- `$VASP_TUTORIALS_ROOT/tutorials/17-advanced-electronic.md`

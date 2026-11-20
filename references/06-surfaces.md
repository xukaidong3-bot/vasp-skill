# Slab 建模、表面能、偶极校正与功函数

> 什么时候读这个文件（触发场景）：要切表面建 slab、扫层数/真空/k 点收敛、算表面能与表面弛豫、加偶极校正、算功函数或把带边对齐真空能级；吸附能见 07。

## 1. 先判断走哪条路线

| 情形 | 做法 | 关键区别 |
|---|---|---|
| 对称化学计量 slab（金属、洁净半导体面） | 一个 slab + 体相 ε，`γ = (E_slab(N) − N·ε_bulk)/(2A)` | 两侧等价、无净偶极，最省；2A 是上下两个面 |
| 只弛豫一侧（固定底层） | `γ_unrel = (E_unrel − Nε)/2A`，`γ_rel = (E_rel − Nε)/A − γ_unrel` | 必须扣掉未弛豫面的贡献 |
| 层数扫描线性漂移 / 极性面 | Fiorentini–Methfessel：拟合 `E_slab(N) = 2Aγ + Nε` | 一次拟合得 γ 与 ε，不依赖体相参考 |
| 氧化物不同终止（非化学计量） | `γ(μ)`，μ_O 由 `½E_O2 + ½ΔH_f ≤ μ_O ≤ ½E_O2` 限定 | 多终止面互相竞争，须画 γ–μ 图 |
| 不对称、带电、单侧吸附 | `LDIPOL = .TRUE.` + `IDIPOL = 3` + `DIPOL` | 否则有寄生电场，能量/功函数都错 |

## 2. 标准流程

**1) 先优化体相拿参考能量**（Cu 例 `INCAR_bulk` + Γ 中心 `16 16 16`）

```text
PREC = Accurate
ENCUT = 500
EDIFF = 1E-7
LREAL = .FALSE.
ISMEAR = 1
SIGMA = 0.10
IBRION = 2
ISIF = 3
NSW = 40
EDIFFG = -1E-3
```

查 `grep "reached required accuracy" OUTCAR`、`grep "F=" OSZICAR | tail -1`，`ε_bulk = E0/N_atom`。该能量须与 slab 同 ENCUT/smearing/`LREAL`。

**2) 切表面建 slab**（记录生成了哪些终止面）

```python
from ase.build import fcc111            # vacuum 是每侧厚度，自动居中
slab = fcc111("Pt", size=(2, 2, 4), a=3.97, vacuum=10.0)
# pymatgen: SlabGenerator(bulk, (1,1,1), min_slab_size=4, min_vacuum_size=7,
#   in_unit_planes=True, center_slab=True).get_slabs()
```

判据：`is_symmetric()` True、`is_polar()` False（Pt(111) 输出 `Pt16 | symmetric: True | polar: False`、`a, b, c = 5.614 5.614 25.213 A, area = 27.30 A^2`）。极性面（岩盐(111)、纤锌矿(0001)）属 Tasker III，能量随厚度发散，须重构/空位/吸附补偿。VASPKIT：803 切表面、801 加真空、804 异质结、402 按层固定。

**3) 固定底层**：`Selective dynamics` 行放在坐标行前，逐原子 `T T T` / `F F F`；4 层 slab 固定底部 2 层，对称 slab 固定中心层（`make_slabs.py` 默认如此）。ASE：`FixAtoms(mask=z < layers[1] + 0.1)`。固定层不超过一半，且固定的层要像体相。

**4) 弛豫 slab**

```text
PREC = Accurate
ENCUT = 400
EDIFF = 1E-6
NELM = 200
LREAL = .FALSE.
ISMEAR = 1
SIGMA = 0.10
IBRION = 2
ISIF = 2            # 一律 ISIF = 2；ISIF = 3 会让真空塌缩
NSW = 100
EDIFFG = -0.02
LDIPOL = .TRUE.     # 仅不对称/有净偶极时
IDIPOL = 3
DIPOL = 0.5 0.5 0.5
```

法向 k 点取 1、面内与体同密度（Pt(111) p(2×2) 用 `8 8 1`）。

**5) 扫收敛**：层数 3 5 7 9 11（金属 4–6 层起、表面能 5–11 层），真空 15/20/25 Å，k 点加密一档。判据：γ 变化 < 0.01–0.02 J/m²、Φ 变化 < 0.02–0.05 eV、总能 < 1 meV/atom；真空 ≥ 15 Å（函数 15–20 Å）。

**6) 算表面能**

```bash
python make_slabs.py bulk/CONTCAR --layers 3 5 7 9 11 --vacuum 15 --kin 16 --kz 12
bash run.sh bulk && bash run.sh slabs && bash run.sh collect
```

输出含 `oriented bulk: 2 atoms, eps_bulk = ...`、`linear fit over N = 7, 9, 11: eps_fit = ... (diff +0.0x meV)`、`Fiorentini-Methfessel gamma = 0.09xxx eV/A^2 = 1.4xxx J/m^2`（1 eV/Å² = 16.02 J/m²）。取向体相 `obulk/`（1×1、2 层、面内 16×16 法向 12 k）须与 slab 同 k 点/ENCUT/smearing。

**7) 偶极校正**：`LDIPOL = .TRUE.` 同时修势、受力、能量；`IDIPOL = 3` 沿第三晶格矢量；`DIPOL` 给 slab 质心分数坐标（居中 `0.5 0.5 0.5`；贴底 z ≈ 0.05–0.35 时约 0.2），用 `atoms.get_center_of_mass(scaled=True)` 取，台阶在与 `DIPOL` 相隔半晶胞处。难收敛加 `AMIN = 0.01`、`NELM = 200`、`EDIFF = 1E-6`；仍不收敛先关 `LDIPOL` 收敛再 `ISTART = 1` 续算。

**8) 功函数**

```text
IBRION = -1
NSW = 0
LDIPOL = .TRUE.
IDIPOL = 3
DIPOL = 0.5 0.5 0.5
LVHAR = .TRUE.       # LOCPOT 只含离子势 + Hartree，并把 LVTOT 置 .FALSE.
LVACPOTAV = .TRUE.   # >= 6.4.3，自动找无场区
```

`grep upper OUTCAR` 取 `vacuum level on the upper side and lower side of the slab`；`LVACPOTAV` 只支持 `IDIPOL` 1–3，偶极体系须同时开偶极校正。`Φ = V_vac − E_F`；平面平均：

```python
pot = pot.reshape((ngx, ngy, ngz), order="F")
v_z = pot.mean(axis=(0, 1))
z = np.arange(ngz) * c / ngz
```

`python plot_locpot_planar.py LOCPOT --outcar OUTCAR`（可加 `--vbm --cbm`）打印 `V_vac upper side: x.xxxx eV (std 0.00xx, slope +0.000x eV/A)`、`Phi upper`、`Delta Phi`。半导体报 `IP = E_vac − E_VBM`、`EA = E_vac − E_CBM`。偶极两侧差 `ΔΦ = e·p_z/(ε0·A) ≈ 180.95·p_z[e·Å]/A[Å²]` eV。

## 3. 关键标签表

| 标签 | 取值 | 何时用 / 注意 |
|---|---|---|
| `ISIF` | slab 一律 `2` | `3` 会优化 c 使真空塌缩 |
| `IBRION` / `EDIFFG` / `NSW` | `2` / `-0.02` eV/Å / `100` | 静态算 `IBRION = -1`、`NSW = 0` |
| `ISMEAR` / `SIGMA` | 金属 `1`（或 `2`）/ `0.10–0.20`；半导体 `0` / `0.05` | 法向 k=1 禁用 `ISMEAR = -5`；取 `energy(sigma->0)` |
| `LREAL` | `.FALSE.` | 实空间投影会给 slab 引入误差 |
| `LDIPOL` / `IDIPOL` / `DIPOL` / `EPSILON` | `.TRUE.` / `1–3`（分子 `4`）/ 质心分数坐标 / `1` | 只写 `IDIPOL` 只修总能，势与受力不变 |
| `LVHAR` | `.TRUE.` | LOCPOT = 离子势 + Hartree；同时把 `LVTOT` 置 `.FALSE.` |
| `LVACPOTAV` / `VACPOTFLAT` | `.TRUE.` / 默认 `0.1` eV/Å | ≥ 6.4.3；按平面平均电场找真空平台；配 `IDIPOL` |
| `WRT_POTENTIAL` | `hartree ionic` | ≥ 6.4.3 写 vaspout.h5；两者各自 G=0 定零，须看和 |
| `EFIELD` | eV/Å | 电子沿 `EFIELD` 方向移动；配 `LDIPOL`；< 6.6.0 不去反演对称 → 加 `ISYM = 0` |

## 4. 易错点

1. **`Tetrahedron method fails (number of k-points < 4)`** → slab 法向只 1 个 k 点却用 `ISMEAR = -5` → 改 `ISMEAR = 0`（半导体）或 `1`（金属）。
2. **`VERY BAD NEWS! internal error in subroutine IBZKPT: Reciprocal lattice and k-lattice belong to different class of lattices`** → 六角/斜方表面胞用了偶数 Monkhorst–Pack 网格 → 改 Γ 中心（如 `8 8 1`）。
3. **γ 随层数线性漂移、加层也不收敛** → 体相参考与 slab 设置不一致（k 点/smearing/ENCUT/LREAL）→ 用取向体相 `obulk` 或改 Fiorentini–Methfessel 拟合。
4. **真空里势不平、`WARNING: potential is not flat`** → 真空太薄或有净偶极 → 真空加到 15–20 Å 并开 `LDIPOL`。
5. **`Did not find any points to average over, which means that no vacuum field-free region was found. Please increase the size of  your cell`** → `LVACPOTAV` 找不到无场区 → 加真空、开偶极校正，或改 `LVHAR` 手动取平台。
6. **`The minimum charge density times volume of the cell along the axis of the dipole correction is larger 1E-1`** → 真空区仍有电荷/`DIPOL` 位置不当 → 居中 slab、`DIPOL` 放质心。
7. **带电 slab 加 `LDIPOL` 直接停** → 势修正只对立方超胞实现 → 非立方改用 `KERNEL_TRUNCATION`（≥ 6.5.0）。
8. **`EFIELD` 下结果与对称性不符** → < 6.6.0 不去反演对称（known issue 50）→ 加 `ISYM = 0`。
9. **平面平均曲线整体偏常数** → ASE `VaspChargeDensity` 读 LOCPOT 会按 CHGCAR 约定除体积 → 乘回体积；`WRT_POTENTIAL` 两项各自定零须相加。

## 5. 算完了怎么验收

- `grep "reached required accuracy" OUTCAR`、`grep "F=" OSZICAR | tail -3`：力已收敛、能量不再降。
- `CONTCAR` 层间距：体相 (111) 面间距 `d0 = 3.97/√3 = 2.292 Å`，未弛豫正是 2.292；弛豫后只 d23/d34 变、d12 不变，金属表面弛豫约 1%。最外层向内塌 > 5% ⇒ 固定层选错或 `ISIF` 用成 3。
- 表面能：γ(N) 在 N ≥ 7 后基本不变；`ε_fit` 与 `ε_bulk` 差 < 1 meV/atom；弛豫 γ 略低于未弛豫。PBE Cu(100) γ ≈ 1.4–1.5 J/m²（Tran et al., Sci. Data 3, 160080 (2016)），实验平均 ≈ 1.8 J/m²。
- 功函数：真空段 slope < 0.01 eV/Å、`grep -i "dipol" OUTCAR | tail -5` 有台阶值；PBE Pt(111) 应 5.6–5.9 eV（实验 5.9–6.1）；无吸附 clean slab 两侧 ΔΦ ≈ 0。

## 6. 本地完整输入

- $VASP_TUTORIALS_ROOT/tutorials/09-surface-workfunction.md
- $VASP_TUTORIALS_ROOT/examples/09-pt111-slab/
- $VASP_TUTORIALS_ROOT/examples/09-cu100-surface-energy/
- $VASP_TUTORIALS_ROOT/examples/09-pt111-workfunction-dipole/

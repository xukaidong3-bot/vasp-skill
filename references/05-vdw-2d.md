# 色散修正与二维材料（IVDW / LUSE_VDW）

> 什么时候读这个文件（触发场景）：算层状体相（石墨、h-BN、MoS2）的层间距与层间结合能、二维单层/双层的几何与带边、双层与异质结堆垛能；或 PBE 优化中 c 轴一直涨、不收敛。

## 1. 先判断走哪条路线

| 情形 | 做法 | 关键区别 |
|---|---|---|
| 层间距/堆垛能，要变胞或跑 MD | `IVDW = 12`（D3-BJ） | 有解析力与应力，可用于 `ISIF = 3`、有限位移声子；DFPT（`IBRION = 7/8`）不支持 |
| 报告级层间结合能，容忍慢 | optB86b-vdW、rev-vdW-DF2（`LUSE_VDW`） | 换写成 6 行 GGA/PARAM/AGGAC，**不再写 `IVDW`** |
| 泛函是 r2SCAN / HSE06 或要 D4 | `IVDW = 13`（需 `-DDFTD4`）或 r2SCAN+rVV10 | 内置 D3 只有 PBE/RPBE/revPBE/PBEsol/TPSS/M06L/SCAN/B3LYP/PBE0/HF 参数 |

## 2. 标准流程

**1) 体相变胞优化（石墨，PBE-D3(BJ)）**

```text
PREC = Accurate
ENCUT = 520          # >= 1.3 x ENMAX(C = 400)，变胞必须有此余量
EDIFF = 1E-6
LREAL = .FALSE.
ISMEAR = 0
SIGMA = 0.05
IVDW = 12            # D3(BJ)，PBE 参数已内置
IBRION = 2
ISIF = 3
NSW = 80
EDIFFG = -0.005      # 层间势能面平坦，判据从严
```

KPOINTS 必须 Γ 中心 `15 15 5`（六角胞用偶数 Monkhorst–Pack 网格报 IBZKPT 错）。查 `grep "reached required accuracy" OUTCAR` 与 `sed -n 3,5p CONTCAR`：c 应 6.6–6.8 Å、a ≈ 2.46–2.47 Å。纯 PBE 对照 c 涨到 8–9 Å 是物理结果，不是输入错。

**2) E–c 扫描复核**（c = 6.2–7.4 Å 取 7–9 点固定原子单点，拟合 E(c) 求极小，排除优化器提前停）

```python
from ase.io import read, write          # 固定 a，只缩放 c
atoms = read("D3BJ/relax/CONTCAR")
for c in np.arange(6.2, 7.41, 0.15):
    a = atoms.copy(); cell = a.cell.copy()
    cell[2] = cell[2] / np.linalg.norm(cell[2]) * c
    a.set_cell(cell, scale_atoms=True)  # 分数坐标不变
    write(f"POSCAR_c{c:.2f}", a, format="vasp", direct=True)
```

**3) 层间结合能**：参考态（单层）必须放在**同一面内晶格**、c = 20 Å 的同一晶胞里，面内 k 网格相同（石墨 15×15×5 → 石墨烯 15×15×1）、同 ENCUT、同 smearing。`E_b = E_mono/N_mono − E_bulk/N_bulk`（meV/atom；石墨 4 原子 + 石墨烯 2 原子即 `E_gr/2 − E_graphite/4`）；面积归一 `(N_L·E_mono − E_bulk)/(N_L·A)`（meV/Å²），1 meV/Å² ≈ 0.016 J/m²。别混「层间结合能」与「剥离能」。

**4) 二维单层：真空 + 只放开面内**（电子步同 1)，改这几行）

```text
ENCUT = 450
LASPH = .TRUE.
IVDW = 12            # 与要比对的双层/异质结完全一致
ISIF = 3
LATTICE_CONSTRAINTS = .TRUE. .TRUE. .FALSE.   # 只放开 a1 a2；需 >= 6.4.3
EDIFFG = -0.01
```

MoS2 取 c = 23.13 Å（≈20 Å 真空）、层放 z ≈ c/2，KPOINTS Γ 中心 `12 12 1`。检查 `head -5 relax/CONTCAR`：第三矢量仍为 `0 0 23.13`。

**5) 真空收敛**：c = 15/20/25/30 Å，原子笛卡尔坐标不动、重新居中（`set_vacuum.py`）。判据：相邻 c 的总能差 < 1 meV/atom；E_F 随晶胞变，跨 c 不可比。

**6) 带边对齐真空能级**：静态算 `LVHAR = .TRUE.`（LOCPOT 只含离子势 + Hartree 势），再 `plot_locpot_planar.py c_20/LOCPOT --outcar c_20/OUTCAR --vbm <VBM> --cbm <CBM>`，得 `IP = E_vac − E_VBM`、`EA = E_vac − E_CBM`。参考量级：单层 MoS2 面内 a ≈ 3.18–3.19 Å、K 点直接带隙 ~1.6–1.7 eV（无 SOC）。

**7) 双层结合能曲线**（刚性扫描单点：`IBRION = -1`、`NSW = 0`，其余同 4)；`make_bilayer.py --d 6.2 --a 3.19 --t 3.13 --c 30 > POSCAR` 建双层，`--mono` 建同晶胞单层参考，`bash run.sh collect` 汇总）。`E_b(d) = (E_bilayer(d) − 2E_mono)/A`，`A = (√3/2)a²`。判据：极小 d ≈ 6.1–6.3 Å（2H 体相 c/2 ≈ 6.15 Å），大 d 处 → 0，删 `IVDW` 的对照几乎无结合；得 d_eq 后在该距离做一次完整弛豫确认。

**8) 异质结失配**：`δ = (L_film − L_sub)/L_sub`，目标 |δ| ≲ 1–3%，写明应变加在哪层。实测例：MoS2 4×4 / 石墨烯 3(√3×√3)R30° → L = 12.76 Å、失配 0.46%、102 原子。通用搜索用 pymatgen `ZSLGenerator` / `CoherentInterfaceBuilder` 或 VASPKIT 804 / 825；不对称异质结必须偶极校正（见 06）。

**9) 二维库仑截断（选读，≥ 6.5.0）**

```text
KERNEL_TRUNCATION { LTRUNCATE = T   IDIMENSIONALITY = 2
                    ISURFACE  = 3   LCOARSEN = T }
```

`LTRUNCATE` 是总开关；原子不能落在法向晶胞边界（必须居中）。6.6.0 的二维截断有噪声势（known issue 72，6.6.1 修复）；中性对称二维材料用足够真空即可。

## 3. 关键标签表

| 标签 | 取值 | 何时用 / 注意 |
|---|---|---|
| `IVDW` | 1/10 = D2；11 = D3(0)；12 = D3(BJ)；13 = D4；15 = simple-DFT-D3；3 = DFT-ulg；4 = dDsC；2/20 = TS；21 = TS-HI；202 = MBD@rsSCS；263 = MBD@rSC/FI；14 = libMBD | 13 需 dftd4 + `-DDFTD4`（≥6.2；dftd4 ≥4.0 新 API 需 ≥6.6.0，≤3.7 用 `-DDFTD4_API_V3`）；15 需 ≥6.6.0 + `-DSDFTD3`；14 需 ≥6.4.3 + `-DLIBMBD`。`LVDW=.TRUE.` 已废弃（=`IVDW=1`） |
| `LUSE_VDW` + `GGA` + `LASPH` | `RE`/`ML`/`OR`/`BO`/`MK`/`CX`/`PE`/`BF`；`LASPH = .TRUE.` | 打开非局域关联；optB86b 用 `MK`；BEEF-vdW 用 `BF`（需 `-Dlibbeef`）；所有 vdW-DF 都建议开 `LASPH` |
| `PARAM1` `PARAM2` `AGGAC` `ZAB_VDW` | optB86b：0.1234 / 1.0；optB88：0.1833333333 / 0.22 | vdW-DF 一律 `AGGAC = 0.0`；vdW-DF2 / rev-vdW-DF2 / BEEF 另加 `ZAB_VDW = -1.8867`（默认 -0.8491） |
| `IVDW_NL` `BPARAM` `CPARAM` | 1 = Dion，2 = rVV10；默认 `BPARAM = 6.3`、`CPARAM = 0.0093` | r2SCAN+rVV10 → `BPARAM = 11.95`；SCAN+rVV10 → 15.7；PBE+rVV10L → 10 |
| `VDW_S6/S8`、`SDFTD3_DAMPING`、`LIBMBD_METHOD` | 保持默认；`SDFTD3_DAMPING` = zero / rational(BJ) / mzero / mrational / optimizedpower | 只在复现文献自定义参数时改；库版本会改结果（平滑截断），项目内固定 |
| `LATTICE_CONSTRAINTS` | `.TRUE. .TRUE. .FALSE.` | 二维变胞只放开面内；VASP ≤ 6.4.2 与 `IBRION = 1/2` 同用会被静默忽略 |
| `ENCUT` / `EDIFFG` | ≥ 1.3×ENMAX；−0.005 ~ −0.01 eV/Å | c 轴对 Pulay 应力与平坦势能面都敏感 |

## 4. 易错点

1. **c 一直涨、优化不收敛** → 半局域泛函没有层间 $r^{-6}$ 吸引 → 加 `IVDW` 或换 vdW-DF，别靠加大 `NSW` 硬等。
2. **vdW-DF 卡在第一步数小时** → ≤ 6.4.2 缺 `vdw_kernel.bindat`，正现场生成核函数 → 把核文件复制进目录，或升 ≥ 6.4.3（8 进程约 2 分钟）；rVV10 不需要核文件。
3. **`IVDW = 13` / `15` 报错退出或提示未编译** → 没链接 dftd4 / simple-dftd3 → 重编译，或退回内置 `IVDW = 12`。
4. **r2SCAN、optB86b 上 D3 结果异常** → 内置 D3 无这些泛函参数 → 换 D4 / simple-DFT-D3 或 vdW-DF。
5. **`LUSE_VDW` 与 `IVDW` 同时写** → 色散重复计入 → 二选一。
6. **二维变胞后真空塌缩** → 没约束 c，或 ≤ 6.4.2 中 `LATTICE_CONSTRAINTS` 被静默忽略 → 升 ≥ 6.4.3（推荐 6.5.1+），或固定 c 扫 a。
7. **`VERY BAD NEWS! internal error in subroutine IBZKPT: Reciprocal lattice and k-lattice belong to different class of lattices`** → 六角/二维胞用了偶数 Monkhorst–Pack 网格 → 改 Γ 中心、法向取 1。
8. **`Tetrahedron method fails (number of k-points < 4)`** → 二维/分子体系用了 `ISMEAR = -5` 且 k 点太少 → 改 `ISMEAR = 0`、`SIGMA = 0.05`。
9. **DFPT 声子不含 vdW** → `IBRION = 7/8` 不支持 `IVDW` → 改 `IBRION = 5/6` 或 phonopy。
10. **libMBD 优化中色散能不随结构变** → known issue 55（≤ 6.5.1 不传新坐标/晶胞）→ 升 6.6.0 或打补丁；known issue 91：k 网格按 `KSPACING = 0.5` 生成而非 KPOINTS。

## 5. 算完了怎么验收

- `grep -i "IVDW" OUTCAR | head -3`、`grep -i -E "edisp|dispersion|vdw" OUTCAR | head`：标签被读入、有色散能行、无 "not compiled"。
- `sed -n 3,5p CONTCAR`：石墨 c 落 6.6–6.8 Å、a 2.46–2.47 Å；二维单层第三矢量仍等于初始真空（未塌缩）。
- 结合能为正（分层参考减体相）；若为负、或大 d 处不趋于 0 ⇒ 参考态晶胞 / k 点 / ENCUT 不一致。
- 石墨 E_b 应 40–60 meV/atom（RPA ≈ 48）；双层 MoS2 极小 ≈ −20 meV/Å²、d_eq ≈ 6.1–6.3 Å；真空收敛表相邻 c 的 E0/atom 差 < 1 meV。参与相减的能量须同 `IVDW`/泛函/库版本，换过库版本就全部重算。

## 6. 本地完整输入

- $VASP_TUTORIALS_ROOT/tutorials/08-vdw-2d.md
- $VASP_TUTORIALS_ROOT/examples/08-graphite-vdw/
- $VASP_TUTORIALS_ROOT/examples/08-mos2-monolayer/
- $VASP_TUTORIALS_ROOT/examples/08-mos2-bilayer-binding/

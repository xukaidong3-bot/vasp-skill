# 吸附与催化：吸附能、ZPE/熵、CHE 台阶图

> 用户要表面吸附能、位点/覆盖度比较、ZPE 与熵校正、HER/OER/ORR 自由能台阶图与过电位、d 带中心或隐式水溶剂时读这里。

## 1. 先判断走哪条路线

| 情形 | 做法 | 关键区别 |
|---|---|---|
| 只要吸附强弱、位点排序 | 三个 OUTCAR 取最后一个 `energy(sigma->0)`：`E_ads=E(slab+ads)-E(slab)-n*E(gas)` | 不用频率；三处电子参数须逐项一致 |
| 要有限温度自由能/电催化台阶图 | 上式 + `IBRION=5` 频率 → ZPE、TS → `G=E_DFT+G_corr` → CHE | 吸附物谐振子，气相理想气体 |
| 体系在水/电解液里 | VASPsol 两步：真空 `LWAVE=.TRUE.` 存 WAVECAR → `ISTART=1`+`LSOL=.TRUE.` | 需打补丁重编译 |
| 关心覆盖度 | p(2×2)=1/4 ML → p(3×3)=1/9 ML → p(4×4)，E_ads 变化 <0.02–0.05 eV 即收敛 | 换催化剂必须同覆盖度比较 |

## 2. 标准流程

**步骤 1 清洁 slab**（底部 2 层固定，`ISIF=2` 绝不变胞）

```bash
PREC=Accurate   ENCUT=450   EDIFF=1E-6   NELM=150
ISMEAR=1        SIGMA=0.10          # 金属 Methfessel-Paxton
LREAL=.FALSE.   LASPH=.TRUE.
LDIPOL=.TRUE.   IDIPOL=3   DIPOL=0.5 0.5 0.5
IBRION=2        ISIF=2     NSW=200   EDIFFG=-0.02
```
查：`grep "reached required accuracy" OUTCAR`；`grep "entropy T" OUTCAR | tail -1` 每原子 <1 meV；`grep "dipolmoment" OUTCAR` 有输出。六角胞 KPOINTS 用 `Gamma`。

**步骤 2 气相参考**（15.0×15.5×16.0 Å 盒子、Γ 点、`VASP_EXE=vasp_gam`）

```bash
ISMEAR=0   SIGMA=0.01       # 有 HOMO-LUMO 隙
ISPIN=2    MAGMOM=2*1       # 仅 O2 / 自由基
ISIF=2                      # 分子盒绝不变胞
```
查：C–O ≈1.14 Å（PBE）；`grep magnetization OUTCAR` O2 为 2.0。

**步骤 3 吸附优化 + 吸附能**（初始高度 C 到表面层：top 1.85 Å、bridge 1.50 Å、fcc/hcp 1.30 Å；INCAR 同步骤 1）

```bash
python ads_energy.py --slab 1_slab --ads 3_co_top 4_co_fcc --gas 2_co_gas
python ads_energy.py --slab slab --ads H_fcc --gas H2_gas -n 0.5   # H2->2H*，按每 H
```
查：top 位 1 个 Pt–C ≈1.85 Å、fcc 位 3 个 ≈2.1 Å；top 位出现 2–3 个近邻 Pt 说明滑到 bridge/hollow。

**步骤 4 频率（只放开吸附质）**

```bash
python set_freq_selective.py 3_co_top/CONTCAR -e C O -o 5_freq_co_top/POSCAR
```
其余同步骤 1，改：
```bash
EDIFF=1E-7      # 力的差分对 SCF 噪声敏感
IBRION=5        # 唯一支持 Selective dynamics 的模式
POTIM=0.015     NFREE=2     NSW=1
```
查：`grep "THz" OUTCAR | head -6` 得 6 个模，C–O 伸缩 ~2000–2100 cm⁻¹。频率在 OUTCAR 打印两遍，只取前 3N_free 行。

**步骤 5 ZPE / 熵 / G_corr**

```bash
python calc_dG.py thermo 5_freq_co_top -T 298.15
python calc_dG.py thermo freq_H2 --gas --geometry linear --symmetry 2 --spin 0 --structure freq_H2/POSCAR
```
输出 ZPE、∫Cv dT、T\*S，`G_corr = ZPE + ∫Cv dT - TS`；吸附物默认把 <50 cm⁻¹ 低频抬到 50 cm⁻¹（`--cutoff 0` 关闭）。298 K 量级：H2 ZPE ~0.27 / TS ~0.40–0.41 eV；H2O(l) ~0.56 / ~0.67 eV；OH\* ~0.35、O\* 0.05–0.07、OOH\* 0.40–0.45 eV。不同来源校正值不可混用。

**步骤 6 CHE 自由能与过电位**

HER：`ΔG_H* = ΔE_H* + ΔE_ZPE - TΔS`，`ΔE_H* = E(slab+H) - E(slab) - ½E(H2)`；常用 `ΔG_H* ≈ ΔE_H* + 0.24 eV`（Nørskov 2005），|ΔG_H\*|→0 为火山顶。

```bash
python calc_dG.py her --slab ../10-co-pt111-adsorption/1_slab --slabH H_fcc --h2 H2_gas --corr 0.24
python calc_dG.py her ... --gcorr-h <G_corr(H*)> --gcorr-h2 <G_corr(H2)>
python calc_dG.py oer my_oer.json --U 1.23 --plot oer.png
```
OER 四步（U=0 V_RHE）：`ΔG1=G_OH*-G_*-G_H2O+½G_H2`、`ΔG2=G_O*-G_OH*+½G_H2`、`ΔG3=G_OOH*-G_O*-G_H2O+½G_H2`、`ΔG4=4.92-(ΔG1+ΔG2+ΔG3)`。`U_L=max ΔG_i`，`η=U_L-1.23`；ORR 反之 `U_L=min ΔG_i`，`η=1.23-U_L`。标度关系 `ΔG_OOH*-ΔG_OH*≈3.2 eV`（理想 2.46）→ 最小 η≈0.37 V。JSON 的 `E` 可写数字或目录路径，`G_corr` 逐个给。

**步骤 7 d 带中心 / 隐式溶剂**

```bash
python dband_center.py 7_slab_dos/DOSCAR -a 13-16 -w -10 2
bash run_vaspsol.sh        # vac -> sol；`bash run_vaspsol.sh lpb` 再加电解质
```
DOS 段要静态：`IBRION=-1`、`NSW=0`、`LORBIT=11`、`NEDOS=3001`、`EMIN=-15`、`EMAX=10`、`ISMEAR=-5`、k 6×6×1→12×12×1。VASPsol 末行会打印 `E(solvent) - E(vacuum)`。

## 3. 关键标签表

| 标签 | 取值 | 何时用 / 注意 |
|---|---|---|
| `ENCUT` | 450（≥1.1×最大 ENMAX） | slab/吸附/气相三处必须相同 |
| `ISMEAR`/`SIGMA` | 金属 1/0.10–0.20；分子 0/0.01；DOS −5 | DOS 用 −5 时 k 点须 ≥4 |
| `EDIFF` | 优化 1E-6；频率 1E-7 | 力不准会让优化器和频率虚振 |
| `EDIFFG` | 优化 −0.02；频率前 −0.01 | 负值 = 力判据 |
| `IBRION` | 2 优化；5 频率；−1 静态 | 6 忽略 Selective dynamics |
| `POTIM`/`NFREE` | 0.015 / 2（中心差分） | 轻原子可试 0.01–0.02 |
| `LREAL`/`LASPH` | `.FALSE.`/`.TRUE.` | ≤50 原子用 `.FALSE.`；大 slab 可 `Auto`，但全程一致 |
| `LDIPOL`/`IDIPOL`/`DIPOL` | `.TRUE.`/3/`0.5 0.5 0.5` | 不对称 slab 必开 |
| `AMIN` | 0.01 | 偶极校正引起 SCF 振荡时 |
| `ISPIN`/`MAGMOM` | 2/`2*1`（O2）；`NUPDOWN=2` 备选 | 自由基、磁性表面同样要开 |
| `IVDW` | 12（D3(BJ)） | 物理吸附必加；一旦用就全流程一起用 |
| `LORBIT`/`NEDOS`/`EMIN`/`EMAX` | 11/3001/−15/10 | d 带中心静态计算 |
| `LSOL`/`EB_K`/`TAU` | `.TRUE.`/78.4/0.000525 | VASPsol 水默认值 |
| `NC_K`/`SIGMA_K`/`LAMBDA_D_K` | 0.0025/0.6/3.0（≈1 M 电解质） | 对水拟合，换溶剂不能只改 `EB_K` |

## 4. 易错点

1. **现象** 吸附能为正且很大（>+1 eV）→ **原因** 三个计算的 POTCAR/ENCUT/PREC/LASPH/LREAL/GGA/IVDW 有一项不同，或气相没优化/自旋错，或取了 TOTEN → **处置** `grep TITEL POTCAR`、`grep ENCUT OUTCAR` 逐项对照后重算。
2. **现象** 报 `VERY BAD NEWS! internal error in subroutine IBZKPT: Reciprocal lattice and k-lattice belong to different class of lattices` → **原因** 六角表面胞用了偶数 MP 网格 → **处置** KPOINTS 改 `Gamma`。
3. **现象** 报 `Tetrahedron method fails (number of k-points < 4)` → **原因** 分子/Γ 点算用了 `ISMEAR=-5` → **处置** 改 `ISMEAR=0`、`SIGMA=0.01`。
4. **现象** O2 的 magnetization 为 0 → **原因** 未设自旋或 `MAGMOM` 太小 → **处置** `ISPIN=2`+`MAGMOM=2*1`，必要时 `NUPDOWN=2`。
5. **现象** 频率出现 2–3 个 50–150 cm⁻¹ 虚频 → **原因** 结构没优化到位、`EDIFF` 太松、`POTIM` 不合适 → **处置** `EDIFFG=-0.01` 重优化、`EDIFF=1E-7`、`POTIM` 试 0.01/0.02。
6. **现象** 频率多出一堆近零模 → **原因** slab 原子也被放开，或用 `IBRION=6`（放开所有原子）→ **处置** `set_freq_selective.py` + `IBRION=5`。
7. **现象** 开 `LDIPOL` 后能量在两值间振荡 → **原因** 电荷在 slab 两侧晃动 → **处置** `AMIN=0.01`、`NELM=200`。
8. **现象** VASPsol 溶剂化能量与真空完全相同 → **原因** 可执行文件没打补丁，或没走两步法 → **处置** `grep -i -m 5 -E "solv|EB_K|TAU|Ediel" sol/vasp.out sol/OUTCAR` 确认。
9. **现象** d 带中心与文献差 0.5 eV 以上 → **原因** 积分窗口/选的原子/NBANDS 不同 → **处置** 统一窗口（如 [−10, E_F]）并注明。

## 5. 算完了怎么验收

- 每个优化目录的 `OUTCAR` 都要有 `reached required accuracy`，否则吸附能不可用。
- 相减只用最后一个 `energy(sigma->0)`（不是 TOTEN）；`grep "entropy T"` 每原子应 <1 meV。
- `ads_energy.py` 三次能量都标 `relaxed`，E_ads 为负且量级合理（PBE 下 1/4 ML CO/Pt(111) 约 −1.6 到 −1.9 eV）；出现 `NOT CONVERGED` 或正值即算错。
- 频率：吸附物 6 个实模无虚频，C–O 伸缩低于气相（气相 PBE 约 2100 cm⁻¹）；气相线形分子只有 1 个真模、其余 5 个近零（热化学要剔除）。
- `calc_dG.py` 的 `G_corr` 与量级表差 >0.1 eV 时先查频率是否完整、`--cutoff`/`--gas` 模型是否选错。
- CHE：`η` 与 `--plot` 台阶图逐段对得上；四步之和须为 4.92 eV。
- VASPsol：`sol/` 的 SCF 步数明显多于 `vac/`，且`E(solvent)-E(vacuum)` 非零；两者几乎相同 = 溶剂没生效。

## 6. 本地完整输入

$VASP_TUTORIALS_ROOT/examples/10-co-pt111-adsorption/
$VASP_TUTORIALS_ROOT/examples/10-her-oer-free-energy/
$VASP_TUTORIALS_ROOT/examples/10-vaspsol-implicit-solvent/
$VASP_TUTORIALS_ROOT/tutorials/10-adsorption-catalysis.md
$VASP_TUTORIALS_ROOT/scripts/

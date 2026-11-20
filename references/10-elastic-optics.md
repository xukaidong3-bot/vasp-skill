# 弹性常数与介电/光学响应

> 什么时候读这个文件：用户要 Cij、B、G、E、ν、德拜温度，或 ε∞/ε0、Born 电荷 Z*、n、k、α、ε(ω)、Tauc 带隙、GW/BSE。晶格振动去 `09-phonon.md`。

## 1. 先判断走哪条路线

| 情形 | 做法 | 关键区别 |
|---|---|---|
| 三维弹性常数（默认） | 弛豫到零压力后 `IBRION=6`、`ISIF=3`、`NFREE=2`、`POTIM=0.015` → 读 `TOTAL ELASTIC MODULI (kBar)` | 一次给全部 Cij，最快；看应力精度 |
| 胞小 / 应力噪声大 / 二维 | 能量–应变：造 ±0.5%–±1.5% 应变拟合 ΔE–δ（或 VASPKIT `201`）；二维出面内刚度 **N/m** | 慢，但对 `ENCUT`/k 点不敏感；不能直接读 kBar 的 C33；1 eV/Å² = 16.0218 N/m |
| 静态 ε∞ + Z* | LDA/GGA：`LEPSILON=.TRUE.`+`LRPA=.FALSE.`+`IBRION=-1`；杂化：`LCALCEPS=.TRUE.` | 前者不占空带、不支持 vdW、必须无 `NCORE/NPAR`；后者 k 网格必须含 Γ |
| ε(ω)、n、k、α | 两步：SCF 留 `WAVECAR` → `LOPTICS=.TRUE.` 非自洽 | 独立粒子近似：无激子、无间接跃迁的声子辅助 |
| 准粒子带隙 / 激子 | `ALGO=EVGW0` → `ALGO=BSE` | 成本是 LOPTICS 的 100–1000 倍 |

## 2. 标准流程

### A. 弹性常数（IBRION=6）

**1) 弛豫到"零应力"**：

```
PREC = Accurate
ENCUT = 400        # ≥1.3×ENMAX（Al 240 eV）
EDIFF = 1E-8
EDIFFG = -1E-3
ISMEAR = 1
SIGMA = 0.10
LREAL = .FALSE.
IBRION = 2
ISIF = 3
NSW = 100
```
查：`grep "external pressure" OUTCAR|tail -1` 收到 **<1 kB**；`grep "entropy=" OSZICAR|tail` 的 T·S 应 <1 meV/atom。`cp CONTCAR POSCAR` 再跑一轮直到压力不变。

**2) 弹性计算**：

```
# PREC/ENCUT/EDIFF/EDIFFG/ISMEAR/SIGMA/LREAL 与弛豫完全一致
IBRION = 6
ISIF = 3
NFREE = 2          # 中心差分 ±POTIM；4 再加 ±2×POTIM，更准但代价翻倍
POTIM = 0.015
NSW = 1
LWAVE = .FALSE.
LCHARG = .FALSE.
```
查：OUTCAR 三块 `SYMMETRIZED ELASTIC MODULI (kBar)`（钳位离子）、`ELASTIC MODULI CONTR FROM IONIC RELAXATION (kBar)`（内应变贡献）、`TOTAL ELASTIC MODULI (kBar)`（两者之和，**通常引用的就是这组**）。

**3) 后处理**（Voigt 重排 + kBar→GPa 都不能省）：
```bash
python voigt_vrh.py -f OUTCAR --poscar CONTCAR   # Born 判据+VRH+德拜；--clamped 取钳位组
python elastic_post.py -f OUTCAR --poscar CONTCAR --matrix Cij_GPa.txt
```
VASP 顺序 **XX YY ZZ XY YZ ZX**，标准 Voigt 是 **xx yy zz yz zx xy**，须按 `[0,1,2,4,5,3]` 重排再 `/10.0`（1 kBar = 0.1 GPa）；直接抄 OUTCAR 矩阵算 VRH 是错的。

**4) 判据**：6 个本征值全 >0。立方：`C11−C12>0`、`C11+2C12>0`、`C44>0`；六方：`C11>|C12|`、`2C13²<C33(C11+C12)`、`C44>0`、`C66>0`。

### B. 能量–应变（交叉验证 / 二维）

`ISIF=2`（二维固定 c；三维用 `3`）、`IBRION=-1`+`NSW=0` 单点、`EDIFF=1E-8`，每个应变都从同一 `CONTCAR` 出发：`(δ,δ,δ)`→`ΔE/V0=(3/2)(C11+2C12)δ²`；`(δ,−δ,0)`→`(C11−C12)δ²`；`(0,0,0,δ,δ,δ)`→`(3/2)C44δ²`。应变 **±0.5%–±1.5%**，别超 2%；拟合抛物线奇次项应≈0，否则参考结构不在极小。

二维跑批与拟合：`bash run_strain.sh`（内部 `MODES="xx yy biaxial shear"`、`STRAINS="-0.015 -0.010 ... 0.015"`）→ `python fit_2d_elastic.py`。
换算：`C2D(N/m)=C3D(GPa)×L_z(Å)/10`；六方 `Y=(C11²−C12²)/C11`、`ν=C12/C11`，正交 `Y=(C11C22−C12²)/C22`、`ν=C12/C22`；判据 `C11>0, C22>0, C11C22>C12², C66>0`。MoS2 量级 C11≈130–150 N/m、Y≈120–140、ν≈0.25。

### C. 介电与光学

**1) SCF**：`PREC=Accurate`、`ENCUT=350`、`EDIFF=1E-8`、`ISMEAR=0`、`SIGMA=0.01`、`LREAL=.FALSE.`、`LWAVE=.TRUE.`、`LCHARG=.TRUE.`；k 点要比常规密（示例 Γ 16×16×16）。

**2) 非自洽 + LOPTICS**（上一段的 `PREC/ENCUT/EDIFF/ISMEAR/SIGMA/LREAL` 照搬）：

```
ISTART = 1
ALGO = Exact
NELM = 1           # 只重算光学矩阵元
NBANDS = 32        # Si2 默认仅 8；需 2-3 倍（24/32/48 做收敛测试）
LOPTICS = .TRUE.
CSHIFT = 0.100
NEDOS = 2000
```
查：`vasprun.xml` 有 `<dielectricfunction>`；`grep -i NBANDS OUTCAR` 确认没被截断。

**3) 后处理**：

`python optics_post.py optics --kind current`（ε1,ε2,n,k,α,R,L）→ 带隙用 `--tauc direct --fit 3.0 3.6` → `python plot_optics.py -f run/vasprun.xml` 出 `optics.dat` 与图（各向异性加 `--component zz`）。

**4) 静态 ε + Born 电荷**：单胞 `LEPSILON=.TRUE.`、`LRPA=.FALSE.`、`IBRION=-1`、`NSW=0`。读 `MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects in DFT)`（=ε∞）与 `BORN EFFECTIVE CHARGES`（每原子 3×3，ΣZ*≈0）；离子贡献在 `IBRION=8` 的 `... IONIC CONTRIBUTION`，满足 ε0=ε∞+ε_ion，并可用 LST `ε0/ε∞=(ω_LO/ω_TO)²` 交叉验证。

**5) 杂化**：`LHFCALC=.TRUE.`、`HFSCREEN=0.2`、`AEXX=0.25`、`LCALCEPS=.TRUE.`（ZnO 用 `ENCUT=520`、`NBANDS=24`、`EDIFF=1E-7`）。

**6) GW/BSE**：GW 用 `ALGO=EVGW0`（不写 `NBANDS`，≥6.3 自动取）、**保留 `LOPTICS=.TRUE.`** 写 `WAVEDER`、收敛 `ENCUTGW`/`NOMEGA`，看 `grep -A 20 "QP shifts" OUTCAR`；BSE 用 `ALGO=BSE`+`NBANDSO=4`、`NBANDSV=8`、`OMEGAMAX=20`、`ANTIRES=0`（TDA），需先有 `Wxxxx.tmp`。自检：关 `LADDER`/`LHARTREE` 应重现 IPA 谱。

## 3. 关键标签表

| 标签 | 取值 | 何时用 / 注意 |
|---|---|---|
| `IBRION` | `6`/`8`/`-1` | `6` 靠对称性、**不支持 Selective dynamics**（那是 `5`）；`LHFCALC=.TRUE.` 时 `ISIF` 默认变 0 |
| `LEPSILON`/`LRPA` | `.TRUE.`/`.FALSE.` | `LRPA=.FALSE.` 才给含局域场的 Z*；仅 LDA/GGA，无 vdW |
| `LCALCEPS`/`EFIELD_PEAD` | `.TRUE.`/`0.01 0.01 0.01` | 杂化唯一途径；k 网格必须含 Γ |
| `LOPTICS`/`CSHIFT`/`NEDOS` | `.TRUE.`/`0.100`/`2000` | `SIGMA` 必须 ≪ `CSHIFT`；`NBANDS` 取默认 2–3 倍 |

## 4. 易错点

1. **Cij 随 `ENCUT`/k 点漂移几 GPa** → 应力未收敛 → `ENCUT` 取 1.3×ENMAX 以上再做 ±15% 步长测试，k 点加密到总能变化 <1 meV/atom。
2. **`TOTAL` 与 `SYMMETRIZED` 差很远** → 结构没到零压力或内应变贡献异常 → 先 `ISIF=3` 跑到 |P|<1 kB 再算。
3. **六角/三角体系 Cij 不合理** → 已知 issue 7（<6.4.3）→ 升级，或临时 `KBLOWUP = .FALSE.`。
4. **`LEPSILON` 报 `PEAD routines do not work for NCORE`（或提示 `remove the tag NPAR from the INCAR file`）** → 线性响应不能做能带并行 → 删 `NCORE`/`NPAR`，并行用 `KPAR`。
5. **Born 电荷之和不为零** → k 点太稀 / 网格不含 Γ / 体系带电（issue 92 把 Z* 置零）→ 加密 k 点并含 Γ。
6. **`no <dielectricfunction> found: was LOPTICS = .TRUE. set and did VASP finish?`** → 没设 `LOPTICS`、第二步没读到 `WAVECAR`（`ISTART=0` 会重跑 SCF 而不写谱）、或作业没跑完。
7. **光谱低能端假峰 / 带隙算成 0** → `SIGMA` 接近或超过 `CSHIFT`，或用了不支持 LOPTICS 的 `ISMEAR=-2` → 改 `ISMEAR=0` 且 `SIGMA ≪ CSHIFT`。
8. **Si 上用间接 Tauc 拟出带隙** → IPA 无声子辅助间接跃迁、无激子 → 间接带隙的 `(αhν)^2` 拟合无意义；直接带隙用 `(αhν)^(1/2)`。
9. **`LOPTICS` 与 `LEPSILON` 混在一个 INCAR 里互相污染** → 两者独立 → 拆两步（ε∞ 走 `LEPSILON`，ε(ω) 走 `LOPTICS`）。

## 5. 算完了怎么验收

- 弹性：`TOTAL ELASTIC MODULI (kBar)` 的 6×6 矩阵，重排 + `/10` 后 6 个本征值全正；立方核 `C11−C12>0`、`C11+2C12>0`、`C44>0`。派生量看 `voigt_vrh.py` 的 VRH B/G、`E=9BG/(3B+G)`、`ν=(3B−2G)/(2(3B+G))`、Pugh `B/G`（>1.75 延性）、`A^U`（0=各向同性）、立方 Zener（1=各向同性）。与能量–应变路线应差 5% 以内。
- 二维：`fit_2d_elastic.py` 的 R²>0.999、奇次项≈0；MoS2 的 C11（N/m）落在 130–150。
- 介电：`BORN EFFECTIVE CHARGES` 存在且 ΣZ*≈0；ε∞ 对实验（Si≈11.9、GaAs≈10.9、ZnO≈3.7）；ε0=ε∞+ε_ion 自洽。
- 光学：`optics.dat` 的 ε2 峰位对已知临界点（Si 的 E1/E2 约 3.4 与 4.2 eV）；`α` 带隙以上为正、以下≈0；`NBANDS` 32→48、`NEDOS` 加倍后峰移 <0.05 eV。锯齿噪声大 = `NEDOS` 太小或 k 点太稀；0 eV 处非物理巨值 = 未收敛或已金属化。

## 6. 本地完整输入

- `$VASP_TUTORIALS_ROOT/examples/13-elastic-al-ibrion6/`
- `$VASP_TUTORIALS_ROOT/examples/13-elastic-mos2-2d/`
- `$VASP_TUTORIALS_ROOT/examples/14-optics-si-loptics/`
- `$VASP_TUTORIALS_ROOT/examples/14-dielectric-gaas-lepsilon/`
- `$VASP_TUTORIALS_ROOT/examples/14-dielectric-zno-hse-lcalceps/`
- `$VASP_TUTORIALS_ROOT/scripts/elastic_post.py`（另见同目录 `optics_post.py`）
- `$VASP_TUTORIALS_ROOT/tutorials/13-elastic.md`
- `$VASP_TUTORIALS_ROOT/tutorials/14-optics-dielectric.md`

# 电荷与成键分析（Bader / 差分电荷 / ELF / 部分电荷与 STM / LOBSTER COHP）

> 什么时候读这个文件：用户问"每个原子带多少电荷、电荷往哪流、到底成不成键、键有多强、某个带或缺陷态长什么样、模拟一张 STM 图"，或 Bader/ELF/COHP 报错要排错时。力学与热力学性质不读这里。

## 1. 先判断走哪条路线

| 用户的问题 | 做法 | 关键区别 |
|---|---|---|
| 每个原子净电荷、离子性 | Bader：`LAECHG=.TRUE.` → `chgsum.pl` → `bader -ref` | 必须有全电子参考，缺了金属原子会算出 −3 e |
| 电荷往哪流、成键区聚集/流失 | 差分电荷：AB、片段、分子三次单点相减 | 三次同晶胞同网格，片段绝不重优化 |
| 共价/离子/金属、孤对、断键 | ELF：`LELF=.TRUE.` | 0–1 无量纲，等值面 0.75–0.85，`KPAR` 必须 1 |
| 某能量窗口/某几条带的空间分布、STM | `LPARD` + `EINT`/`IBAND` → PARCHG | 两步：先 `LWAVE` 自洽，再用 WAVECAR 后处理 |
| 成键/反键、键强排序 | LOBSTER COHP/COOP/COBI | 未约化 k 集（`ISYM=0/-1`）、基组覆盖、spilling < 5% |

前提：都在已优化结构上做单点（`NSW=0`、`IBRION=-1`）。

## 2. 标准流程

### A. Bader

1. 自洽单点，网格取默认的 1.5–2 倍：

```
PREC=Accurate  ENCUT=400  EDIFF=1E-6  ALGO=Fast  LREAL=Auto  ISMEAR=1  SIGMA=0.1
NELM=120  NSW=0  IBRION=-1  LCHARG=.TRUE.  LAECHG=.TRUE.  LWAVE=.FALSE.  NCORE=4  KPAR=2
NGXF=160  NGYF=160  NGZF=768
```

2. POTCAR 顺序与 POSCAR 一致；跑完 `ls -lh CHGCAR AECCAR0 AECCAR2` 都在。
3. 合并全电子密度再积分：

```
chgsum.pl AECCAR0 AECCAR2      # → CHGCAR_sum
bader CHGCAR -ref CHGCAR_sum   # → ACF.dat BCF.dat AVF.dat
tail -3 ACF.dat                # NUMBER OF ELECTRONS
```

4. 检查：`NUMBER OF ELECTRONS` = NELECT（示例 170 = 16×10+4+6），偏差 < 0.01 e；`MIN DIST` 大于原子芯半径；净电荷 = ZVAL − `CHARGE`。

### B. 差分电荷

5. AB / slab / co 三份目录共用同一 INCAR：与 A 相同，只把 `ISPIN=1`、网格写成 `NGXF=120 NGYF=120 NGZF=576`（同 `ENCUT`/`PREC` 才同网格）。
6. 片段坐标原样切：`python split_fragments.py POSCAR --keep Pt --out slab/POSCAR`、`--keep C O --out co/POSCAR`。
7. 相减：`chgdiff.pl` 连做两次，或 `chgsum.pl A B 1.0 -1.0`，或用仓库脚本（先校验网格/晶胞/位移，并给平面平均 `cdd_planar.dat`）：

```
python $VASP_TUTORIALS_ROOT/scripts/chgdiff.py diff AB/CHGCAR slab/CHGCAR mol/CHGCAR -o CHGCAR_diff
```

8. 看图：VESTA 等值面 ±0.002 e/Å³，黄=聚集、青=流失。

### C. ELF

9. 在 A 的单点上改三处：`LELF=.TRUE.`、`KPAR=1`、`NPAR=1`，`PREC` 保持 `Accurate`。VESTA 打开 ELFCAR，等值面 0.75–0.85。

### D. 部分电荷与 STM

10. 先自洽并留波函数（同 A，改）：`PREC=Accurate ENCUT=400 EDIFF=1E-7 ISMEAR=0 SIGMA=0.05 NBANDS=24 NSW=0 IBRION=-1 LORBIT=11 LWAVE=.TRUE. LCHARG=.TRUE.`。KPOINTS 用 Γ 12×12×1；`grep -E "NELECT|E-fermi|val. band max" OUTCAR` 核对带边（Mo_pv 12 + S 6×2 = 24 → 12 条占据带、CBM 第 13 条）。
11. 复制 POSCAR/KPOINTS/POTCAR/WAVECAR，`ENCUT`/`ISPIN`/`ISYM`/`NBANDS` 与自洽一致，只改：

```
ISTART=1  LPARD=.TRUE.  NBMOD=-3  EINT=-1.0 0.0   # 占据态窗口，相对 E_F
LSEPB=.FALSE.  LSEPK=.FALSE.                      # 合成单个 PARCHG
```

空态用 `EINT=1.0 3.0`（按自己 CBM 调）。查 `awk 'NR>15 && NF==5' PARCHG | head` 不全为 0。
12. 出图：`python stm_from_parchg.py PARCHG --mode height --height 3.0 --repeat 6 -o stm_h.png`（恒流 `--mode current`）。

### E. LOBSTER

13. VASP 侧：`ISYM=0`、`LWAVE=.TRUE.`、`NBANDS` ≥ 局域基函数总数且能被能带组数 (ranks/KPAR/NCORE) 整除、必须 `vasp_std`。GaAs 中 Ga_d(3d4s4p)=9 + As(4s4p)=4 = 13 → `NBANDS=16`：

```
PREC=Accurate  ENCUT=400  EDIFF=1E-7  ISMEAR=0  SIGMA=0.05  NSW=0  IBRION=-1
ISYM=0  NBANDS=16  LWAVE=.TRUE.  LCHARG=.TRUE.  NEDOS=2001  LORBIT=11  NCORE=4
```

14. 同目录写 lobsterin，`export OMP_NUM_THREADS=16` 后 `lobster > lobster.out`：

```
COHPstartEnergy -20.0
COHPendEnergy 5.0
gaussianSmearingWidth 0.05
basisSet pbeVaspFit2015
basisfunctions Ga 3d 4s 4p
basisfunctions As 4s 4p
cohpGenerator from 2.0 to 2.6 type Ga type As orbitalwise
saveProjectionToFile
```

## 3. 关键标签表

| 标签 | 取值 | 何时用 / 注意 |
|---|---|---|
| `LAECHG` | `.TRUE.` | Bader 前提；出 AECCAR0（芯）与 AECCAR2（自洽后全电子价） |
| `NGXF/NGYF/NGZF` | 默认 1.5–2 倍（示例 160/160/768） | 三份必须相同；太粗电子数对不上 |
| `LCHARG` | 默认 `.TRUE.` | CHGCAR 存 ρ·V_cell，积分 = ZVAL 之和（非原子序数）；`ISPIN=2` 第二块是 ρ↑−ρ↓ |
| `LELF` / `KPAR` | `.TRUE.` / `1` | `KPAR>1` 不支持；`NPAR` 也设 1 |
| `NBMOD` | `n>0`：IBAND 的 n 条带；`0`：全部带；`-1`：全部占据带且只写 CHGCAR（默认！）；`-2`：EINT 绝对能量；`-3`：EINT 相对 E_F | 只要新 CHGCAR 没 PARCHG = 忘改它 |
| `IBAND`/`KPUSE` | 带号/k 点号列表 | 设 IBAND 后 NBMOD 自动等于列表长度；KPUSE 把 k 权重重置为 1 |
| `ISYM` | `0`（LOBSTER、重算 PARCHG） | LOBSTER 要未按空间群约化的 k 集 |
| `NBANDS` | ≥ 基函数总数，且能被能带组数整除 | 不整除时 VASP 自动加大并打印 `The number of bands has been changed` |

## 4. 易错点

1. **现象**：金属原子 Bader 电荷是 −3 e 之类怪值。→ 原因：直接积分价电子 CHGCAR，没做全电子参考。→ 处置：`LAECHG=.TRUE.` + `chgsum.pl AECCAR0 AECCAR2` 后 `bader CHGCAR -ref CHGCAR_sum`。
2. **现象**：ACF.dat 电子数 ≠ NELECT。→ 原因：网格太粗、`-vac`、积分错的文件。→ 处置：加 `NGXF/NGYF/NGZF`；去掉 `-vac`（1.03 起默认 off）；确认 `-ref CHGCAR_sum`。
3. **现象**：加密网格后电荷还在变。→ 原因：网格未收敛或权重算法。→ 处置：继续加密，用 `-b weight` 检验。
4. **现象**：差分电荷核附近正负伪影、不守恒。→ 原因：片段被重新优化（移 0.01 Å 就出伪影），或三次 `ENCUT`/`PREC`/晶胞不同。→ 处置：坐标原样复制；三份 INCAR 写同一组 `NGXF/NGYF/NGZF`。
5. **现象**：`WARNING: dimensions on CHGCAR file are different`。→ 原因：续算或 `ICHARG=1/11` 时网格与旧 CHGCAR 不符。→ 处置：统一网格或删旧 CHGCAR。
6. **现象**：ELF 报 `ELF: KPAR>1 not implemented`；非共线报 `WARNING: ELF not implemented for non collinear case`。→ 处置：`KPAR=1`、`NPAR=1`。
7. **现象**：LPARD 跑完只有新 CHGCAR 没 PARCHG，或 PARCHG 全 0。→ 原因：前者 `NBMOD` 默认 −1；后者窗口落在带隙、EINT 当绝对能量、`NBANDS` 不足。→ 处置：显式 `NBMOD=-3`；按带边定窗口；加 `NBANDS`。
8. **现象**：LPARD 读 WAVECAR 报错或结果与自洽不符。→ 原因：`ENCUT`/`NBANDS`/晶胞/k 点与自洽不同，或 `KPUSE` 重置权重。→ 处置：全部一致；固定 k 时自洽与后处理都设 `ISYM=-1`。
9. **现象**：lobsterout 的 `abs. total spilling` > 5%（原文形如 `abs. total spilling:   1.xx%`、`abs. charge spilling:  0.xx%`）。→ 原因：基组不全、`NBANDS` 不足、POTCAR 不匹配。→ 处置：补 `basisfunctions`、加 `NBANDS`（经验 <1–2%）。

## 5. 算完了怎么验收

- **Bader**：ACF.dat 末 `NUMBER OF ELECTRONS` 与 NELECT 差 < 0.01 e；两套网格目标原子电荷差 < 0.01–0.02 e。
- **差分电荷**：真空区 Δρ ≈ 0；`cdd_planar.dat` 的 ΔQ(z) 积分到真空回到 ≈ 0；0.002 e/Å³ 等值面下核附近无正负交替斑点。
- **ELF**：共价晶体键中点 ≈ 0.9 以上；离子晶体高值集中在阴离子周围；金属间隙 ≈ 0.5 平台。ELFCAR 未乘体积（ASE 读入已除体积）；只解读价层，不能积分。
- **STM**：改 k 网格后图像不变才算收敛；恒高针尖 2–4 Å、恒流 ρ_iso 取 1E-5–1E-3 e/Å³；只比亮点位置，不比绝对值。PBE 单层 MoS2 带隙约 1.6–1.8 eV 可用来对照窗口。
- **COHP**：ICOHP 为负 = 净成键，比较键强看 −ICOHP（对称等价原子应相同）；E_F 处仍明显占据反键态 = 键过度填充；`gaussianSmearingWidth` 与 `SIGMA` 一致。

## 6. 本地完整输入

- `$VASP_TUTORIALS_ROOT/examples/16-bader-co-pt111/`
- `$VASP_TUTORIALS_ROOT/examples/16-cdd-co-pt111/`
- `$VASP_TUTORIALS_ROOT/examples/16-elf-si/`
- `$VASP_TUTORIALS_ROOT/examples/16-stm-mos2/`
- `$VASP_TUTORIALS_ROOT/examples/16-lobster-gaas/`
- `$VASP_TUTORIALS_ROOT/tutorials/16-charge-bonding.md`

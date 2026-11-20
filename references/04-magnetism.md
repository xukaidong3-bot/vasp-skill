# 磁性、DFT+U、交换耦合 J 与 SOC/MAE

> 什么时候读这个文件：要比 FM/AFM/NM 能量、加 Hubbard U、算 J 与 T_C、MAE、或约束磁矩做能量映射时。

## 1. 先判断走哪条路线

| 情形 | 做法 | 关键区别 |
|---|---|---|
| 比三种磁性、定磁序 | `ISPIN = 2` + `MAGMOM`，同参数跑三份 | 只有同 `ENCUT`/k 点/泛函/U 值的总能可比 |
| 强关联氧化物带隙偏小 | DFT+U，`LDAUTYPE = 2` | 先 U ramping；必设 `LMAXMIX` |
| 要 J、T_C / T_N | 算 FM、AFM-I、AFM-II 三个构型 | 只用能量差；巡游磁体需 TB2J |
| 磁各向异性 MAE | `vasp_ncl` + `LSORBIT` + `ICHARG = 11` | **不写 `ISPIN`**、`ISYM = -1`、`NBANDS` ×2 |
| 约束磁矩能量映射 | `I_CONSTRAINED_M` + `M_CONSTR` + `LAMBDA` | 罚能 E_p 必须远小于目标能量差 |

## 2. 标准流程

**1) 磁性金属基座（bcc Fe）**

```text
ENCUT   = 400        # ~1.35 x ENMAX(Fe_pv)=293
EDIFF   = 1E-6
LASPH   = .TRUE.
ISMEAR  = 1          # 磁性金属用 1；SIGMA <= 0.1
SIGMA   = 0.10
ISPIN   = 2
MAGMOM  = 2*3.0      # 2 个 Fe
LORBIT  = 11
LREAL   = .FALSE.
```

KPOINTS `Gamma` `12 12 12`（a=2.83）。查 `grep -A8 "magnetization (x)" OUTCAR | tail -9`（≈2.2 μB/Fe）、`grep "T\*S" OUTCAR` < 1 meV/atom。FM/AFM/NM 三份用 `bash run_fm_afm_nm.sh`（加 `collect` 汇总）：`# AFM : MAGMOM = 3.0 -3.0`、`# NM : 删 ISPIN 与 MAGMOM`；AFM 若 `grep "F=" OSZICAR | tail -1` 的 `mag=` 回到 3.0 量级即塌回 FM。

**2) AFM 超胞 + DFT+U（NiO）** AFM-II 菱方胞：缩放 4.17、矢量 `1 .5 .5 / .5 1 .5 / .5 .5 1`、`Ni O`/`2 2`、坐标 `(0,0,0)(.5,.5,.5)(.25,.25,.25)(.75,.75,.75)`（[111] 交替）；AFM-I 四方胞 a(½,½,0)/a(−½,½,0)/a(0,0,1)（[001] 交替）。

```text
ENCUT     = 520      # 1.3 x ENMAX(O)=400
ISPIN     = 2
MAGMOM    = 2.0 -2.0 2*0.0
LDAU      = .TRUE.
LDAUTYPE  = 2        # Dudarev：只有 U-J 有效
LDAUL     = 2 -1     # Ni d=2；O 不加 U=-1
LDAUU     = 6.2 0.0  # MP 对氧化物中 Ni 的值
LDAUJ     = 0.0 0.0
LDAUPRINT = 1        # 打印占据矩阵
LMAXMIX   = 4        # d=4、f=6；必设
```

KPOINTS `Gamma` `10 10 10`（AFM-I 用 `9 9 6`）；POTCAR 取 `potpaw_PBE/{Ni_pv,O}`；其余同基座。经验 U（MP，eV）：Co 3.32、Cr 3.7、Fe 5.3、Mn 3.9、Mo 4.38、Ni 6.2、V 3.25、W 6.2。加 U 前先 ramping（每步 `ISTART=1`、`LWAVE=.TRUE.`）：

```bash
for U in 0.0 2.0 4.0 6.2; do sed "s/^LDAUU .*/LDAUU = $U 0.0/" INCAR.template > INCAR
  mpirun -np 32 vasp_std > vasp_U$U.out 2>&1; done
grep -A12 "occupancies and eigenvectors" OUTCAR | tail -14   # 本征值应接近 0 或 1
```

**3) 线性响应求 U** 32 原子 2×2×2 AFM-II 超胞，第 1 个 Ni 写成独立物种（POSCAR `Ni Ni O` / `1 15 16`，POTCAR `{Ni_pv,Ni_pv,O}`）。`INCAR.gs` 用 `LWAVE = .TRUE.`，`INCAR.lr` 改：

```text
LDAUTYPE= 3          # 线性响应
LDAUL   = 2 -1 -1
LDAUU   = VPERT 0.00 0.00
LDAUJ   = VPERT 0.00 0.00
ISTART  = 1
#ICHARG = 11         # NSCF 那步由脚本取消注释
```

`VLIST="-0.20 -0.10 -0.05 0.05 0.10 0.20"`，KPOINTS `Gamma` `2 2 2`；`bash run_lr.sh gs|response|fit`。`U = chi^-1 - chi0^-1`（`chi0=dN_NSCF/dV`、`chi=dN_SCF/dV`，用 `python fit_u.py --gs gs/OUTCAR --dirs nscf_* scf_* --plot`）；NiO 参考 `U=5.6 eV`（`chi0=0.49`、`chi=0.13 1/eV`）。

**4) Heisenberg J 与 T_C** 约定 `E = E0 - sum_<ij> J_ij e_i·e_j`（每对算一次，J>0 铁磁）：

```text
NiO（每 Ni）：FM = E0-6J1-3J2、AFM-I = E0+2J1-3J2、AFM-II = E0+3J2
  J1 = (E_AFM-I - E_FM)/8 ,  J2 = (E_AFM-II - E_FM)/6 - J1
bcc Fe：J1_eff = (E_AFM - E_FM)/8（每 Fe）
平均场：k_B T_C^MFA = (1/3) sum_j J_0j
```

`bash run_exchange.sh` / `collect`（=`python exchange_J.py --fm FM_R/OUTCAR --afm2 AFM2/OUTCAR --afm1 AFM1/OUTCAR`）。

**5) SOC 与 MAE** 先共线跑一次（`vasp_std`：`ENCUT=450`、`EDIFF=1E-8`、`ISMEAR=-5`、`ISPIN=2`、`MAGMOM=3.0`、`LMAXMIX=4`、`NBANDS=16`），再用 `vasp_ncl` 跑三方向：

```text
LSORBIT = .TRUE.
# 不要与 ISPIN = 2 同用
SAXIS   = 0 0 1      # _100 用 1 0 0；_111 用 1 1 1
MAGMOM  = 0 0 3.0
ICHARG  = 11         # 力定理要求，读共线 CHGCAR
ISYM    = -1
NBANDS  = 32         # 约共线的 2 倍
LORBMOM = .TRUE.
GGA_COMPAT = .FALSE.
```

`bash run_mae.sh col|soc|collect`：`MAE(dir)=E(dir)-E(001)` 乘 1e6 得 μeV/atom。

**6) 约束磁矩**（非共线，不需 SOC）

```text
LNONCOLLINEAR   = .TRUE.
MAGMOM          = 0 0 2.2  1.1 0 1.905
I_CONSTRAINED_M = 1        # 1 方向、2 方向+大小、4 方向+符号（>=6.4.0）
M_CONSTR        = 0 0 1  0.5 0 0.866
LAMBDA          = 10       # 10 -> 20 -> 50
RWIGS           = 1.20     # < 最近邻距离一半
```

## 3. 关键标签表

| 标签 | 取值 | 何时用 / 注意 |
|---|---|---|
| `ISPIN` / `MAGMOM` | `2`；共线每原子 1 个、非共线 3 个，支持 `2*3.0`、`2*-3.0` | 初猜取预期 1.2–1.5 倍；续算别删；SOC 时不写 `ISPIN` |
| `ISYM` | `0`（关对称）、`-1`（SOC/MAE 必用） | AFM 被对称化成零磁矩时用 0 |
| `LDAUTYPE` | `2` Dudarev（默认）、`3` 线性响应 | type 1 是 Liechtenstein |
| `LDAUL`/`LDAUU`/`LDAUJ` | d=`2`、f=`3`、不加 U=`-1`；`LMAXMIX` d=`4`、f=`6` | 与 POTCAR 物种顺序对应；加 U 必设 `LMAXMIX`，`LDAUPRINT=1` 打印占据矩阵 |
| `SAXIS` / `GGA_COMPAT` | `0 0 1`（默认）、`1 0 0`、`1 1 1`；`.FALSE.` | `LHFCALC=.TRUE.` 强制 `ISYM=3`，只支持默认 `SAXIS`；SOC/MAE 必设 `GGA_COMPAT=.FALSE.` |
| `LORBMOM` | `.TRUE.` | 输出轨道磁矩（笛卡尔系，与自旋系不同） |

## 4. 易错点

1. **`MAGMOM` 个数写错** → `Error reading item 'MAGMOM' from file INCAR.`。原因：6.5.0 起 `ISPIN=2` 与 `LNONCOLLINEAR`/`LSORBIT` 同现按非共线规则查个数。处置：删掉 `ISPIN` 行（≥6.5.1 由 `LNONCOLLINEAR` 优先、`ISPIN` 回退 1）；`vasp_std` 跑 SOC 会直接退出（需不带 `-DNGXhalf`/`-DNGZhalf` 编译），换 `vasp_ncl`。
2. **不设 `LMAXMIX` 就加 U** → `ICHARG=11` 的能带/自洽 DOS 与 SCF 对不上。原因：CHGCAR 只按默认 `LMAXMIX=2` 写。处置：SCF 设 `LMAXMIX = 4`（d）/`6`（f）重跑；续算别删 `MAGMOM`（参与对称性判断，删掉会把 AFM 磁矩对称化成零）。
3. **不同 U 的总能直接比较** → 反应能/相稳定性全错。原因：Dudarev 总能随 U−J 单调变。处置：所有体系（含参考态）同 U；GGA 与 GGA+U 混用要加 MP2020 修正。
4. **AFM 跑成 FM 或落进亚稳态** → `mag=` 回到饱和值，或同结构两次差几十 meV（d 占据序不同）。处置：核对 `MAGMOM` 顺序、试 `ISYM = 0`、加大初猜；DFT+U 用多起点 + ramping + `LDAUPRINT = 1`，必要时上 OMC 补丁。
5. **MAE 随 k 点跳动甚至变号** → 能量差只有 μeV 量级。处置：`ISYM = -1`、`EDIFF = 1E-8`、k 网格 16³→24³→32³ 验、`NBANDS` 加倍。
6. **版本坑** → 负 `MAGMOM` 时 `NBANDS` 偏少（issue 51）；≤6.6.1 杂化+SOC 的 PAW 单中心交换项对部分元素有误（issue 116）。处置：用 6.6.1 起并显式写 `NBANDS`。
7. **磁矩振荡不收敛** → 处置：`AMIX_MAG = 0.8`、`BMIX_MAG = 0.0001`（`BMIX = 0` 会让某些版本崩）；或先 `ISPIN = 1` 再转 `ISPIN = 2` + `ICHARG = 1`；或用 `NUPDOWN` 卡住总磁矩再放开。

## 5. 算完了怎么验收

- **磁矩与能量序**：`grep -A8 "magnetization (x)" OUTCAR | tail -9` 的逐原子 `tot`。AFM 总磁矩须 0、子格 ±；bcc Fe ≈2.2 μB/Fe，NiO 的 Ni ≈±1.7 μB、O 近 0；子格变 0 或全同号 → 磁序设错。NiO 应 AFM 最低、bcc Fe 应 FM 最低；差 < 10 meV 先确认 `T*S` < 1 meV/atom。
- **占据矩阵与带隙**：`grep -A12 "occupancies and eigenvectors" OUTCAR | tail -14` 本征值应接近 0 或 1（分数占据 → 亚稳态）；Ni²⁺ d⁸ 应多数自旋 5 个 d 全满、少数 t2h 满 eg 空。NiO 带隙 PBE < 1 eV、PBE+U(6.2) ≈3 eV（实验 ~4 eV）；PBE+U 仍 0.5 eV → U 没生效。
- **线性响应与 J**：`|chi0|` 明显大于 `|chi|`、残差 ≲1E-3、U 与文献同量级（NiO ~5.6 eV）。AFM-II 必须最低、`J2 < 0`、`|J1| << |J2|`；`FM_T` 与 `FM_R` 每磁性原子差 > 1 meV → k 点/ENCUT 不收敛，J 与 T_N 不可信。MFA 高估 T_N 数十个百分点（NiO 实验 523 K）。
- **MAE 与约束磁矩**：`grep "free  energy   TOTEN" soc_*/OUTCAR` 核对；`[001]` 与 `[100]` 的能量差即噪声本底，真 MAE 须显著大于它（bcc Fe 的 `E[111]-E[001]` ≈1 μeV/atom，轨道矩 0.0x μB）；`grep "NBANDS=" col/OUTCAR soc_001/OUTCAR` 确认已翻倍。约束磁矩时 `M_int` 应贴合 `M_CONSTR`，`E_p` 远小于要分辨的能量差。

## 6. 本地完整输入

- `$VASP_TUTORIALS_ROOT/tutorials/07-magnetism-u-soc.md`
- `$VASP_TUTORIALS_ROOT/examples/07-fe-bcc-fm/`
- `$VASP_TUTORIALS_ROOT/examples/07-nio-afm-dftu/`
- `$VASP_TUTORIALS_ROOT/examples/07-linear-response-u/`
- `$VASP_TUTORIALS_ROOT/examples/07-fe-bcc-soc-mae/`

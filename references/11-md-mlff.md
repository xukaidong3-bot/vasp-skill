# 从头算分子动力学（AIMD）与机器学习力场（MLFF）

> 什么时候读这个文件：用户要跑有限温动力学——熔化/液态、扩散、RDF/MSD、VDOS、NVT/NpT、约束 MD，或用 VASP 内置 MLFF（train→refit→run）加速 AIMD 时。静态性质不读这里。

## 1. 先判断走哪条路线

| 情形 | 做法 | 关键区别 |
|---|---|---|
| 几 ps、几十到一两百原子 | `IBRION = 0` + `MDALGO = 2/4`（NVT），或 `MDALGO = 3` + `ISIF = 3`（NpT） | 每步一次 SCF，最贵；`POTIM` 必须显式给 |
| 几十到几百 ps、数百原子、要大统计量 | MLFF：train → refit → run | 训练期在线做 DFT，生产期纯预测（refit 后约快 20–100 倍） |
| NVT 且要 MSD/扩散系数 | 必须避开 `MDALGO = 0` | 旧 Nosé–Hoover 总把坐标折回盒子，MSD 无法计算 |
| 液体等无长程有序体系要 NpT | `MDALGO = 3` + `ISIF = 3` + `ICONST` 约束盒子形状 | 完全自由的 NpT 会把盒子压扁成一层 |
| 只做预测 | GRACE（≥6.6，实验性，须 `mpirun -np 1`）或外部 MLIP | 精度依体系，须自行验证 |

## 2. 标准流程

### A. AIMD：Si64 的 NVT（3000 K 熔化 → 1800 K 生产）

1. **建模**：2×2×2 金刚石惯用胞 64 原子，盒长 10.511 Å（近熔点液态 Si 约 2.57 g/cm³）；`KPOINTS` 只 Γ → `vasp_gam`。2. **第 1 段熔化**（1500 步 × 2 fs = 3 ps）：

```
PREC=Normal  ENCUT=250  EDIFF=1E-5  ALGO=Fast  LREAL=Auto
ISMEAR=0  SIGMA=0.10  NELMIN=4  MAXMIX=40  ISYM=0
IBRION=0  MDALGO=2  SMASS=1.0  ISIF=2
TEBEG=3000  TEEND=3000  POTIM=2.0  NSW=1500  NBLOCK=1
RANDOM_SEED=248489752 0 0   LWAVE=.FALSE.  LCHARG=.FALSE.
```

检查 `grep "T=" OSZICAR | tail -5`：温度应围绕 3000 K 涨落。平均电子步数 `grep -cE "^(DAV|RMM):" OSZICAR` ÷ MD 步数应为个位数。

3. **第 2 段生产**：`cp CONTCAR POSCAR`，只改 `TEBEG=TEEND=1800`、`NSW=5000`（10 ps），丢弃前 1–2 ps 作平衡。

### B. AIMD：32 H₂O 的 NpT（Langevin + Parrinello–Rahman）

4. **INCAR 关键行**（电子步同 A，另加色散，`ENCUT` 必须提高）：

```
MDALGO=3  ISIF=3  TEBEG=300  TEEND=300  POTIM=0.5  NSW=6000
LANGEVIN_GAMMA=10.0 10.0    # 每个元素种类一个值
LANGEVIN_GAMMA_L=10.0  PMASS=1000  PSTRESS=0.0  IVDW=12
ENCUT=550                   # ≥1.3×ENMAX(O)=400：变胞时基组固定
```

5. **`ICONST` 保持立方盒**（`LA`/`LR` 是原始坐标，`S` 是复合坐标，末列 `0` = 约束）：

```
LA 1 2 0
LA 1 3 0
LA 2 3 0
LR 1 0
LR 2 0
LR 3 0
S  1  0  0  0  0  0  0
S  0  1  0  0  0  0  0
S  0  0  1  0  0  0  0
S  0  0  0  1 -1  0  0
S  0  0  0  1  0 -1  0
S  0  0  0  0  1 -1  0
```

检查 `grep "volume of cell" OUTCAR`：密度在平衡值附近涨落（32 分子约 ±0.03 g/cm³），非单调漂移。

### C. MLFF：Si64 的 train → refit → run

6. **在线训练**（DFT 段在续训时**一字不改**）：

```
PREC=Normal  ENCUT=350  EDIFF=1E-6  ALGO=Fast  LREAL=Auto
ISMEAR=0  SIGMA=0.05  NELM=100  ISYM=0        # 不要设 MAXMIX
IBRION=0  MDALGO=3  ISIF=3  TEBEG=200  TEEND=1000  POTIM=2.0  NSW=10000
LANGEVIN_GAMMA=10.0  LANGEVIN_GAMMA_L=10.0  PMASS=1000
ML_LMLFF=.TRUE.  ML_MODE=train
```

检查 `grep STATUS ML_LOGFILE | grep -v "#" | awk '{print $3}' | sort | uniq -c`：开头几乎每步 `critical`/`learning`，随后 DFT 间隔迅速变稀。
7. **续训**：`cp ../1_train/ML_ABN ML_AB`，照抄 `INCAR_train`，只改温度与步数。
8. **refit**（生产前必做）：`cp 1_train/ML_ABN 2_refit/ML_AB`，INCAR 只留 `ML_LMLFF=.TRUE.` + `ML_MODE=refit`（自动设 `ML_LFAST=.TRUE.`、`NSW=1`、SVD 拟合）。检查 `grep ERR ML_LOGFILE | grep -v "#" | tail -1`。
9. **生产运行**：`cp 2_refit/ML_FFN 3_run/ML_FF`、`cp 1_train/CONTCAR 3_run/POSCAR`，MD 段同训练但 `TEBEG=TEEND=800`、`NSW=100000`、`LANGEVIN_GAMMA=1.0`，加 `ML_LMLFF=.TRUE.`、`ML_MODE=run`、`ML_OUTBLOCK=10`、`ML_ESTBLOCK=50`。检查 ML_LOGFILE 溢出因子应接近 0（→1 表示在外推）。

## 3. 关键标签表

| 标签 | 取值 | 何时用 / 注意 |
|---|---|---|
| `POTIM` | 含 H 0.5；O/C/N 1–1.5；Si 等 2–3（fs） | `IBRION = 0` 时**无默认值**，不写必崩 |
| `MDALGO` | 0 旧 N–H（折回坐标，不能算 MSD）；1 Andersen；2 Nosé–Hoover；3 Langevin（NVT/NpT/NpH）；4 NHC；5 CSVR（≥6.4.3）；13 多恒温器 | 全部要求编译带 `-Dtbdyn` |
| `ISIF` | 0/2 定容（2 算应力）；3 变胞（**仅** `MDALGO = 3`） | MD 中不支持 4 和 7 |
| `LANGEVIN_GAMMA` | 平衡 5–10、生产 ≤1（ps⁻¹） | 个数必须 = POSCAR 元素种类数；摩擦大的人为减慢扩散 |
| `LANGEVIN_GAMMA_L` / `PMASS` | 10.0 / 1000（amu） | 不设摩擦默认 0（不同系综）；盒长加倍 `PMASS` 也加倍 |
| `ML_MODE` | train/select/refit/refitbayesian/run/delta/none（≥6.4.0，取代 6.3 的 `ML_ISTART`） | `train` 只要目录有 `ML_AB` 就**总是**续训，从零训练须删掉 |
| `ML_RCUT1`/`ML_RCUT2` | 8 / 5 Å（默认） | 描述符截断，refit 时调；液体可试 `ML_RCUT2 = 4` |
| `ML_CTIFOR` | 自动更新；300–500 K 约 0.02、1000–2000 K 约 0.06 | 太低每步做 DFT，太高很早停止学习 |
| `ML_ESTBLOCK` | 20–100（旧名 `ML_IERR`） | 只能在 `run` 模式自由设置 |

## 4. 易错点

1. **现象**：MD 一开始就崩溃。→ 原因：`IBRION = 0` 时 `POTIM` 没有默认值。→ 处置：显式写 `POTIM`。
2. **现象**：MSD 很快饱和像固体。→ 原因：`MDALGO = 0` 折回坐标（官方明确 MSD 无法计算），或温度/时长不足仍是晶体。→ 处置：改 `MDALGO = 2/4`；先 3000 K 熔化再降温；分析前按最小像展开 `Δs = s(t+1) − s(t) − round(Δs)`。
3. **现象**：NpT 中液体盒子被压扁、严重剪切，或体积持续漂移。→ 原因：液体无抗剪切能力；变胞时基组的 Pulay 应力；`PMASS` 不当。→ 处置：`ICONST` 约束形状；`ENCUT` 比定容高 ≥30% 并分段重启；试 `PMASS = 100–1000`。
4. **现象**：`MDALGO` 不生效或标签不被识别。→ 原因：编译无 `-Dtbdyn`。→ 处置：按官方模板重编译。
5. **现象**：续训后力场反而更差。→ 原因：忘 `cp ML_ABN ML_AB`（丢上一轮新数据），或改了 DFT 参数/POTCAR。→ 处置：改名复制；DFT 段与 POTCAR 完全不变（一个力场只拟合一个势能面）。
6. **现象**：`ML_MODE = run` 报找不到 `ML_FF` 或类型不匹配。→ 原因：忘 `cp ML_FFN ML_FF`。→ 处置：复制 refit 的 `ML_FFN`。
7. **现象**：训练中 SCF 频繁不收敛、电子步异常多。→ 原因：设了 `MAXMIX > 0`。→ 处置：删掉 `MAXMIX`。
8. **现象**：很早就全 `accurate`（几乎不做 DFT），或一直每步都学习。→ 原因：前者 `ML_CTIFOR` 沿用了液体阈值；后者变形/相变、温度过高、DFT 设置不一致。→ 处置：降阈值回默认；分阶段训练（体相→表面→分子+表面）。
9. **现象**：报错要求增大 `ML_MCONF` 或参考构型达 `ML_MB`；train 模式设 `ML_ESTBLOCK ≠ 1` 报错。→ 处置：增大上限（从已有 `ML_AB` 续训不丢数据）；`ML_ESTBLOCK` 只在 `run` 设。
10. **现象**：生产运行中止、溢出因子 > 0.9。→ 原因：6.6 起超 0.9 直接终止。→ 处置：把该温度/压力/相加入续训后重新 refit。

## 5. 算完了怎么验收

- **OSZICAR MD 行**：`T=` 围绕 `TEBEG` 涨落无漂移（64 原子相对涨落 √(2/(3N)) ≈ 10% 正常）；`E=`（守恒量）长时间不漂移，漂移即 `POTIM` 太大或 `EDIFF` 太松。
- **RDF**：第一峰后无晶体的尖锐二、三峰；配位数取到第一极小，截断半径 ≤ 最短晶胞高度一半（液态 Si 第一峰约 2.4–2.5 Å，水 O–O 约 2.7–2.8 Å）。
- **扩散**：只拟合 MSD 线性区（跳过 MSD ∝ t² 的弹道区与过长尾部），D = MSD/(6t)；1 Å²/fs = 0.1 cm²/s，液体量级 10⁻⁴ cm²/s；固体 MSD 饱和正常。
- **MLFF 精度**：`ERR` 行 300–1000 K 下能量 < ~1 meV/atom、力 30–100 meV/Å 可接受；`BEEF` 一般小于 `ERR`，只是相对指标；`LCONF` 随学习增长（液体约需 2000–4000，简单固体 500–1000）。训练误差小 ≠ 可靠：抽轨迹帧用与训练**完全相同**的参数做 DFT 单点（py4vasp `MLFFErrorAnalysis`），换 `RANDOM_SEED` 跑多条轨迹，时间平均应一致。
- **算错信号**：温度失控上升、守恒能单调漂移、`ERR` 跳升；压力看 `grep -E "external pressure|total pressure" OUTCAR`。

## 6. 本地完整输入

- `$VASP_TUTORIALS_ROOT/examples/15-aimd-si-liquid-nvt/`（`INCAR.melt`、`INCAR.prod`）
- `$VASP_TUTORIALS_ROOT/examples/15-aimd-water-npt/`（`INCAR`、`ICONST`）
- `$VASP_TUTORIALS_ROOT/examples/15-mlff-si-train-run/`（`INCAR_train`、`INCAR_refit`、`INCAR_run`）
- `$VASP_TUTORIALS_ROOT/scripts/md_analysis.py`
- `$VASP_TUTORIALS_ROOT/tutorials/15-md-mlff.md`

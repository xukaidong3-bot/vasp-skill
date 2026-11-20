# 结构优化：离子位置与晶胞（02-structure）

> 什么时候读这个文件：用户要「优化结构 / 算晶格常数 / 算体模量 / 二维材料只放面内 / 分子放真空盒 / 体积偏小或真空层塌了」，或判断一次优化是否真收敛、要不要续算。

## 1. 先判断走哪条路线

| 情形 | 做法 | 关键区别 |
|---|---|---|
| 体相晶体，晶格常数未知 | `ISIF = 3` 变胞 → `cp CONTCAR POSCAR` 再优化 1–2 轮到体积不变 | 变胞必须 `ENCUT ≥ 1.3×max(ENMAX)`，否则 Pulay 应力让体积偏小 |
| 已知晶格（表面/吸附/缺陷超胞） | `ISIF = 2` 只动原子；固定底层用 POSCAR `Selective dynamics` | 固定体积时 Pulay 应力近似各向同性常数，影响小 |
| 二维材料单层 | `ISIF = 3` + `LATTICE_CONSTRAINTS = .TRUE. .TRUE. .FALSE.` | slab 直接 `ISIF = 3` 会压塌真空层；与 `IBRION = 1/2` 组合需 VASP ≥ 6.4.3（推荐 ≥ 6.5.1） |
| 孤立分子/团簇 | 12–15 Å 盒子 + `Gamma 1 1 1` + `vasp_gam` + `ISIF = 2` | 盒子人为，绝不优化晶胞；不能用 `ISMEAR = -5` |
| 要 B0、B0′ 或想避开 Pulay | E–V 扫描 `V/Vref = 0.94–1.06`（7 点）→ Birch–Murnaghan 拟合 | 每体积重建基组，Pulay 误差被平均掉 |

## 2. 标准流程

1. **第一阶段：手搭结构（受力 > 1 eV/Å）**：`PREC = Normal`、`EDIFF = 1E-5`、`EDIFFG = -0.05`、`IBRION = 2`、`NSW = 100`；跑完 `cp CONTCAR POSCAR`。
2. **第二阶段：生产优化**
   ```text
   PREC = Accurate ; ENCUT = 400   # >= 1.3 x ENMAX(Si ~245 eV)：变胞/应力必须
   EDIFF = 1E-6 ; ISMEAR = 0 ; SIGMA = 0.05 ; LREAL = .FALSE.
   IBRION = 2 ; ISIF = 3 ; NSW = 60
   EDIFFG = -0.005                 # 负值 = 按力；ISIF>=3 时应力也受检
   LWAVE = .FALSE. ; LCHARG = .FALSE.
   ```
   检查：`grep "reached required accuracy" OUTCAR` 有输出。
3. **消 Pulay 应力：再跑一轮**（不要设 `ISTART = 2`，那会把旧基组连同 Pulay 应力一起带过来）：
   ```bash
   mkdir -p run1 && cp INCAR KPOINTS POTCAR POSCAR OUTCAR OSZICAR CONTCAR run1/
   cp CONTCAR POSCAR && sbatch job.slurm
   ```
   判据：第二轮 1–3 个离子步内结束；两轮体积差 < 0.1–0.2%，能量差 < 1 meV/atom；`external pressure` 约 1 kB 内。
4. **二维材料只放面内**（1H-MoS₂，c = 20 Å 真空，KPOINTS `12 12 1`）：
   ```text
   ENCUT = 450 ; LASPH = .TRUE.    # >= 1.3 x max(ENMAX) of Mo_pv/S；d 元素建议 LASPH
   IBRION = 2 ; ISIF = 3 ; NSW = 100 ; EDIFFG = -0.01
   LATTICE_CONSTRAINTS = .TRUE. .TRUE. .FALSE.   # a1,a2 放开；a3 固定
   ```
   检查：`head -5 CONTCAR` 第三行仍是 `0 0 20`；`grep "in kB" OUTCAR | tail -1` 中 XX/YY/XY ≈ 0，ZZ 非零正常。
5. **分子放盒子**（CO，15 Å，用 `vasp_gam`）：`ENCUT = 520`、`ISMEAR = 0`、`SIGMA = 0.01`、`ISPIN = 1`（O2/自由基用 `ISPIN = 2` + `MAGMOM`）、`IBRION = 2`、`ISIF = 2`、`NSW = 50`、`EDIFFG = -0.01`。CO 键长 PBE 约 1.13–1.14 Å。
6. **E–V 扫描 + BM 拟合**：按体积缩放 POSCAR 第 2 行（体积 s → 长度乘 s^(1/3)），每点单点算（`NSW = 0`、`IBRION = -1`，Si 可用 `ISMEAR = -5`）：
   ```bash
   for d in V_*; do echo "$(awk '/volume of cell/{v=$NF} /energy\(sigma->0\)/{e=$NF} \
   END{print v, e}' $d/OUTCAR)"; done | sort -n > eos.dat
   python fit_bm.py eos.dat --lattice fcc --plot
   ```
   判据：V0 落在扫描区间内（最好中间三分之一）；残差 meV 量级；BM/Murnaghan/Vinet 的 V0 差 < 0.1%；与 `ISIF = 3` 体积差 < 0.5%。

## 3. 关键标签表

| 标签 | 取值 | 何时用 / 注意 |
|---|---|---|
| `IBRION` | 2 CG（首选，对 POTIM 不敏感）/ 1 RMM-DIIS（近极小时快，需精确受力，配 `NELMIN = 4–8`）/ 3 阻尼 MD（`SMASS ≈ 0.4`，必须给 POTIM） | `IBRION = 3` + `POTIM = 0` + `IOPT` 即 VTST 优化器（1 LBFGS / 7 FIRE） |
| `ISIF` | 2 只动原子 / 3 全放开 / 4 固定体积 / 5 只变形状 / 6 只变晶胞 / 8 原子+体积（≥ 6.4.1） | 杂化泛函默认 `ISIF = 0` 且不算应力，必须显式写 |
| `EDIFFG` | 负值按力：常规 −0.01~−0.02，声子/弹性前 −1E-3；默认 = `EDIFF×10`（正值，偏松） | `ISIF ≥ 3` 时应力和力都受检，但仍要自己看压强 |
| `EDIFF` | 1E-6 生产 / 1E-7–1E-8 声子与弹性 | 太松会让 CG 线搜索失败（`ZBRENT`） |
| `NSW` | 50–200 | 跑满 NSW 不报错，但**不会**打印 `reached required accuracy` |
| `POTIM` | 默认 0.5；`IBRION = 1` 振荡时降到 0.1–0.3 | `IBRION = 2` 会打印建议的 `trialstep`，乘当前 POTIM 即合适步长 |
| `LATTICE_CONSTRAINTS` | 三个逻辑值对应 a1 a2 a3 | 非正交晶胞即把应力张量对应行列置零；旧版本**静默忽略** |
| `PSTRESS` | kB（10 GPa = 100 kB） | 加压优化配 `ISIF = 3`；比较相稳定性看 H = E + PV |

## 4. 易错点

1. **`ISIF = 3` 优化出晶格明显偏小** → Pulay 应力（固定平面波数目、基组不完备，几乎总让体积偏小）。→ `ENCUT ≥ 1.3×max(ENMAX)`；CONTCAR → POSCAR 重优化到体积不变；或改用 E–V 拟合。
2. **二维材料真空层被压缩** → 对 slab 直接用了 `ISIF = 3`，真空区没有应力撑着。→ 加 `LATTICE_CONSTRAINTS = .TRUE. .TRUE. .FALSE.`；旧版本该标签在 `IBRION = 1/2` 下不读，只能改 E–a 扫描（固定 c，a 缩放 ±4% × 7 点，每点 `ISIF = 2` 后拟合）。
3. **写了 `LATTICE_CONSTRAINTS` 但 c 轴照旧在变** → 版本不足：< 6.5.0 有「约束在收敛判据检查前未施加」的 bug，与 `ISIF = 4/5` 组合需 ≥ 6.5.1（known issue 34）。→ 升级并在 OUTCAR 核对晶格变化。
4. **杂化泛函优化时晶格完全不动** → `LHFCALC = .TRUE.` 下 `ISIF` 默认 0、不算应力。→ 显式 `ISIF = 3`。
5. **出现 `reached required accuracy` 但晶胞仍偏** → `ISIF ≥ 3` 时负 EDIFFG 也检查应力分量（内部单位，不直观）。→ 另外单独看 `external pressure`。
6. **`ZBRENT: fatal error in bracketing`（或 `bracketing interval incorrect`）** → CG 线搜索失败：离极小太远、受力噪声大（EDIFF 太松）、能量面很平。→ `cp CONTCAR POSCAR` 续算；`EDIFF` 收到 1E-6 或更小；减小 `POTIM`；改 `IBRION = 1` 或 VTST FIRE（`IBRION=3, POTIM=0, IOPT=7`）。
7. **`WARNING: Sub-Space-Matrix is not hermitian in DAV` / `Error EDDDAV: Call to ZHEGV failed`** → 原子靠太近或步长过大。→ 查 CONTCAR 原子间距、减小 POTIM、删旧 WAVECAR 重来。
8. **分子用 `ISMEAR = -5` 报 `Tetrahedron method fails (number of k-points < 4)`** → 四面体法要 ≥ 4 个 k 点。→ 孤立体系 `ISMEAR = 0`、`SIGMA = 0.01–0.05`。
9. **E–V 曲线锯齿状** → 基组/FFT 网格随体积不连续变化、k 点太稀。→ 提高 ENCUT、`PREC = Accurate`、加密 k 点；6.6.0 + 自旋极化 r2SCAN 另有受力错误（issue 94），升 6.6.1 或打 `pot.F.patch`。

## 5. 算完了怎么验收

- `grep "reached required accuracy" OUTCAR`：出现 `reached required accuracy - stopping structural energy minimisation` 才算收敛；没有 = 跑满 NSW 或中断。
- 自由自由度上的最大力 < |EDIFFG|（固定住的原子/方向要排除，否则误判不收敛）。
- `grep "external pressure" OUTCAR | tail -1`：|P| 约 1–2 kB 以内；`grep "volume of cell" OUTCAR` 最后两步体积不再变。
- 晶格从 `head -5 CONTCAR` 读：fcc 原胞 `|a1| = a/√2`，故 `a = √2·|a1|`。Si 的 PBE 值约 5.46–5.47 Å（实验 5.431 Å，GGA 高估约 1%）属正常；Si 的 B0 PBE 约 88–90 GPa（实验 98–99 GPa）。

## 6. 本地完整输入

- `$VASP_TUTORIALS_ROOT/tutorials/04-relax-eos.md`
- `$VASP_TUTORIALS_ROOT/examples/04-relax-si-isif3/INCAR`
- `$VASP_TUTORIALS_ROOT/examples/04-relax-mos2-2d/INCAR`
- `$VASP_TUTORIALS_ROOT/examples/04-relax-co-molecule/INCAR`
- `$VASP_TUTORIALS_ROOT/examples/04-eos-si/fit_bm.py`

# 声子：phonopy 有限位移 / DFPT / NAC / 虚频诊断

> 什么时候读这个文件：用户要声子色散、声子 DOS、热力学量（F/S/Cv）、LO–TO 劈裂、结构动力学稳定性，或问"这个虚频是真的吗"。只管晶格振动；弹性与介电/光学去 `10-elastic-optics.md`。

## 1. 先判断走哪条路线

| 情形 | 做法 | 关键区别 |
|---|---|---|
| 默认（含 vdW/杂化/DFT+U/meta-GGA） | phonopy 有限位移：`phonopy-init -d` → 每个 `POSCAR-00N` 单点（`IBRION=-1`）→ `phonopy-init -f` | 泛函不限，各位移独立可并行；位移幅度需自选 |
| 纯 LDA/GGA，想少花机时 | DFPT：`SPOSCAR` 一次 `IBRION=8` → `phonopy-init --fc vasprun.xml` | 仅 LDA/GGA，不支持 `IVDW`，必须无 `NCORE/NPAR`；不用选幅度 |
| 极性材料（MgO、ZnO、钙钛矿） | 任一路线 **+ 单胞 `LEPSILON=.TRUE.` → `phonopy-vasp-born` → `BORN`** | 不加 NAC 就漏 Γ 点 LO–TO 劈裂，光学支频率整体偏低 |
| 不想装 phonopy / 大体系要拆分续算 | VASP 原生 `IBRION=6`（+`LPHON_DISPERSION`） | 需 ≥6.3.2；`CHECKPOINT_FD` 续算需 ≥6.6 且编译开 HDF5 |
| 二维材料 | 有限位移，超胞只在面内扩 | `--dim 5 5 1`、c ≥ 15–20 Å、KPOINTS z 取 1 |

`IBRION=5–8` 本身只给**超胞的 Γ 点**频率；整条色散只能靠 phonopy 或 `LPHON_DISPERSION` 做傅里叶插值。

## 2. 标准流程（有限位移）

**1) 高精度弛豫**——声子极依赖"是否在能量极小点"：

```
PREC = Accurate
ENCUT = 400        # ≥1.3×ENMAX，与后续超胞必须完全一致
EDIFF = 1E-8
EDIFFG = -1E-3
ISMEAR = 0
SIGMA = 0.01
LREAL = .FALSE.
IBRION = 2
ISIF = 3
NSW = 100
```
查：`grep "reached required accuracy" OUTCAR`；`grep "external pressure" OUTCAR|tail -1` 应 |P|<0.5 kB。`cp CONTCAR POSCAR` 再跑一轮（PAW 基组不随变胞更新）。

**2) 生成位移超胞**：

```bash
phonopy-init -d --dim 2 2 2 -c POSCAR   # ≥4.x；2/3.x 用 phonopy -d
ls  # SPOSCAR POSCAR-001 phonopy_disp.yaml（后处理要用，别删）
```
每方向 ≥10–12 Å；软模体系再大。

**3) 每个位移做单点（绝不弛豫）**：

```
PREC = Accurate
ENCUT = 400        # 同弛豫
EDIFF = 1E-8
NELMIN = 5
ISMEAR = 0
SIGMA = 0.01
LREAL = .FALSE.    # 否则力噪声会造虚频
ISYM = 0
IBRION = -1
NSW = 0
NCORE = 4
```
KPOINTS 用等效密度：原胞 12×12×12 → 2×2×2 超胞 6×6×6，4×4×4 超胞 3×3×3（示例为 8 原子惯用胞扩 2×2×2 ≈10.9 Å，取 Γ 中心 3×3×3）。

```bash
bash run_disp.sh prepare|submit|check|collect
```
查：每个 `OSZICAR` 只有 1 个离子步（`grep -c 'F=' = 1`）；SCF 步数 ≠ `NELM`；`vasprun.xml` 以 `</modeling>` 结尾。

**4) 收集力**：`phonopy-init -f disp-*/vasprun.xml` → `FORCE_SETS`。查 phonopy 打印的 `Max drift of force constants`，到 0.01 eV/Å² 量级说明力不准；`disp-*` 字典序须与 `POSCAR-001/002…` 对应。

**5) 色散与热力学**（`band.conf`）：

```
ATOM_NAME = Si
DIM = 2 2 2
PRIMITIVE_AXES = 0 1/2 1/2  1/2 0 1/2  1/2 1/2 0
BAND = 0 0 0  1/2 0 1/2  1/2 1/4 3/4  3/8 3/8 3/4  0 0 0  1/2 1/2 1/2
BAND_LABELS = $\Gamma$ X W K $\Gamma$ L
BAND_POINTS = 101
BAND_CONNECTION = .TRUE.
```

```bash
phonopy --config band.conf -p -s
phonopy-bandplot --gnuplot band.yaml > band.dat
```
`mesh.conf` 写 `MESH = 24 24 24`、`GAMMA_CENTER = .TRUE.`、`TPROP = .TRUE.`、`TMIN/TMAX/TSTEP = 0/1000/10`；投影 DOS 加 `PDOS = 1, 2`。输出 `total_dos.dat`、`thermal_properties.yaml`（**按原胞**计，kJ/mol 与 J/K/mol，1 kJ/mol ≈ 0.010364 eV）。

**6) DFPT**：`phonopy-init -d --dim 2 2 2 -c POSCAR-unitcell` 只留 `SPOSCAR`；超胞 `IBRION = 8`（**不写 `NCORE/NPAR`**）；`grep -c 'name="hessian"' vasprun.xml` 须非零；`phonopy-init --fc vasprun.xml` → `FORCE_CONSTANTS`；`band.conf` 加 `FORCE_CONSTANTS = READ`。

**7) NAC**：在**单胞**（非超胞）上 `LEPSILON=.TRUE.`+`LRPA=.FALSE.`+`IBRION=-1`，k 点约为超胞等效密度 2 倍；然后
```bash
phonopy-vasp-born --pa "0 1/2 1/2 1/2 0 1/2 1/2 1/2 0" > ../BORN
```
`band.conf` 写 `NAC = .TRUE.`（≥4.x 有 BORN 自动开启，`--nonac` 关）。要定方向用 `QPOINTS = 0 0 0` + `Q_DIRECTION = 1 0 0`。

**8) QHA / phono3py（入门）**：QHA 取 **7–11 个体积、−5%…+5%**（至少 5），各体积 `ISIF=2`（非立方 `ISIF=4`）+ 一次声子 → `thermal_properties.yaml-{1..11}` 与 `e-v.dat` → `phonopy-qha -p -s e-v.dat thermal_properties.yaml-{...}`（新版用 `run_qha`）。三声子：`phono3py-init -d --dim 2 2 2 --pa F -c POSCAR-unitcell` → 各 `disp-0000n` 单点 → `phono3py-init --cf3 disp-{00001..00755}/vasprun.xml` → `phono3py --mesh 11 11 11 --br`。

## 3. 关键标签表

| 标签 | 取值 | 何时用 / 注意 |
|---|---|---|
| `IBRION` | `-1` / `2` / `6` / `8` | 位移超胞只能 `-1`+`NSW=0`；`6`=对称性约化有限差分（`5`=全部位移）；`8`=DFPT |
| `POTIM` | `0.015` | `IBRION=5/6` 的位移步长（Å）；过大超出谐区，VASP 5.1+ 会自动重置为此值 |
| `NFREE` | `2` | `1` 强烈不推荐；`3/5` 等非法值被重置为 2 |
| `EDIFF` | `1E-8` | 有限差分/DFPT 必须；放宽直接变虚频 |
| `LREAL` | `.FALSE.` | 全流程禁用，实空间投影的力误差不可接受 |
| `LEPSILON`/`LCALCEPS` | `.TRUE.` | 前者 DFPT（仅 LDA/GGA），后者有限电场（杂化唯一选择） |
| `PRIMITIVE_AXES` | 显式 9 个数 | v4 默认 `AUTO`，升级后 q 路径会变；显式写才可复现 |
| `LPHON_DISPERSION` | `.TRUE.` | ≥6.3.2，配 `QPOINTS`（KPOINTS 格式）插值色散 + `PHON_NWRITE` |
| `PHON_DOS`/`PHON_SIGMA`/`PHON_NEDOS` | `1`（高斯）/`2`（四面体） | 声子 DOS，配 `QPOINTS` 均匀网格 |
| `LPHON_POLAR`+`PHON_DIELECTRIC`+`PHON_BORN_CHARGES` | `.TRUE.` + 9 数 / 每原子 9 数 | VASP 原生 NAC，数值取自单胞 `LEPSILON`/`LCALCEPS` |
| `CHECKPOINT_FD` | `RESET`/`CONTINUE`/`NONE`/`PREPARE`/`SINGLE` | ≥6.6、开 HDF5、**仅 `IBRION=6`**；续算日志出现 `Continuing from previous run` |

## 4. 易错点

1. **整条支大虚频** → 把超胞弛豫了（`NSW>0`/`IBRION≠-1`）或参考结构没到极小 → 改回 `IBRION=-1`+`NSW=0`，先用 `ISIF=3`、`EDIFFG=-1E-3` 重做结构再生成超胞。
2. **`PEAD routines do not work for NCORE` / 提示 `remove the tag NPAR from the INCAR file`** → 线性响应不能能带并行 → 删 `NCORE`、`NPAR`，要并行用 `KPAR`。
3. **DFPT 力常数不对 / `--fc` 失败** → 计算没跑完，或体系带 vdW → `grep -c 'name="hessian"'` 须非零；DFPT 不支持 `IVDW`，有 vdW 走有限位移。
4. **`"band.conf" could not be read as a phonopy.yaml-like file. If "band.conf" is a configuration file, give it by --config.`** → 4.x 不再接受 conf 作位置参数 → `phonopy --config band.conf -p -s`。
5. **`phonopy -d` / `-f` 报迁移错误** → 4.x 前处理拆到 `phonopy-init` → 用 `phonopy-init -d` / `-f` / `--fc`。
6. **`phonopy-vasp-born --pa auto` 报 `could not convert string to float`** → 只吃数值矩阵 → 写全 `--pa "0 1/2 1/2 1/2 0 1/2 1/2 1/2 0"`，且与 `PRIMITIVE_AXES` 一致。
7. **升级后谱线形状/q 路径变了** → 默认 `PRIMITIVE_AXES=AUTO` → conf 里显式写矩阵。
8. **Γ 附近极小虚频（<0.3 THz）** → 力精度/ASR/超胞不足 → `EDIFF=1E-8`、`PREC=Accurate`、`LREAL=.FALSE.`、加大超胞。二维 ZA 支 Γ 邻域的小虚频多为旋转不变性未满足，通常不代表真实不稳定。
9. **只在高对称点（X/M/R）明显虚频** → 真实软模（倾转、CDW） → 造能把该 q 折到 Γ 的超胞沿本征矢位移，或用 IRC（`IBRION=40`）找新极小。
10. **`thermal_properties` 报 `GammaAcousticWarning`** → Γ 近零声学模进入求和 → 4.7.1+ 设 `EXCLUDE_GAMMA_ACOUSTIC = .TRUE.`。

## 5. 算完了怎么验收

- `band.yaml`/`band.pdf`：Γ 三支声学模 ≈0（偏离 >0.1 THz 说明力常数不满足 ASR）；Si 光学模 Γ 约 15 THz，MgO 的 TO/LO 约 12/21 THz（量级参考）。
- 虚频扫描：`bash phonopy_workflow.sh imag --imag-tol 0.1`（`--imag-tol` 默认 0.1 THz，有虚频退出码 1）；或直接看 `band.yaml` 负频率。
- 收敛：超胞从 2×2×2 加到 3×3×3 后关键频率变化 <0.1 THz；`ENCUT` 扰动使 Γ 光学模移动 <0.05 THz；`MESH` 24³→32³ 使 F/Cv 变化 <0.1%。
- NAC：Γ 附近出现 LO–TO 劈裂（谱线在两侧不连续是正常的）；所有原子 Z* 之和应≈0，偏差大 = k 点不够。
- 热力学：`thermal_properties.yaml` 按原胞计；高温极限 Cv 应趋近 3Nk_B（N = 原胞原子数），否则 MESH/PDOS 设错。

## 6. 本地完整输入

- `$VASP_TUTORIALS_ROOT/examples/12-si-phonopy-fd/`
- `$VASP_TUTORIALS_ROOT/examples/12-mgo-dfpt-nac/`
- `$VASP_TUTORIALS_ROOT/examples/12-graphene-phonon/`
- `$VASP_TUTORIALS_ROOT/scripts/phonopy_workflow.sh`
- `$VASP_TUTORIALS_ROOT/tutorials/12-phonon.md`

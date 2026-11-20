# 四个输入文件、续算与收敛测试（01-inputs）

> 什么时候读这个文件：要写或审 INCAR/POSCAR/KPOINTS/POTCAR；怀疑「标签写了却没生效」；要从上次算续算（ISTART/ICHARG）；要给新体系定 ENCUT 与 k 点并做收敛测试。

## 1. 先判断走哪条路线

| 情形 | 做法 | 关键区别 |
|---|---|---|
| 新体系第一个自洽/单点 | `ISTART = 0`、`ICHARG = 2`，先删旧 WAVECAR/CHGCAR | 读到不匹配的旧文件是 `BRMIX` 报错主因 |
| 优化中断后续算 | `cp CONTCAR POSCAR`；要省首步 SCF 就留 WAVECAR + `ISTART = 1` | 变胞续算**不要** `ISTART = 2`：基组不变，Pulay 应力会被原样带走 |
| 能带/DOS 非自洽 | 自洽段 `LCHARG = .TRUE.`，第二段 `ICHARG = 11` + Line-mode | Line-mode 不能自洽；≥6.3 可用 `KPOINTS_OPT` 一步完成 |
| 系列小改动（改 U、改 MAGMOM） | `ICHARG = 1` 读上一个 CHGCAR | 官方建议，显著省电子步 |
| 定参数 | ENCUT → k 网格 → ISMEAR/SIGMA | 一次只改一个；测试时 `ISTART = 0`、`EDIFF = 1E-7` |

## 2. 标准流程

1. **净化格式**（头号杀手是 Tab 与 Windows CRLF）：
   ```bash
   dos2unix INCAR POSCAR KPOINTS
   grep -n -E '\^I|\^M' <(cat -A INCAR)   # 应无输出
   ```
   检查：有 Tab/CR 会出现 `Error reading item 'MAGMOM' from file INCAR.` 或标签被静默忽略。
2. **按 POSCAR 第 6 行顺序拼 POTCAR 并核对**：
   ```bash
   cat $VASP_PP_PATH/potpaw_PBE/{Ga_d,As}/POTCAR > POTCAR
   grep TITEL POTCAR          # 应依次 Ga_d、As
   sed -n '6p' POSCAR         # 应写 Ga As
   grep -E "TITEL|ZVAL|ENMAX" POTCAR
   ```
   检查：顺序一致；记下 `ENMAX`（定 ENCUT）与 `ZVAL`（定 NELECT）。
3. **最小 INCAR**（半导体单点）：
   ```text
   ISTART = 0 ; ICHARG = 2
   PREC   = Accurate
   ENCUT  = 320        # 1.3 x ENMAX(PAW_PBE Si = 245.3 eV)
   LREAL  = .FALSE. ; ALGO = Normal
   EDIFF  = 1E-6 ; NELM = 100
   ISMEAR = 0          # 默认 1 是给金属的
   SIGMA  = 0.05
   IBRION = -1 ; NSW = 0
   LWAVE  = .FALSE. ; LCHARG = .TRUE.
   ```
   检查：`ISMEAR`/`SIGMA` 必须显式写（默认 `ISMEAR = 1`、`SIGMA = 0.2`）。
4. **跑前自检 + 运行**：
   ```bash
   for f in INCAR POSCAR POTCAR; do [ -s $f ] || echo "ERROR: $f missing"; done
   [ -s KPOINTS ] || grep -qi "^ *KSPACING" INCAR || echo "WARN: 将用默认 KSPACING = 0.5"
   export OMP_NUM_THREADS=1; ulimit -s unlimited
   mpirun -np 4 vasp_std > vasp.out 2>&1
   ```
5. **核对 VASP 真读到了什么**（拼错标签不报错，只静默忽略）：
   ```bash
   grep -A 3 "Startparameter for this run" OUTCAR
   grep -E "ENCUT|PREC|ISMEAR|SIGMA|EDIFF |NELM|LREAL|ISPIN" OUTCAR | head -20
   ```
   检查：`ISMEAR` 显示 1 而你想写 0 → 拼成了 `ISMAER` 之类。
6. **ENCUT 收敛测试**（Si：250 275 300 325 350 375 400 450 500 eV，即 1.02–2.04×ENMAX）：
   ```text
   ISTART = 0 ; PREC = Accurate
   ENCUT  = __ENCUT__
   LREAL  = .FALSE.
   EDIFF  = 1E-7        # 比生产严一档，否则 SCF 噪声掩盖趋势
   ISMEAR = 0 ; SIGMA = 0.05
   NSW = 0 ; IBRION = -1 ; ISIF = 2   # 单点也算应力，观察 Pulay 应力
   LWAVE = .FALSE. ; LCHARG = .FALSE.
   ```
   判据：ΔE < 1 meV/atom（以最高 ENCUT 为参考）、ΔP < 1 kB；压强通常比能量收敛慢。
7. **k 点收敛**：固定上一步 ENCUT，扫 Γ 中心 N×N×N（Si 原胞 N = 2…16 对应 NKPTS = 3/8/16/29/47/72/104/145）；或不写 KPOINTS 而在 INCAR 写 `KSPACING`（0.5→0.15 Å⁻¹ 对应 Si 5³→14³）。判据同为 ΔE < 1 meV/atom。

## 3. 关键标签表

| 标签 | 取值 | 何时用 / 注意 |
|---|---|---|
| `ISTART` | 0 新算 / 1 读 WAVECAR 重建基组 / 2 基组不变 / 3 另读 TMPCAR | WAVECAR 缺失或版本不符时自动退回 0，ICHARG 一并变 2 |
| `ICHARG` | 2 原子叠加 / 1 读 CHGCAR 并外推 / 11 固定电荷非自洽 / 12 Harris–Foulkes | 11 时 ENCUT/PREC/ISPIN/POTCAR 必须与自洽段一致；DFT+U 还要 `LMAXMIX = 4`(d)/`6`(f) |
| `LCHARG` / `LWAVE` | `.TRUE.`/`.FALSE.` | 续算与能带才需要；每个离子步都写会占满磁盘 |
| `KSPACING` | 0.15–0.35 Å⁻¹（金属取小） | **只在没有 KPOINTS 文件时生效**，默认 0.5；网格会「跳格」，EOS 比较改用显式网格 |
| `LORBIT` | 11 | 写 PROCAR 与分轨道 DOSCAR；要 DOS/能带才开 |
| `EFERMI` | `MIDGAP` | VASP ≥ 6.4；不确定有无带隙时让费米能级落在带隙中间 |
| `MAGMOM` | 个数 = NIONS（`vasp_ncl` 为 3×NIONS） | 个数不符报 `Error reading item 'MAGMOM' from file INCAR.`；可写 `2*3.0 4*0` |
| `NWRITE` / `LH5` | 0–3（默认 2）/ `.TRUE.` | NWRITE 控 OUTCAR 详细度（长 MD 用 0/1）；LH5 改输出 vaspout.h5/vaspwave.h5，此时不再写 WAVECAR/CHGCAR |

## 4. 易错点

1. **结果离谱但不报错** → POSCAR 第 6 行元素符号只是注释，原子身份由 POTCAR 顺序 + 第 7 行数目定。→ `grep TITEL POTCAR` 对照 `sed -n '6p' POSCAR`。
2. **标签「没起作用」** → 拼错（`ISMAER`）、Tab、CRLF 都被静默忽略。→ 看 OUTCAR 的 `Startparameter for this run` 回显，`dos2unix` 后重传。
3. **k 点很少、能量/费米能级偏差大** → 既没 KPOINTS 也没写 KSPACING，VASP 静默用默认 `KSPACING = 0.5`。→ 补 KPOINTS 文件。
4. **`VERY BAD NEWS! internal error in subroutine IBZKPT: Reciprocal lattice and k-lattice belong to different class of lattices`** → fcc/六方晶格用了偶数 Monkhorst–Pack 网格。→ 改 Γ 中心网格（fcc、六方、二维、四面体法都只用 Γ 中心）。
5. **`Tetrahedron method fails`（`number of k-points < 4`）** → `ISMEAR = -5` 但 k 点少于 4 个。→ 改 `ISMEAR = 0`、`SIGMA = 0.05`。
6. **半导体带隙/占据数异常** → 忘写 ISMEAR，用上默认 Methfessel–Paxton。→ 显式 `ISMEAR = 0`（或 `-5`，需 ≥ 4 k 点）。
7. **`WARNING: dimensions on CHGCAR file are different`** → `ICHARG = 11` 的 ENCUT/PREC 与自洽段不一致。→ 两段用同一套 ENCUT/PREC/网格。
8. **`BRMIX: very serious problems the old and the new charge density differ`** → 读入不匹配的 CHGCAR/WAVECAR。→ 删旧文件重算；必要时调 `SYMPREC` 或 `ISYM = 0`。
9. **`internal error in subroutine PRICEL` / `SGRCON`** → POSCAR 有效数字不足（官方建议 ≥ 7 位）导致对称性判断不一致。→ 对称化后重写，或调 `SYMPREC`。
10. **`The distance between some ions is very small`** → Direct/Cartesian 混用或缩放因子写错。→ VESTA/ASE 打开 POSCAR 目检第 2 行与坐标模式行。

## 5. 算完了怎么验收

- `OUTCAR` 出现 `General timing and accounting informations for this job` = 正常结束；没有就别读结果。
- `tail -3 OSZICAR`：最后离子步电子步数应远小于 `NELM`；**恰好等于 NELM 说明该步没收敛**。
- `grep "energy(sigma->0)" OUTCAR | tail -1` 与 `grep -E "free +energy +TOTEN" OUTCAR | tail -1` 之差即 T*S；`grep "EENTRO" OUTCAR | tail -1` 的每原子 `|T*S|` 应 < 1 meV，超了减小 SIGMA 并加密 k 点。
- `grep "NELECT" OUTCAR` 应等于 ΣZVAL；`grep "NKPTS" OUTCAR` 与预期网格一致（Si 8×8×8 Γ → 29）。
- `grep "total drift" OUTCAR | tail -1`：力漂移须远小于关心的力；>0.1 eV/Å 说明 ENCUT/FFT 网格不足。
- 收敛判据：ΔE < 1 meV/atom、ΔP < 1 kB、力变化 < 5–10 meV/Å；能量随 ENCUT 不单调 = EDIFF 太松或读了旧 WAVECAR。

## 6. 本地完整输入

- `$VASP_TUTORIALS_ROOT/tutorials/02-input-output-files.md`
- `$VASP_TUTORIALS_ROOT/tutorials/03-convergence.md`
- `$VASP_TUTORIALS_ROOT/examples/02-input-files-si/Si_scf/INCAR`
- `$VASP_TUTORIALS_ROOT/examples/02-input-files-si/INCAR_annotated_template`
- `$VASP_TUTORIALS_ROOT/examples/03-convergence-si-al/encut_test/INCAR.template`

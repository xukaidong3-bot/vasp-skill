# VASP Agent Skill

这是一个给 AI agent 使用的 VASP 参数决策、标准流程与排障 skill。`SKILL.md` 是决策层和路由索引；`references/` 保存按专题按需加载的细节，帮助 agent 选择参数、组织计算并验收结果。

维护者：`xukaidong3-bot`。

## 目录结构

- `SKILL.md`：skill 元数据、路由表、通用骨架与硬约束。
- `references/00-build-env.md` … `14-postprocess-errors.md`：15 个专题 reference。
- `ATTRIBUTION.md`、`LICENSE`：来源、许可与软件边界声明。

## 使用

1. 将本 skill 放入 agent 的 skills 目录并按需加载 `SKILL.md`。
2. 根据 §1 路由表只读取与任务相关的一个 `references/*.md`。
3. 若需要外部讲义中的示例，先设置 `VASP_TUTORIALS_ROOT=/path/to/DFT-Tutorials/VASP`。

## 覆盖范围

编译与并行环境；输入文件与收敛测试；结构优化与 EOS；能带、DOS/PDOS 与高级泛函；磁性、DFT+U、SOC/MAE；vdW 与二维材料；slab、表面能与功函数；吸附、催化与自由能；NEB、CI-NEB、Dimer 与过渡态；声子与热力学；弹性常数、介电与光学；AIMD 与 MLFF；电荷、ELF、STM 与 COHP；Wannier、有效质量、反折叠、费米面、拓扑、缺陷及后处理报错。

## 不包含

本仓库不包含 VASP 本体、源代码、POTCAR/赝势或计算输出，也不代替 VASP 官方手册、许可证或软件供应商的支持。VASP 是商业软件，使用者须自行取得合法授权。

---
name: pr-review-lab
description: Use when reviewing SwanLab Python SDK PRs, running gh pr diff or gh pr review, checking public API/type stubs/run lifecycle/settings/API/CLI/transport/concurrency/integrations/media/probe/protobuf/packaging/Python compatibility changes, or gating SDK merge readiness. Reviews locally first and only posts to GitHub after the user confirms.
---

# PR 审核 (SwanLab)

## 概述

按照 SwanLab Python SDK 的仓库约定审核 Pull Request，并用 `gh` CLI 提交结构化结论。重点确认变更不引入公共 API、运行生命周期、数据上传、配置认证、跨平台兼容、安全、性能、构建和可维护性回归。

技术栈：**Python 3.9-3.14**、**uv**、**Hatch**、**pytest**、**Ruff**、**basedpyright**、**Pydantic**、**requests**、**Protobuf**。

适用：SwanLab SDK PR 合并前审核、`gh pr diff` / `gh pr review`、检查 SDK/API/CLI/集成/Proto/依赖/测试变更。

不适用：SwanLab-Cloud 前端 PR、其他仓库 PR、无 PR 上下文的单文件编辑、只要求解释代码。

**CRITICAL — 审核 PR 前检查文件列表。当包含 `core/**`、`protos/**` 或 `scripts/generate_protos.py` 时，MUST 先用 Read 工具读取 [`../pr-review-core/SKILL.md`](../pr-review-core/SKILL.md) 及其 [`references/domain-contracts.md`](../pr-review-core/references/domain-contracts.md)，按其中的 Go core 领域契约审核对应文件。审核结论需同时覆盖两侧范围。**

## 审核流程

审核分为两个阶段，**默认只执行阶段一**。未获得用户明确确认前，不得向 GitHub 提交任何 review、comment 或审批状态。

### 阶段一：审查（默认）

1. 先确认 PR 意图、base 分支、head 分支、文件列表和 Actions 状态，再阅读实现。

```bash
gh pr view <PR>
gh pr view <PR> --json files,additions,deletions,baseRefName,headRefName
gh pr diff <PR>
gh pr checks <PR>
```

2. 把 head 分支拉到本地隔离工作区再审核，不要污染用户当前工作区：

```bash
git fetch origin pull/<PR>/head:pr-<PR>
git worktree add ../swanlab-pr-<PR> pr-<PR>
```

无法建立 worktree 时可用 `gh pr checkout <PR>`，但必须先确认当前工作区干净，并在结束后切回原分支。

3. 根据文件列表判断适用领域，并按“范围参考”读取相关契约。所有应用代码都要检查横切质量；涉及凭据、网络、文件或用户数据时额外检查安全与隐私。
4. 对关键发现检查 head 分支中的完整文件和调用关系，不只依据 diff 片段。
5. 对照仓库现有实现、测试和公开契约定级。个人风格偏好不能作为阻塞问题。
6. 复杂 PR 可并行从“仓库约定与影响范围”和“正确性、安全、性能与验证”两个视角审核；小 PR 无需拆分。
7. 在对话中输出完整审查意见草稿：结论、逐条 finding（含文件与行号）、验证结果。**到此停止**，明确告知用户尚未发布，并询问是否发布。

### 阶段二：发布（需用户显式确认）

仅当用户在阶段一之后明确要求发布（如“发布”“提交 review”“post”）时才执行“提交审核”一节的命令。

- 用户只是要求审查、看意见或讨论时，一律停在阶段一。
- 用户要求修改草稿时，更新草稿后重新征求确认，不要顺带发布。
- 发布前把最终结论和 `event` 类型（`APPROVE` / `REQUEST_CHANGES` / `COMMENT`）复述给用户确认，避免审批状态与用户意图不符。
- 发布后回报 review 链接，并清理临时 worktree 和本地分支：

```bash
git worktree remove ../swanlab-pr-<PR>
git branch -D pr-<PR>
```

## 独立审核原则

- 每次审核以当前 base、head、diff、实际文件和 PR Actions 为准。
- 历史 review、评论和审批状态只作为线索，finding 必须在当前代码上重新验证。
- 已修复的问题不重复提出；仍存在的问题按当前影响重新定级。
- PR 作者和外部 reviewer 使用相同标准。平台不允许提交 `APPROVE` 或 `REQUEST_CHANGES` 时，使用 `COMMENT`，正文仍写明实际结论。

## 范围参考

详细规则位于 [`references/domain-contracts.md`](references/domain-contracts.md)。先读取 PR 文件列表，再按下表加载适用章节。路径一律写完整前缀；同一文件命中多行时，读取所有命中章节。

| 变更路径或内容                                                                                                                                                                                       | 读取章节                                 |
| ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------- |
| `swanlab/__init__.py`、`swanlab/__init__.pyi`、`swanlab/sdk/__init__.py`、`swanlab/sdk/typings/**`、`swanlab/sdk/protocol/**`、`swanlab/exceptions.py`                                               | 公共 API 与类型兼容                      |
| `swanlab/sdk/cmd/**`、`swanlab/sdk/internal/run/**`、`swanlab/sdk/internal/context/**`、`swanlab/sdk/internal/bus/**`、`swanlab/sdk/internal/impl.py`                                                | Run 生命周期与运行时组件                 |
| `swanlab/sdk/internal/core_python/**`（含 `transport/**`、`store/**`、`metrics/**`、`watcher/**`、`sync.py`）                                                                                        | Record、Core 与传输                      |
| `swanlab/sdk/internal/settings/**`、`swanlab/sdk/internal/pkg/client/**`、`swanlab/sdk/internal/pkg/nrc/**`、`swanlab/sdk/internal/core_python/client/**`、`swanlab/sdk/internal/core_python/api/**` | 配置、认证与 HTTP 客户端                 |
| `swanlab/api/**`、`swanlab/cli/**`                                                                                                                                                                   | Public API 与 CLI                        |
| `swanlab/sdk/internal/run/transforms/**`、文件保存和媒体处理                                                                                                                                         | 媒体、转换与文件                         |
| `swanlab/sdk/internal/probe_python/**`（含 `hardware_vendor/**`、`monitor/**`、`environment/**`）、`swanlab/vendor/**`                                                                               | 硬件与环境探测                           |
| `swanlab/integration/**`、`swanlab/converter/**`、`swanlab/plugin/**`                                                                                                                                | 集成、转换器与插件                       |
| `swanlab/sdk/internal/pkg/**`（`safe/`、`fs/`、`console/`、`fork/`、`timer/`、`executor/`、`scope/`、`adapter/` 等）、`swanlab/utils/**`                                                             | 共享工具与内部 pkg                       |
| `protos/**`、`swanlab/proto/**`、`core/proto/**`                                                                                                                                                     | Protobuf 与生成代码                      |
| `pyproject.toml`、`uv.lock`、`Makefile`、`swanlab/package.json`、`.github/workflows/**`                                                                                                              | 构建、依赖与发布                         |
| `tests/**`                                                                                                                                                                                           | 测试与验证                               |
| `swanlab/deprecated/**`                                                                                                                                                                              | 公共 API 与类型兼容 + 横切质量与代码规范 |
| 所有 Python 代码、配置和依赖                                                                                                                                                                         | 横切质量与代码规范                       |
| 凭据、网络、日志、文件路径、用户数据                                                                                                                                                                 | 安全与隐私                               |

上表未命中的路径仍需按“横切质量与代码规范”审核，并在 review 中说明该文件按通用标准处理。

## 验证要求

优先运行与改动直接相关的测试，再根据范围执行完整门禁：

```bash
uv run pytest <相关测试路径>
uv run ruff check .
uv run basedpyright
uv run pytest
```

涉及打包、依赖或版本逻辑时补充：

```bash
uv build
```

涉及 `.proto` 时，确认源定义与 Python/Go 生成代码同步；具备完整工具链时在干净或临时工作区运行：

```bash
make proto
```

CI 覆盖 Ubuntu、Windows、macOS 和 Python 3.9-3.14。无法在本地覆盖的环境以 PR Actions 为准；Action 未运行或失败时必须如实记录，不能视为通过。

## 严重性

| 前缀          | 含义                 | 示例                                                         |
| ------------- | -------------------- | ------------------------------------------------------------ |
| **阻塞:**     | 阻塞合并             | 凭据泄露、数据丢失、公共 API 严重破坏、发布产物不可用        |
| **必须修改:** | 合并前必须处理       | 类型契约不同步、生命周期错误、竞态、兼容性回归、必要测试缺失 |
| **建议:**     | 值得改进但通常不阻塞 | 局部可维护性、命名、适度简化                                 |
| **细节:**     | 次要问题             | 小范围风格或文档问题                                         |
| **说明:**     | 信息记录，不要求行动 | 验证限制、后续风险或上下文                                   |

每条 finding 必须使用显式前缀，并说明问题、影响和建议。不要用无前缀文本表达必须修改的问题。

## 提交审核（仅阶段二）

本节所有命令都属于发布动作，**只有在用户明确确认后才执行**。阶段一只在对话中给出草稿。

优先把可行动 finding 提交为 PR diff 上的 inline review comment。顶层 review body 只放总体结论、审核范围、验证结果和无法挂到 diff 行的问题，不重复 inline finding。

### 选择提交方式

- 无 finding 或只需总评：使用 `gh pr review --approve/--comment --body-file`。
- 有可定位 finding：通过 Reviews API 在一次 review 中提交 `comments[]`。
- 有阻塞或必须修改 finding：使用 `REQUEST_CHANGES`。
- 只有建议、细节或说明：使用 `COMMENT`。
- 无可行动 finding 且验证通过：使用 `APPROVE`。

### 定位 inline comment

```bash
gh pr diff <PR> --patch --color=never
```

- 新增或修改后的代码使用 `side: "RIGHT"` 和新文件行号。
- 删除导致的问题使用 `side: "LEFT"` 和旧文件行号。
- 只评论 diff 中存在的行；跨文件问题或目标行不在 diff 时写入顶层正文。
- 每条 inline comment 只包含一个 finding。

### 提交 inline review

```bash
gh api --method POST repos/{owner}/{repo}/pulls/<PR>/reviews --input /tmp/pr-review.json
```

Review JSON 使用以下结构：

```json
{
  "event": "REQUEST_CHANGES",
  "body": "## 审核结论: 要求修改\n\n## 本次 PR 解决的问题\n<一句话概述>\n\n## 审核范围\n<范围清单>\n\n## 未挂行发现\n<没有则写：无>\n\n## 验证\n<命令与 Actions 状态>",
  "comments": [
    {
      "path": "swanlab/path/file.py",
      "line": 42,
      "side": "RIGHT",
      "body": "**必须修改:** <问题>\n\n影响：<影响>\n\n建议：<建议>"
    }
  ]
}
```

### 顶层正文结构

```markdown
## 审核结论: 通过 / 要求修改 / 评论

## 本次 PR 解决的问题

<一句话概述>

## 审核范围

- [x] 公共 API 与类型兼容: <详情>
- [x] Run 生命周期与运行时: <详情>
- [x] Record、Core 与传输: <详情>
- [x] 配置、认证、API 与 CLI: <详情>
- [x] 媒体、集成与插件: <详情>
- [x] Protobuf、构建与依赖: <详情>
- [x] 安全、性能与跨平台兼容: <详情>
- [x] 测试与验证: <详情>

## 未挂行发现

<没有则写：无>

## 验证

- [ ] 相关测试
- [ ] uv run ruff check .
- [ ] uv run basedpyright
- [ ] uv run pytest
- [ ] PR Actions 已通过，或已记录未运行/失败状态
```

没有 inline finding 时使用正文文件提交：

```bash
gh pr review <PR> --request-changes --body-file /tmp/pr-review-body.md
gh pr review <PR> --approve --body-file /tmp/pr-review-body.md
gh pr review <PR> --comment --body-file /tmp/pr-review-body.md
```

## 完成验证

阶段一（审查）：

- [ ] 当前 PR 意图、base、head 和文件范围已确认
- [ ] head 分支已拉到本地隔离工作区，关键发现在完整文件上验证
- [ ] 所有适用领域和横切质量已审核
- [ ] 阻塞和必须修改 finding 已解决或明确记录
- [ ] 本地验证结果和 PR Actions 状态已记录
- [ ] 审查意见草稿已输出，并已明确告知用户尚未发布

阶段二（发布，需用户确认）：

- [ ] 用户已明确要求发布，且已确认 `event` 类型
- [ ] 可定位 finding 已提交为 inline comments
- [ ] 最终 review 已通过 Reviews API 或 `gh pr review` 提交
- [ ] 临时 worktree 和本地分支已清理

## 参考

- [pr-review-core](../pr-review-core/SKILL.md) — 当 PR 包含 `core/**` 变更时必读其 [`references/domain-contracts.md`](../pr-review-core/references/domain-contracts.md)

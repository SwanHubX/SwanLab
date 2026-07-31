---
name: pr-review-lab
description: Use when reviewing SwanLab Python SDK PRs, running gh pr diff or gh pr review, checking public API/type stubs/run lifecycle/settings/API/CLI/transport/concurrency/integrations/media/protobuf/packaging/Python compatibility changes, or gating SDK merge readiness.
---

# PR 审核 (SwanLab)

## 概述

按照 SwanLab Python SDK 的仓库约定审核 Pull Request，并用 `gh` CLI 提交结构化结论。重点确认变更不引入公共 API、运行生命周期、数据上传、配置认证、跨平台兼容、安全、性能、构建和可维护性回归。

技术栈：**Python 3.9-3.14**、**uv**、**Hatch**、**pytest**、**Ruff**、**basedpyright**、**Pydantic**、**requests**、**Protobuf**。

适用：SwanLab SDK PR 合并前审核、`gh pr diff` / `gh pr review`、检查 SDK/API/CLI/集成/Proto/依赖/测试变更。

不适用：SwanLab-Cloud 前端 PR、其他仓库 PR、无 PR 上下文的单文件编辑、只要求解释代码。

## 审核流程

1. 先确认 PR 意图、base 分支、head 分支、文件列表和 Actions 状态，再阅读实现。

```bash
gh pr view <PR>
gh pr view <PR> --json files,additions,deletions,baseRefName,headRefName
gh pr diff <PR>
gh pr checks <PR>
```

2. 根据文件列表判断适用领域，并按“范围参考”读取相关契约。所有应用代码都要检查横切质量；涉及凭据、网络、文件或用户数据时额外检查安全与隐私。
3. 对关键发现检查 head 分支中的完整文件和调用关系，不只依据 diff 片段。
4. 对照仓库现有实现、测试和公开契约定级。个人风格偏好不能作为阻塞问题。
5. 复杂 PR 可并行从“仓库约定与影响范围”和“正确性、安全、性能与验证”两个视角审核；小 PR 无需拆分。

## 独立审核原则

- 每次审核以当前 base、head、diff、实际文件和 PR Actions 为准。
- 历史 review、评论和审批状态只作为线索，finding 必须在当前代码上重新验证。
- 已修复的问题不重复提出；仍存在的问题按当前影响重新定级。
- PR 作者和外部 reviewer 使用相同标准。平台不允许提交 `APPROVE` 或 `REQUEST_CHANGES` 时，使用 `COMMENT`，正文仍写明实际结论。

## 范围参考

详细规则位于 [`references/domain-contracts.md`](references/domain-contracts.md)。先读取 PR 文件列表，再按下表加载适用章节。

| 变更路径或内容 | 读取章节 |
|---|---|
| `swanlab/__init__.py`、`swanlab/__init__.pyi`、`swanlab/sdk/__init__.py`、`swanlab/sdk/typings/**` | 公共 API 与类型兼容 |
| `swanlab/sdk/cmd/**`、`swanlab/sdk/internal/run/**`、`context/**`、`bus/**` | Run 生命周期与运行时组件 |
| `swanlab/sdk/internal/core_python/**`、`transport/**`、`store/**` | Record、Core 与传输 |
| `swanlab/sdk/internal/settings/**`、`pkg/client/**`、`pkg/nrc/**` | 配置、认证与 HTTP 客户端 |
| `swanlab/api/**`、`swanlab/cli/**` | Public API 与 CLI |
| `swanlab/sdk/internal/run/transforms/**`、文件保存和媒体处理 | 媒体、转换与文件 |
| `swanlab/integration/**`、`swanlab/converter/**`、`swanlab/plugin/**` | 集成、转换器与插件 |
| `protos/**`、`swanlab/proto/**`、`core/proto/**` | Protobuf 与生成代码 |
| `pyproject.toml`、`uv.lock`、`Makefile`、`.github/workflows/**` | 构建、依赖与发布 |
| `tests/**` | 测试与验证 |
| 所有 Python 代码、配置和依赖 | 横切质量与代码规范 |
| 凭据、网络、日志、文件路径、用户数据 | 安全与隐私 |

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

| 前缀 | 含义 | 示例 |
|---|---|---|
| **阻塞:** | 阻塞合并 | 凭据泄露、数据丢失、公共 API 严重破坏、发布产物不可用 |
| **必须修改:** | 合并前必须处理 | 类型契约不同步、生命周期错误、竞态、兼容性回归、必要测试缺失 |
| **建议:** | 值得改进但通常不阻塞 | 局部可维护性、命名、适度简化 |
| **细节:** | 次要问题 | 小范围风格或文档问题 |
| **说明:** | 信息记录，不要求行动 | 验证限制、后续风险或上下文 |

每条 finding 必须使用显式前缀，并说明问题、影响和建议。不要用无前缀文本表达必须修改的问题。

## 提交审核

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

- [ ] 当前 PR 意图、base、head 和文件范围已确认
- [ ] 所有适用领域和横切质量已审核
- [ ] 阻塞和必须修改 finding 已解决或明确记录
- [ ] 本地验证结果和 PR Actions 状态已记录
- [ ] 可定位 finding 已提交为 inline comments
- [ ] 最终 review 已通过 Reviews API 或 `gh pr review` 提交

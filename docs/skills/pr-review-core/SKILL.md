---
name: pr-review-core
description: Use when reviewing SwanLab Go core PRs under core/**, running gh pr diff or gh pr review for Go changes, checking cross-platform build tags/goroutine lifecycle/parent process monitoring/console logging/gRPC-protobuf contracts/golangci-lint/go test -race/Go module changes, or gating Go core merge readiness. Reviews locally first and only posts to GitHub after the user confirms.
---

# PR 审核 (SwanLab Go core)

## 概述

按照 SwanLab Go core 的仓库约定审核 Pull Request，并用 `gh` CLI 提交结构化结论。重点确认变更不引入跨平台兼容、进程生命周期、goroutine 安全、gRPC/protobuf 契约、诊断日志和工具链回归。

技术栈：**Go 1.26.5**（`core/.go-version` 与 `core/go.mod`）、**golangci-lint v2**、**gRPC**、**protobuf**、**golang.org/x/sys**。

适用：`core/**` Go core PR 合并前审核、`gh pr diff` / `gh pr review`、检查 Go core 进程管理/控制台日志/proto/模块/工具链变更。

不适用：SwanLab Python SDK PR（使用 `pr-review-lab`）、SwanLab-Cloud 前端 PR、无 PR 上下文的单文件编辑、只要求解释代码。

**CRITICAL — 审核 PR 前检查文件列表。当包含 `swanlab/**`、`protos/**`、`scripts/generate_protos.py`、`pyproject.toml` 或 `Makefile` 时，MUST 先用 Read 工具读取 [`../pr-review-lab/SKILL.md`](../pr-review-lab/SKILL.md) 及其 [`references/domain-contracts.md`](../pr-review-lab/references/domain-contracts.md)，按其中的 Python SDK 领域契约审核对应文件。审核结论需同时覆盖两侧范围。**

## 审核流程

审核分为两个阶段，**默认只执行阶段一**。未获得用户明确确认前，不得向 GitHub 提交任何 review、comment 或审批状态。两阶段结构、独立审核原则、严重性前缀和提交审核格式（inline comment 定位、Review JSON、顶层正文模板）均同 [`pr-review-lab`](../pr-review-lab/SKILL.md)；以下仅列出 Go core 专属内容。

### 阶段一：审查（默认）

1. 先确认 PR 意图、base 分支、head 分支、文件列表和 Actions 状态，再阅读实现。

```bash
gh pr view <PR>
gh pr view <PR> --json files,additions,deletions,baseRefName,headRefName
gh pr diff <PR>
gh pr checks <PR>
```

2. 直接在当前工作区 checkout 到 PR head 审核，不创建单独的 worktree 和分支。checkout 前确认当前工作区干净；存在未提交改动时停下来询问用户如何处理，不要继续：

```bash
git status
gh pr checkout <PR>
```

3. 根据文件列表判断适用领域，按下方"范围参考"读取 [`references/domain-contracts.md`](references/domain-contracts.md) 中对应章节。所有 Go 代码都要检查"Go core 横切质量与代码规范"；涉及凭据、网络、日志或文件路径时额外检查安全（同 `pr-review-lab` 的"安全与隐私"）。
4. 对关键发现检查 head 分支中的完整文件和调用关系，不只依据 diff 片段。
5. 对照 `core/` 现有实现、测试和平台实现定级。个人风格偏好不能作为阻塞问题。
6. 在对话中输出完整审查意见草稿：结论、逐条 finding（含文件与行号）、验证结果。**到此停止**，明确告知用户尚未发布，并询问是否发布。

### 阶段二：发布（需用户显式确认）

仅当用户在阶段一之后明确要求发布时才执行。发布前把最终结论和 `event` 类型（`APPROVE` / `REQUEST_CHANGES` / `COMMENT`）复述给用户确认。inline comment 提交方式、Review JSON 结构和顶层正文模板见 [`pr-review-lab` "提交审核"](../pr-review-lab/SKILL.md) 一节——提交时把审核范围复选框替换为下方的 Go core 版本。

发布后回报 review 链接，并切回原分支：

```bash
git checkout <原分支>
```

## 范围参考

详细规则位于 [`references/domain-contracts.md`](references/domain-contracts.md)。先读取 PR 文件列表，再按下表加载适用章节。路径一律写完整前缀；同一文件命中多行时，读取所有命中章节。

| 变更路径或内容                                                                                                         | 读取章节                         |
| ---------------------------------------------------------------------------------------------------------------------- | -------------------------------- |
| `core/internal/pkg/process/**`                                                                                         | 平台与进程生命周期               |
| `core/internal/pkg/console/**`（含 `log/`）                                                                            | 控制台与诊断日志                 |
| `core/proto/**`、`protos/**`、`scripts/generate_protos.py`                                                             | Protobuf 与生成代码              |
| `core/go.mod`、`core/go.sum`、`core/.golangci.yml`、`core/.go-version`、`core/Makefile`、`Makefile`（`core-*` target） | Go 模块、工具链与构建            |
| 所有 `core/**` Go 代码、配置和依赖                                                                                     | Go core 横切质量与代码规范       |
| 凭据、网络、日志、文件路径、用户数据                                                                                   | 安全与隐私（同 `pr-review-lab`） |

上表未命中的 `core/**` 路径仍需按"Go core 横切质量与代码规范"审核，并在 review 中说明该文件按通用标准处理。

## 验证要求

在 `core/` 目录下运行（与 CI `test-core-pr.yml` 行为一致）：

```bash
cd core && go test -race $(go list ./... | grep -v /proto/)
cd core && golangci-lint run
cd core && go build ./...
```

涉及 `.proto` 时，确认源定义与 Go 生成代码同步；具备完整工具链时在干净或临时工作区运行：

```bash
make proto
```

CI 覆盖 Ubuntu、Windows、macOS，跨平台构建标签和 `go test -race`。无法在本地覆盖的平台以 PR Actions 为准；Action 未运行或失败时必须如实记录，不能视为通过。

## 严重性

同 [`pr-review-lab` 严重性](../pr-review-lab/SKILL.md)：阻塞 / 必须修改 / 建议 / 细节 / 说明。每条 finding 必须使用显式前缀，并说明问题、影响和建议。

## 顶层正文结构

当 `pr-review-core` 为唯一适用 skill 时使用以下审核范围；跨语言 PR 与 `pr-review-lab` 的正文项合并。

```markdown
## 审核结论: 通过 / 要求修改 / 评论

## 本次 PR 解决的问题

<一句话概述>

## 审核范围

- [x] 平台与进程生命周期: <详情>
- [x] 控制台与诊断日志: <详情>
- [x] Protobuf 与生成代码: <详情>
- [x] Go 模块、工具链与构建: <详情>
- [x] Go core 横切质量与代码规范: <详情>
- [x] 安全、性能与跨平台兼容: <详情>
- [x] 测试与验证: <详情>

## 未挂行发现

<没有则写：无>

## 验证

- [ ] go test -race
- [ ] golangci-lint run
- [ ] go build ./...
- [ ] PR Actions 已通过，或已记录未运行/失败状态
```

## 完成验证

阶段一（审查）：

- [ ] 当前 PR 意图、base、head 和文件范围已确认
- [ ] 已在当前工作区 checkout 到 PR head 分支，关键发现在完整文件上验证
- [ ] 所有适用 Go core 领域和横切质量已审核
- [ ] 跨语言变更已读取 `pr-review-lab` 并覆盖 Python 范围
- [ ] 阻塞和必须修改 finding 已解决或明确记录
- [ ] 本地验证结果和 PR Actions 状态已记录
- [ ] 审查意见草稿已输出，并已明确告知用户尚未发布

阶段二（发布，需用户确认）：

- [ ] 用户已明确要求发布，且已确认 `event` 类型
- [ ] 可定位 finding 已提交为 inline comments
- [ ] 最终 review 已通过 Reviews API 或 `gh pr review` 提交
- [ ] 已切回原分支

## 参考

- [pr-review-lab](../pr-review-lab/SKILL.md) — 审核流程、严重性、提交审核格式（共享）；当 PR 包含 `swanlab/**` 变更时必读其 [`references/domain-contracts.md`](../pr-review-lab/references/domain-contracts.md)

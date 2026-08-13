# SwanLab Go core 领域审核契约

先根据 PR 文件列表判断适用范围，只读取与变更相关的章节。所有 Go 应用代码还需检查"Go core 横切质量与代码规范"；涉及凭据、网络、日志、文件路径或用户数据时，再检查"安全与隐私"（同 `pr-review-lab`）。

未命中具体章节的变更仍需结合调用关系、相邻模块和测试判断。本文没有定义的项目契约不得自行推断为硬性规则。

## 平台与进程生命周期

适用路径：`core/internal/pkg/process/**`。

- 构建标签分区必须为每个目标平台产生且仅产生一个实现（`//go:build linux`、`//go:build windows`、`//go:build !linux && !windows`）；不带构建标签的文件不得重复平台逻辑。
- `NotifyOnParentExit` 是公共入口，平台 `notifyOnParentExit` 为内部实现。`parentPID` 校验、注册后二次确认（race re-check）和 PID 变更检测必须保持完整。
- Linux 路径：`PR_SET_PDEATHSIG` 绑定调用时的父进程（非任意 PID），必须保留 prctl 后的二次 PID 确认；`SIGUSR1` 可能来自其他进程，收到信号后仍须校验父 PID。
- Windows 路径：以 `SYNCHRONIZE` 权限打开父进程句柄，`WaitForSingleObject` 等待 signaled；不能依赖 Windows 上不随父进程退出变化的 PPID。句柄必须由等待 goroutine 独占并在退出或失败后关闭。
- Poll 兜底路径：轮询间隔在响应速度和唤醒开销间取平衡；不能引入忙等或无限循环。
- channel 只能关闭一次；等待 goroutine 不得泄漏——`signal.Stop`、`CloseHandle`、`ticker.Stop` 在成功和失败路径上都要执行。
- 本包只检测父进程退出并关闭通知 channel，不直接终止进程或执行业务清理。退出顺序、清理超时和是否强制退出由调用方决定。
- 平台变更需要同步检查对应测试文件（`*_test.go`、`*_integration_test.go`）的构建标签是否匹配。

## 控制台与诊断日志

适用路径：`core/internal/pkg/console/**`（含 `log/`）。

- 与 Python `swanlab/sdk/internal/pkg/console` 的镜像关系必须保持：`Info`/`Debug`/`Warning`/`Error`/`Trace` 语义、`SWANLAB_DEBUG` 和 `NO_COLOR` 环境变量、loguru 风格格式（时间戳 | 级别 | 位置 - 消息）。
- `console.Init`/`Reset` 和 `log.Init`/`Reset` 必须保持对称。`Logger` 持有文件句柄，Go 无析构函数，重新绑定前必须先 `Reset` 关闭旧句柄。
- 日志文件权限强制 0600（`os.Chmod` 不受 umask 影响）；轮转（`maxBytes` × `backupCount`）保持备份链完整。
- 内存缓冲上限（`memoryCapacity`）超出时丢弃最旧记录，不能无界增长。
- `console.mu`（终端写入互斥）和 `log.Logger.mu`（文件写入互斥）必须序列化输出，且不在快速路径死锁。
- `callerLocation` 必须跳过 `console` 包自身的栈帧；保留 `shortFunc` 的截断逻辑。
- `Trace` 必须处理 nil error；终端显示短消息，诊断日志文件额外记录完整 goroutine 栈（`debug.Stack()`）。
- 修改日志级别方法（`Debug`/`Info`/`Warning`/`Error`/`Critical`）时保持级别名义化设计——当前均原样写入 msg，对应 Python 仅输出消息的 formatter。
- 测试用 `New()` 构造独立实例以隔离；包级 `std` 默认实例的测试用 `Reset`/`Init` 前后清理。

## Protobuf 与生成代码

适用路径：`core/proto/**`、`protos/**`、`scripts/generate_protos.py`。

- 详细的 Protobuf 兼容规则见 [`pr-review-lab` 的"Protobuf 与生成代码"](../../pr-review-lab/references/domain-contracts.md) 章节（已发布 `v1` 消息向后兼容、字段号不复用、append-only 旧记录可解析等）。
- `core/proto/**` 是 `make proto` 生成的 Go 代码，不得手动编辑。变更 `protos/**` 源定义后必须重新生成并同时提交 Python 和 Go 两侧输出。
- `scripts/generate_protos.py` 对 Go 生成代码的导入路径不做改写（仅改写 Python 侧）；不要意外破坏 Go 的 package 声明。
- 保持 `.golangci.yml` 中 `core/proto` 的 generated 排除（`exclusions.generated: lax`）；不要移除。
- Go proto 包路径（`github.com/swanhubx/swanlab/core/proto/swanlab/...`）与 `goimports` local-prefixes 配置一致。

## Go 模块、工具链与构建

适用路径：`core/go.mod`、`core/go.sum`、`core/.golangci.yml`、`core/.go-version`、`core/Makefile`、`Makefile`（`core-*` target）、`.github/workflows/test-core-pr.yml`。

- `core/.go-version` 与 `go.mod` 的 `go` 指令保持一致方向：CI 通过 `go-version-file: core/.go-version` 读取版本，修改时同步检查 `go.mod` 最低版本要求。
- 依赖变更同步 `go.mod` 与 `go.sum`；`modules-download-mode: readonly` 下不能引入未在 `go.sum` 中的依赖。
- golangci-lint 版本需匹配 `.golangci.yml` 的 `version: "2"`；修改 linter 列表或设置时确认 CI 的 golangci-lint-action 版本兼容。
- `.golangci.yml` 的设置（`gocyclo` 30、`lll` 160、`mnd` scope、`goimports` local-prefixes）是审核基准，不另造格式规则。
- Makefile `core-*` target（`core-test`、`core-lint`、`core-build`、`core-tidy`）与 CI 行为一致时优先使用 Makefile；注意 `core-test` 不含 `-race`，需要竞争检测时用显式 `go test -race`。
- CI matrix 覆盖 Ubuntu、Windows、macOS；`core/internal/pkg/process` 的平台相关源文件依赖跨系统编译和测试，不能只在一个平台验证。

## Go core 横切质量与代码规范

适用路径：所有 `core/**` Go 代码、配置和依赖变更。

- 从正确性、兼容性、可读性、架构、安全和性能判断具体回归，不以个人偏好提出 finding。
- 以 `core/.golangci.yml` 启用的 linters（`errcheck`、`govet` 含 `shadow`、`staticcheck`、`unused`、`gocyclo`、`lll`、`mnd`、`revive` 等）为准，不另造格式规则。
- 错误必须显式处理；`_ =` 抑制和宽泛 `panic`/`recover` 需要具体理由。现有代码中 `_ =` 用于尽力而为的清理（如 `os.Remove`/`os.Rename` 失败可忽略），不应扩大到掩盖需要上抛的错误。
- goroutine 必须有可识别的终止路径；channel 发送/关闭必须避免跨发送方关闭或向已关闭 channel 发送。
- 资源（文件句柄、OS 句柄、ticker、signal listener）在成功和异常路径上都要释放。
- 导出的包 API 需要文档注释；保持现有的中文解释性注释（记录清理顺序、平台差异、Python 镜像理由、竞态窗口说明等）。
- `internal/` 下的包不得在模块外部导入；`core` 是 Python SDK 启动的子进程，不是可导入的库。
- 不新增调试输出、临时分支、无用 import 或由本次变更产生的死代码。
- 复杂度是审查信号而非机械行数门禁；优先最小正确改动，不为单一用例增加抽象层。

## 测试与验证

适用路径：`core/**_test.go`、CI workflow。

- 行为变更应在对应包下增加回归测试，测试文件与源码同包。
- 遵循标准 `testing` 包约定：`t.TempDir()` + `t.Cleanup()` LIFO 管理文件句柄（Windows 不允许删除打开的文件）；使用表驱动用例。
- 断言消息包含实际值（如 `got %q`、`got %d, want %d`），不使用无上下文的 `t.Fatal`。
- 平台相关测试文件必须带匹配的构建标签（如 `//go:build linux`）。
- 竞争、并发和生命周期变更同时覆盖成功与失败路径，本地用 `go test -race` 验证。
- 测试不依赖真实外部资源（真实父进程 PID 除外，使用 `os.Getppid()` 等可控值）。
- 本地验证结果不能替代 CI 的跨操作系统矩阵。

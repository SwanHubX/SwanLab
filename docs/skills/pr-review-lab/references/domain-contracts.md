# SwanLab 领域审核契约

先根据 PR 文件列表判断适用范围，只读取与变更相关的章节。所有 Python 应用代码还需检查“横切质量与代码规范”；涉及凭据、网络、日志、文件路径或用户数据时，再检查“安全与隐私”。

未命中具体章节的变更仍需结合调用关系、相邻模块和测试判断。本文没有定义的项目契约不得自行推断为硬性规则。

## 公共 API 与类型兼容

适用路径：`swanlab/__init__.py`、`swanlab/__init__.pyi`、`swanlab/sdk/__init__.py`、`swanlab/sdk/typings/**` 及公开类和函数。

- 顶层运行时导出、`__all__` 和 `swanlab/__init__.pyi` 保持同步。
- 新增或修改公开参数时检查默认值、返回值、异常、文档和类型签名是否一致。
- 已发布 API 的删除、重命名或语义变化需要兼容或明确的弃用路径，不能静默破坏用户代码。
- 保持 Python 3.9-3.14 可用，不引入只在单一 Python 版本成立的行为。
- 修改 `swanlab.run`、`get_run()` 等动态属性时覆盖初始化前、运行中和结束后的行为。

## Run 生命周期与运行时组件

适用路径：`swanlab/sdk/cmd/**`、`swanlab/sdk/internal/run/**`、`swanlab/sdk/internal/context/**`、`swanlab/sdk/internal/bus/**`。

- `init`、`log`、`finish` 和 reinit 的状态转换清晰，失败路径不会遗留半初始化的全局状态。
- `online`、`offline`、`local`、`disabled` 四种模式保持各自边界；disabled 模式不应产生不必要的网络或文件副作用。
- 运行时组件按依赖顺序启动、按反向顺序停止；停止时先完成异步任务和末尾事件，再关闭消费者。
- terminal proxy、后台线程、executor、watcher、signal 和 atexit hook 都有确定的安装与清理路径。
- 主线程、非主线程、fork/spawn、异常退出和重复 finish 场景不能死锁、重复上传或丢失数据。
- 修改全局或上下文状态时检查多线程和 asyncio 隔离。

## Record、Core 与传输

适用路径：`swanlab/sdk/internal/core_python/**`、`transport/**`、`store/**`、Record builder 和 consumer。

- Record 保持可顺序写入、读取、重放和断点续传；序号、批次和完成状态不能倒退或重复计数。
- scalar、media、column、config、log、save 和生命周期记录的转换链路保持完整。
- 上传重试责任不能在多层重复：上传 API 使用既定请求策略，sender 决定可重试错误，transport 负责批次重试。
- 当前契约中 5xx/基础设施错误应允许 transport 重试，4xx 业务拒绝不应无界重试；不能用宽泛异常吞掉需要上抛的错误。
- 成功、部分失败和最终失败时，tracker、buffer、落盘状态与实际上传结果一致。
- 批处理、队列和并发上传应有边界，避免大数据量下无界内存、忙等、锁竞争或线程泄漏。
- 文件与媒体上传检查 Windows/POSIX 路径、跨机器 sync、重复上传和幂等行为。

## 配置、认证与 HTTP 客户端

适用路径：`swanlab/sdk/internal/settings/**`、`swanlab/sdk/internal/pkg/client/**`、`swanlab/sdk/internal/pkg/nrc/**`。

- 保持 Settings 声明的配置源优先级；有意调整时同步兼容测试和用户文档。
- Settings 负责类型、格式和默认值；目录创建、网络请求等业务副作用留在初始化阶段。
- `merge_settings` 只合并用户显式设置字段，嵌套配置不能被意外整块覆盖。
- API host、web host、环境变量、YAML、Secret 和 netrc 的回退关系明确，切换 host 时不能复用错误凭据。
- HTTP 请求设置合理的 timeout、retry 和错误映射；请求级上下文在成功和异常后都要恢复。
- 日志和异常不能暴露 API key、认证 token、cookie、presigned URL 或完整敏感响应。

## Public API 与 CLI

适用路径：`swanlab/api/**`、`swanlab/cli/**`。

- API 输入校验、分页、迭代器、空值、过滤条件和错误类型保持稳定。
- API 与 CLI 复用统一的客户端、host 和认证逻辑，不在命令中复制请求实现。
- CLI 的参数名、别名、默认值、帮助文本、输出和退出码保持兼容。
- 非交互环境不能意外等待输入；密码和 API key 输入不得明文回显。
- self-hosted 与默认云端地址的请求和展示链接不能混用。

## 媒体、转换与文件

适用路径：`swanlab/sdk/internal/run/transforms/**`、媒体类型、`save`、文件 watcher 和上传逻辑。

- Image、Audio、Video、Text、Html、ECharts、Object3D、Molecule 的输入类型、输出文件和 Record 元数据一致。
- 可选依赖应按需加载；未安装媒体或框架依赖时不能影响基础 `import swanlab`。
- 文件名、扩展名、MIME、大小、hash、caption 和 step 在转换与上传链路中不丢失。
- 文件句柄、临时文件、watcher、线程和大对象有释放路径。
- 大文件和批量媒体避免一次性读入内存，并遵守上传大小、分片和批次限制。
- `save` 的 glob、base path、符号链接和跨平台路径不能越过预期边界或上传错误文件。

## 集成、转换器与插件

适用路径：`swanlab/integration/**`、`swanlab/converter/**`、`swanlab/plugin/**`。

- 第三方框架未安装时保持基础包可导入，并提供清晰错误而非导入期崩溃。
- 框架 callback 只在正确进程或 rank 记录，避免重复初始化、重复 step 和重复 finish。
- monkey patch 保留原库行为，重复调用尽量幂等，并明确安装时机和生命周期。
- 第三方配置和指标转换保留关键字段、step、命名空间和数据类型。
- callback 的异步任务在 Python shutdown 和 Run finish 时不会静默丢失。
- 新增集成或插件时补充对应隔离测试，不依赖真实外部服务。

## Protobuf 与生成代码

适用路径：`protos/**`、`swanlab/proto/**`、`core/proto/**`、`scripts/generate_protos.py`。

- `protos/**` 是协议源；变更后同步提交 Python 和 Go 生成代码。
- 已发布 `v1` 消息保持向后兼容，不复用或改变已有字段号和枚举值语义。
- Record 顶层信封和 append-only 日志保持旧记录可解析、可跳过未知字段。
- 媒体和文件仍通过路径与元数据引用，不把大二进制意外内联进 Record。
- 修改生成脚本时检查导入路径修复、`.pyi` 生成和 package `__init__.py`。

## 构建、依赖与发布

适用路径：`pyproject.toml`、`uv.lock`、`Makefile`、`swanlab/package.json`、`.github/workflows/**`。

- 依赖变更同步 `pyproject.toml` 与 `uv.lock`，区分核心依赖、dev dependency 和 optional extra。
- 基础安装不能强制引入只服务于单个媒体类型或第三方框架的重依赖。
- wheel 包含公共 stub、Proto 和运行时所需资源；版本仍由既定单一来源读取。
- 发布和构建修改检查 tag 解析、版本写入、wheel 产物和 PyPI/GitHub Release 流程。
- 保持 Linux、Windows、macOS 与 Python 3.9-3.14 的语法、路径、编码和依赖兼容。

## 测试与验证

适用路径：`tests/**`、pytest 配置和 CI workflow。

- 行为变更应在对应 `tests/unit` 路径增加回归测试，测试目录与源码模块尽量对应。
- 默认使用 mock、`tmp_path` 和 disabled/local 模式隔离真实网络、用户凭据和全局文件。
- 生命周期、并发、重试、跨平台路径和兼容性变更同时覆盖成功与失败路径。
- 测试清理全局 Settings、Run、环境变量、线程、patch 和文件，避免依赖执行顺序。
- 本地验证结果不能替代 CI 的操作系统和 Python 版本矩阵。

## 横切质量与代码规范

适用路径：所有 Python 代码、配置和依赖变更。

- 从正确性、兼容性、可读性、架构、安全和性能判断具体回归，不以个人偏好提出 finding。
- 以 `pyproject.toml` 中 Ruff、basedpyright 和 pytest 配置为准，不另造格式规则。
- 新代码保持类型清晰；`Any`、`cast`、`type: ignore` 和宽泛异常必须有具体理由和最小作用域。
- 避免导入期网络、文件写入、线程启动和重型可选依赖加载。
- 资源获取与释放成对；异常路径不能跳过锁释放、context reset、线程停止或文件关闭。
- 不新增调试输出、临时分支、无用 import 或由本次变更产生的死代码。
- 复杂度是审查信号而非机械行数门禁；优先最小正确改动，不为单一用例增加抽象层。

## 安全与隐私

适用条件：凭据、认证、网络、日志、HTML、文件路径、上传、环境探测或用户数据发生变化。

- API key、token、cookie、密码和 Secret 不写入日志、异常、测试快照或上传元数据。
- netrc、配置文件和 Secret 文件遵守最小权限，不能写入错误用户目录或意外覆盖其他 host 凭据。
- 自定义 host、webhook 和外部 URL 检查协议、凭据发送范围、重定向和 SSRF 风险。
- HTML、终端日志、配置和第三方框架数据视为不可信输入，上传或展示前保持明确的安全边界。
- `save`、sync 和解压/复制逻辑防止路径穿越、符号链接逃逸和任意文件覆盖。
- 环境与硬件探测只采集产品需要的数据，不上传本地敏感路径、环境变量值或用户私有内容。

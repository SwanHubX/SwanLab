# SwanLab Protobuf Definitions

本目录包含 SwanLab 运行时产出的所有结构化数据的 Protobuf 类型定义。

## 设计思路

所有消息以 **append-only** 的方式写入运行日志文件，每条写入均封装为一个顶层 `Record`。
`Record` 内部通过 `oneof record_type` 区分消息类型，`num` 字段保证全局单调递增，支持断点续传与日志重放。

媒体数据（图像、音频、视频）**不内联进 Record**，而是先写入磁盘，Record 只保存文件路径引用和元数据。
这样做的好处是：Record 文件保持轻量可流式读取，文件上传器可与 metric 记录完全解耦。

## 目录结构

```
protos/
└── swanlab/
    ├── record/
    │   └── v1/
    │       └── record.proto            # 顶层信封 Record
    │
    ├── run/
    │   └── v1/
    │       └── run.proto               # Run 生命周期（StartRecord、FinishRecord、RunState、ResumeMode）
    │
    ├── metric/
    │   ├── column/
    │   │   └── v1/
    │   │       └── column.proto        # 列定义（ColumnRecord），描述 metric 的展示与归类配置
    │   └── data/
    │       └── v1/
    │           └── data.proto          # 指标数据（ScalarRecord、MediaRecord）
    │
    ├── config/
    │   └── v1/
    │       └── config.proto            # 实验配置引用（ConfigRecord、UpdateType）
    │
    ├── env/
    │   └── v1/
    │       └── env.proto               # 软件环境与主机元数据引用（MetadataRecord、RequirementsRecord、CondaRecord）
    │
    ├── terminal/
    │   └── v1/
    │       └── log.proto               # 终端代理捕获的 stdout/stderr 输出（LogRecord、LogLevel）
    │
    ├── save/
    │   └── v1/
    │       └── save.proto              # 文件保存记录（SaveRecord、SavePolicy）
    │
    ├── settings/
    │   ├── core/
    │   │   └── v1/
    │   │       └── core.proto          # Core 同步服务配置（CoreSettings）
    │   └── probe/
    │       └── v1/
    │           └── probe.proto         # Probe 硬件采集服务配置（ProbeSettings）
    │
    └── grpc/
        ├── core/
        │   ├── v1/
        │   │   ├── core.proto          # CoreService gRPC 服务（接收实验记录）
        │   │   └── sync.proto          # CoreSyncService gRPC 服务（启动日志同步）
        │   └── ...
        └── probe/
            └── v1/
                └── probe.proto         # ProbeService gRPC 服务（启动硬件采集）
```

## 消息流

```
swanlab.init(...)
  → Record { start:        StartRecord        }   # Run 基本信息，写入一次
  → Record { metadata:     MetadataRecord     }   # 主机与环境快照引用，写入一次
  → Record { requirements: RequirementsRecord }   # Python 依赖引用，写入一次
  → Record { conda:        CondaRecord        }   # Conda 环境引用，写入一次
  → Record { config:       ConfigRecord       }   # 实验 config（INIT），写入一次

# 训练循环中，首次出现某 key 时先发送列定义
  → Record { column: ColumnRecord { key="loss", type=SCALAR, ... } }
  → Record { column: ColumnRecord { key="img",  type=IMAGE, ... } }
# 每个指标独立生成一条 Record，swanlab.log({"loss": 0.5, "img": Image(...)}, step=10) 产生两条：
  → Record { scalar: ScalarRecord { key="loss", step=10, ... } }
  → Record { media:  MediaRecord  { key="img",  step=10, ... } }

# 终端代理（每行）
  → Record { log: LogRecord { line="...", level=INFO } }

# 运行中 config 更新
  → Record { config: ConfigRecord { update_type=PATCH, ... } }

# 文件保存
  → Record { save: SaveRecord { name="...", policy=NOW } }

swanlab.finish()
  → Record { finish: FinishRecord { state=FINISHED } }
```

## Record 信封

`Record` 是所有消息的顶层信封，通过 `oneof record_type` 互斥区分消息类型：

| 字段号 | 字段名 | 消息类型 | 说明 |
|--------|--------|----------|------|
| 10 | `start` | StartRecord | Run 创建 |
| 11 | `finish` | FinishRecord | Run 结束 |
| 12 | `column` | ColumnRecord | 指标列定义 |
| 13 | `scalar` | ScalarRecord | 标量数据 |
| 14 | `media` | MediaRecord | 媒体数据 |
| 15 | `config` | ConfigRecord | 实验配置 |
| 16 | `log` | LogRecord | 终端输出 |
| 17 | `metadata` | MetadataRecord | 主机元数据引用 |
| 18 | `requirements` | RequirementsRecord | Python 依赖引用 |
| 19 | `conda` | CondaRecord | Conda 环境引用 |
| 20 | `save` | SaveRecord | 文件保存 |

`num` 字段保证全局单调递增，用于去重和断点续传。特殊值：Start=-1, Finish=-2, Config=-3, Metadata=-4, Requirements=-5, Conda=-6。

## gRPC 服务

### CoreService

`grpc/core/v1/core.proto` 定义核心业务接口，用于同步/异步接收实验记录：

| RPC 方法 | 请求类型 | 响应类型 | 说明 |
|----------|----------|----------|------|
| `DeliverRunStart` | DeliverRunStartRequest | DeliverRunStartResponse | 实验开始 |
| `UpsertColumns` | UpsertColumnsRequest | Empty | 定义指标列 |
| `UpsertScalars` | UpsertScalarsRequest | Empty | 记录标量值 |
| `UpsertMedia` | UpsertMediaRequest | Empty | 记录媒体值 |
| `UpsertConfigs` | UpsertConfigsRequest | Empty | 更新配置 |
| `UpsertLogs` | UpsertLogsRequest | Empty | 终端输出 |
| `UpsertMetadata` | UpsertMetadataRequest | Empty | 主机元数据更新 |
| `UpsertRequirements` | UpsertRequirementsRequest | Empty | 依赖更新 |
| `UpsertConda` | UpsertCondaRequest | Empty | Conda 环境更新 |
| `UpsertSaves` | UpsertSavesRequest | Empty | 文件保存 |
| `DeliverRunFinish` | DeliverRunFinishRequest | DeliverRunFinishResponse | 实验结束 |

### CoreSyncService

`grpc/core/v1/sync.proto` 定义日志同步服务接口：

| RPC 方法 | 请求类型 | 响应类型 | 说明 |
|----------|----------|----------|------|
| `DeliverSyncStart` | DeliverSyncStartRequest | Empty | 启动本地日志读取与云端同步 |
| `DeliverSyncFinish` | Empty | Empty | 结束日志同步 |

### ProbeService

`grpc/probe/v1/probe.proto` 定义硬件采集服务接口：

| RPC 方法 | 请求类型 | 响应类型 | 说明 |
|----------|----------|----------|------|
| `DeliverProbeStart` | DeliverProbeStartRequest | Empty | 启动硬件信息采集 |
| `DeliverProbeFinish` | Empty | Empty | 结束硬件信息采集 |

## 各包详解

### run — Run 生命周期

- **StartRecord**: 对应 `swanlab.init()`，包含 project、name、description、tags、group、resume 模式等
- **FinishRecord**: 对应 `swanlab.finish()`，包含运行状态和错误信息
- **RunState**: `RUNNING` / `FINISHED` / `CRASHED` / `ABORTED`
- **ResumeMode**: `NEVER` / `ALLOW` / `MUST`

### metric/column — 列定义

- **ColumnRecord**: 描述 metric key 的显示配置和归类规则，首次出现某 key 时发送
- **ColumnType**: `SCALAR` / `IMAGE` / `AUDIO` / `TEXT` / `VIDEO` / `ECHARTS` / `OBJECT3D` / `MOLECULE`
- **ColumnClass**: `CUSTOM`（用户自定义）/ `SYSTEM`（系统列）
- **SectionType**: `PINNED` / `HIDDEN` / `PUBLIC` / `SYSTEM`，控制指标在 UI 中的分组

### metric/data — 指标数据

- **ScalarRecord**: 标量记录，key + step + value（原生支持 IEEE 754 NaN/Infinity）
- **MediaRecord**: 媒体记录，key + step + value（包含 filename、sha256、size、caption）

### config — 实验配置

- **ConfigRecord**: 配置引用，实际数据写入 `files/config.yaml`，Record 只保存更新类型
- **UpdateType**: `INIT`（init 时完整 config）/ `PATCH`（运行中增量更新）

### env — 环境与元数据

- **MetadataRecord**: 主机与环境元数据引用 → `files/swanlab-metadata.json`
- **RequirementsRecord**: Python 包依赖引用 → `files/requirements.txt`
- **CondaRecord**: Conda 环境引用 → `files/conda.yml`

三者均只含 `timestamp` 字段，实际数据存储在对应文件路径中，Record 仅作引用/触发。

### terminal — 终端输出

- **LogRecord**: 终端单行输出，由 `console.proxy_type` 决定捕获 stdout/stderr/两者
- **LogLevel**: `INFO`（普通输出）/ `ERROR`（实验错误输出）

### save — 文件保存

- **SaveRecord**: 由 `swanlab.save()` 产生，包含 name、source_path、target_path
- **SavePolicy**: `NOW`（立即上传）/ `END`（运行结束时上传）/ `LIVE`（持续监听文件变化）

### settings — 服务配置

- **CoreSettings** (`settings/core/v1/core.proto`): Core 同步服务的启动配置
- **ProbeSettings** (`settings/probe/v1/probe.proto`): Probe 硬件采集服务的启动配置

## ColumnRecord 生成策略

默认 column 记录会在 ScalarRecord 和 MediaRecord 上传时由 core 自动生成。core 也允许手动上传 column record，携带更多自定义参数。

## 版本策略

所有 package 路径均含版本号（`v1`），后续不兼容变更通过新增 `v2` 目录演进，不修改已发布版本。

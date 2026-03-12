# SwanLab Protobuf Definitions

本目录包含 SwanLab 运行时产出的所有结构化数据的 Protobuf 类型定义。

## 设计思路

所有消息以 **append-only** 的方式写入运行日志文件，每条写入均封装为一个顶层 `Record`。
`Record` 内部通过 `oneof` 区分消息类型，`num` 字段保证全局单调递增，支持断点续传与日志重放。

媒体数据（图像、音频、视频）**不内联进 Record**，而是先写入磁盘，Record 只保存文件路径引用和元数据。
这样做的好处是：Record 文件保持轻量可流式读取，文件上传器可与 metric 记录完全解耦。

## 目录结构

```
protos/
└── swanlab/
    ├── record/
    │   └── v1/
    │       └── record.proto        # 顶层信封，所有消息的统一入口
    │
    ├── run/
    │   └── v1/
    │       └── run.proto           # Run 生命周期（RunRecord、FinishRecord）
    │
    ├── data/
    │   └── v1/
    │       ├── log.proto           # 用户 swanlab.log() 产生的 metric 记录
    │       ├── column.proto        # 列定义（ColumnRecord），描述 metric 的展示与归类配置
    │       ├── scalar.proto        # 标量值（数字 / 字符串 / bool）
    │       ├── image.proto         # 图像引用
    │       ├── audio.proto         # 音频引用
    │       ├── video.proto         # 视频引用
    │       ├── text.proto          # 文本值
    │       ├── table.proto         # 表格值
    │       └── echarts.proto       # ECharts 可视化配置
    │
    ├── config/
    │   └── v1/
    │       └── config.proto        # 用户传入的实验 config
    │
    └── system/
        └── v1/
            ├── hardware.proto      # 运行时硬件监控（StatsRecord）
            ├── env.proto           # 软件环境与主机元数据引用（MetadataRecord、RequirementsRecord、CondaRecord）
            └── console.proto       # 终端代理捕获的 stdout/stderr 输出
```

## 消息流

```
swanlab.init(...)
  → Record { run:          RunRecord          }   # run 基本信息，写入一次
  → Record { metadata:     MetadataRecord     }   # 主机与环境快照引用，写入一次
  → Record { requirements: RequirementsRecord }   # Python 依赖引用，写入一次
  → Record { conda:        CondaRecord        }   # Conda 环境引用，写入一次
  → Record { config:       ConfigRecord       }   # 实验 config，写入一次

# 训练循环中，首次出现某 key 时先发送列定义
  → Record { column: ColumnRecord { key="loss", type=FLOAT, ... } }
  → Record { column: ColumnRecord { key="img",  type=IMAGE, ... } }
# 每个指标独立生成一条 LogRecord，swanlab.log({"loss": 0.5, "img": Image(...)}, step=10) 产生两条：
  → Record { log: LogRecord { key="loss", step=10, value=ScalarValue{...} } }
  → Record { log: LogRecord { key="img",  step=10, value=ImageValue{...}  } }

# 终端代理（每行）
  → Record { console: ConsoleRecord { line="...", stream=STDOUT } }

# 硬件监控（每 interval 秒）
  → Record { stats: StatsRecord { cpu=..., gpus=[...], ... } }

# 运行中 config 更新
  → Record { config: ConfigRecord { update_type=PATCH, ... } }

swanlab.finish()
  → Record { finish: FinishRecord { state=FINISHED } }
```

## 版本策略

所有 package 路径均含版本号（`v1`），后续不兼容变更通过新增 `v2` 目录演进，不修改已发布版本。

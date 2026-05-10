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
    │       └── record.proto            # 顶层信封，所有消息的统一入口
    │
    ├── run/
    │   └── v1/
    │       └── run.proto               # Run 生命周期（StartRecord、FinishRecord）
    │
    ├── metric/
    │   ├── column/
    │   │   └── v1/
    │   │       └── column.proto        # 列定义（ColumnRecord），描述 metric 的展示与归类配置
    │   └── data/
    │       └── v1/
    │           └── data.proto          # 用户 swanlab.log() 产生的指标数据（ScalarRecord、MediaRecord）
    │
    ├── config/
    │   └── v1/
    │       └── config.proto            # 用户传入的实验 config（ConfigRecord）
    │
    ├── env/
    │   └── v1/
    │       └── env.proto               # 软件环境与主机元数据引用（MetadataRecord、RequirementsRecord、CondaRecord）
    │
    ├── terminal/
    │   └── v1/
    │       └── log.proto               # 终端代理捕获的 stdout/stderr 输出（LogRecord、StreamType）
    │
    └── save/
        └── v1/
            └── save.proto              # 文件保存记录（SaveRecord）
```

## 消息流

```
swanlab.init(...)
  → Record { start:        StartRecord        }   # run 基本信息，写入一次
  → Record { metadata:     MetadataRecord     }   # 主机与环境快照引用，写入一次
  → Record { requirements: RequirementsRecord }   # Python 依赖引用，写入一次
  → Record { conda:        CondaRecord        }   # Conda 环境引用，写入一次
  → Record { config:       ConfigRecord       }   # 实验 config，写入一次

# 训练循环中，首次出现某 key 时先发送列定义
  → Record { column: ColumnRecord { key="loss", type=FLOAT, ... } }
  → Record { column: ColumnRecord { key="img",  type=IMAGE, ... } }
# 每个指标独立生成一条 Record，swanlab.log({"loss": 0.5, "img": Image(...)}, step=10) 产生两条：
  → Record { scalar: ScalarRecord { key="loss", step=10, ... } }
  → Record { media:  MediaRecord  { key="img",  step=10, ... } }

# 终端代理（每行）
  → Record { log: LogRecord { line="...", stream=STDOUT } }

# 运行中 config 更新
  → Record { config: ConfigRecord { update_type=PATCH, ... } }

# 文件保存
  → Record { save: SaveRecord { name="...", policy=NOW } }

swanlab.finish()
  → Record { finish: FinishRecord { state=FINISHED } }
```

## 版本策略

所有 package 路径均含版本号（`v1`），后续不兼容变更通过新增 `v2` 目录演进，不修改已发布版本。


## ColumnRecord生成策略

默认 column 记录会在 ScalarRecord 和 MediaRecord 上传时由 core 自动生成。core 也允许手动上传 column record，携带更多自定义参数。
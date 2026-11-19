# tests/benchmark

SwanLab SDK 性能基准。目录按被测对象分层：端到端场景、多进程共享实验、Go 侧微基准、
Python 侧微基准，各自独立运行。

所有联网基准都跑在**本地 Mock Cloud** 上：它只监听 `127.0.0.1`，只累加计数、不保存明细，
因此不会访问任何真实后端。子进程的 `SWANLAB_API_HOST`、`SWANLAB_API_KEY`、`SWANLAB_ROOT`
和 `cwd` 由 `common.build_child_env()` 统一设置，并在子进程里回读 `Settings` 做确认。

## 目录

```text
tests/benchmark/
├── common.py                    # 公共工具：统计、进程查找、内存/CPU 采样、子进程环境、结果落盘
├── mock_cloud/server.py         # 本地 Mock Cloud（只计数，带 /__mock/stats）
├── e2e/
│   ├── bench_e2e.py             # 场景 A/B/C × No-SDK / Python core / Go core
│   ├── bench_async_log.py       # log vs async_log 对比
│   └── bench_skip_store_e2e.py  # skip_store 与媒体/文件混合写入对比（responses mock）
├── multi_process/
│   └── bench_multi_process.py   # Go 共享 / Go 独占 / Python core 三种多进程方式
├── python_core/                 # 只测 Python 侧实现的微基准（步数去重、存储、格式化）
└── results/                     # 结果 JSON（不入库）
```

Go 侧微基准放在对应的 Go 包里（如 `core/internal/metrics/bench_test.go`），因为它们需要
访问包内私有结构。

```bash
# 存储（写 / 回读 / 异步刷盘），固定 -benchtime=1x -count=3
go test -C core -run '^$' -bench 'Store' -benchtime=1x -count=3 ./internal/store

# 步数去重
go test -C core -run '^$' -bench 'Steps|Accept' -benchtime=1x -count=3 ./internal/metrics

# Python 侧对照（同样 100 万条 × 120 字节，含回读吞吐）
uv run pytest tests/benchmark/python_core/store/bench_store_fsync.py -q -s
```

## 运行

```bash
# 小规模验证（秒级），先确认链路通
uv run python tests/benchmark/e2e/bench_e2e.py --scenario C --steps 100 --keys 100 --core both

# 场景 A/B/C 正式跑（每个 case 重复 3 次取中位数）
uv run python tests/benchmark/e2e/bench_e2e.py --scenario A --core both --repeats 3
uv run python tests/benchmark/e2e/bench_e2e.py --scenario B --core both --repeats 3
uv run python tests/benchmark/e2e/bench_e2e.py --scenario C --core both --repeats 3
uv run python tests/benchmark/e2e/bench_e2e.py --scenario C --skip-store

# 场景 C 混合档（图片 + run.save）
uv run python tests/benchmark/e2e/bench_e2e.py --scenario C --core both --repeats 3 \
  --media-every 10 --save-count 2

# 给 mock 加 20ms 响应延迟，观察上传队列是否堆积
uv run python tests/benchmark/e2e/bench_e2e.py --scenario C --latency-ms 20

# 多进程（DDP 场景），--repeats 3 取墙钟中位数并输出离散度
uv run python tests/benchmark/multi_process/bench_multi_process.py --ranks 8 --steps 100 --keys 1000 --repeats 3
uv run python tests/benchmark/multi_process/bench_multi_process.py --ranks 16 --steps 100 --keys 1000 --sleep 0.2 --repeats 3
```

`make bench` 等价于 `uv run pytest tests/benchmark`，只跑带 `test_*` 入口的历史用例；
`bench_e2e.py` 和 `bench_multi_process.py` 是命令行脚本，不参与 pytest 收集。

## 场景

| 场景 | 负载 | 说明 |
|---|---|---|
| A | 100,000 步 × 9 指标 = 90 万条 | 高频少指标（强化学习、小步迭代） |
| B | 7 步 × 100,000 指标 = 70 万条 | 低频大批量；只报 init/finish/内存/体积，不报延迟分位 |
| C | 1,000 步 × 100 指标 = 10 万条 | 常规规模，支持 `--skip-store`、`--media-every`、`--save-count` |

混合档（`--media-every N --save-count K`）每 N 步记一张 16×16 PNG，并用 `run.save(policy="end")`
登记 K 个 64 KiB 文件（finish 时统一上传）。脚本会额外断言 media 条数与完成文件数。

`--profile` 固定 `core.record_batch` 与 `core.record_interval`，两档结果不要放在同一张表里比：

| 档位 | batch | interval | 用途 |
|---|---|---|---|
| `default`（默认） | 10,000 | 1.5 | 接近真实使用 |
| `isolated` | 10,000 | 3,600 | 仅用于对照：Python 侧会攒到 finish 统一上传，而 Go 侧满 batch 仍立即发批；弱网反压请用 `--latency-ms` |

`record_batch` 会影响多个环节：Python `BackgroundConsumer` 的小批固定是 `min(100, batch)`、
gRPC 批是 `min(2000, batch)`，Go 侧队列容量是 `8 × batch`。两个 core 的发送时机也不一样
（Python `Transport.put()` 不唤醒条件变量，靠定时或 finish 刷批；Go 攒满 `batch` 或等满
`interval` 即发）。

## 校验项（任一不通过则本轮数据不可用）

`bench_e2e.py` 每个 case 都会检查：

| 检查 | 含义 |
|---|---|
| 实际使用的 core 与预期一致 | Go 二进制缺失时 SDK 会静默退回 Python core，必须挡住 |
| 子进程未连真实后端 | 回读 `Settings`，`api_host` / `api_key` / `on_online_error` 必须是注入的值 |
| Go core 进程数 / Python 无 core 进程 | 确认进程模型与预期一致 |
| mock 收到标量条数 | 应等于 `steps × keys`（按 record 类型分别计数） |
| 按类型核对 | `scalar` / `log` / `media` / `histogram` 四类，mock 收到的条数必须等于本地回读条数 |
| 上传进度核对 | `uploaded_records` 不超过 `mock house 总数 + column + save`，详见下节；`--skip-store` 时跳过 |
| 列定义覆盖 | `/series` 收到的 key 必须覆盖 `run.log` 用过的全部 key，缺失即失败，多出的 key 记入结果 |
| 本地文件标量条数 | `DataStoreReader` 回读后应等于 `steps × keys`（`--skip-store` 时跳过） |
| 上传全程未被暂停 | 仅 `--latency-ms 0` 时断言 `upload.paused` 不能出现 true；弱网档保留该字段作观测 |
| core 日志无报错 / 子进程无异常堆栈 | 排除失败但被吞掉的场景 |

多进程脚本额外检查「共享模式下只有 1 个 core 进程」「每个 rank 各启动 1 个」「各 rank 用的
core 进程互不相同」，依据来自子进程自己上报的 pid（配合父进程全生命周期累积的 `core_ever_seen`）。
各 rank 先按 `--stagger`（默认 `0.15s`）错峰完成 `init()`，经启动屏障同步放行后才开始计时与打点，
以剥离启动错峰对稳态并发 `wall_s` 的稀释；本地文件回读统一在采样停止后由父进程离线完成。
`--repeats N` 时按墙钟中位数选代表运行，并输出所有轮次与相对离散度（`(max-min)/中位数`）；
控制台还会打印每个拓扑的逐 rank 队列峰值（`max_buffered_bytes` / `max_pending_records`）以及 Mock Cloud CPU，
用来区分「队列打满」、「CPU/磁盘争用」与「Mock 自身开销」。

上传进度核对的上界来自两个 core 的实际分发路径（`core/internal/uploader/dispatch.go` 与
`swanlab/sdk/internal/core_python/transport/dispatch/`）：`uploaded_records` 统计已处理完的
记录，包含

- 上报到 `/house/metrics` 的 `scalar` / `log` / `media` / `histogram`；
- 上报到 `/series` 的 `column`（mock 不按条计数，用本地回读条数补上）；
- `save` 记录：只要文件传完（或被本地规则拒绝）就计数，与 mock 的 `files_completed`
  **不是一回事**（配置 / requirements / metadata 这类 save 走 profile 接口，不产生
  `/files/complete` 调用，mock 计到 0 条）。

`start` / `finish` 只落盘，不进上传队列，因此不计入。`uploaded_records` 是每 0.2 秒
采样一次的快照，收尾阶段可能晚于实际完成，所以只校验上界，不要求相等。

## 测量方式

- `run.log()` 延迟用 `perf_counter` 逐次采样，报 mean/P50/P95/P99（微秒）。它反映的是主线程
  是否被后台消费拖慢，而不是 core 的绝对吞吐。
- macOS 上优先通过 `proc_pid_rusage`（50ms 间隔 + 阶段边界精确打点）同步采集各进程
  `RSS`（`ri_resident_size`）、扣除跨进程共享只读页后的物理内存 `phys_footprint`（`ri_phys_footprint`）
  以及经 `mach_timebase_info` 换算的纳秒级累计 CPU 时间；其余平台回退 `ps -o pid=,rss=,time=`。
- 总内存报**同一采样时刻的并发总和峰值**（`peak_rss` / `peak_phys`），避免短命进程退出或
  错峰导致的漏计；Mock Cloud 进程独立采样统计，不计入被测总值。
- 结果写入 `results/*.json`，含 git sha、core 二进制 sha256、`GOGC` 和实际使用的参数。

## 注意事项

1. 每组开跑前会先清掉残留的 core 进程（`detach=True` 默认要空闲 300s 才退出）。残留实例会
   带着上一轮的指标状态和内存继续服务，会污染多进程结果。
2. Mock Cloud 的 `/api/health`、`/api/house/health` 和未注册路径都必须返回 200 或 4xx；
   返回 5xx 或断连会让 Go core 认为站点不可达并暂停上传。
3. `--interpreter` 可以指定另一个解释器（例如装了 `swanlab==0.10.6` 的 venv），用来对比上一版正式发布。
4. `--keep` 保留子进程工作目录和 mock 日志，便于排查失败 case。

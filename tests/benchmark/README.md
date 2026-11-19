# tests/benchmark

SwanLab SDK 性能基准。目录按被测对象分层，方案与背景见
`playground/core-plan/bench-go-core-plan.md`。

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

## 运行

```bash
# 小规模验证（秒级），先确认链路通
uv run python tests/benchmark/e2e/bench_e2e.py --scenario C --steps 100 --keys 100 --core both

# 场景 A/B/C 正式跑（每个 case 重复 3 次取中位数）
uv run python tests/benchmark/e2e/bench_e2e.py --scenario A --core both --repeats 3
uv run python tests/benchmark/e2e/bench_e2e.py --scenario B --core both --repeats 3
uv run python tests/benchmark/e2e/bench_e2e.py --scenario C --core both --repeats 3
uv run python tests/benchmark/e2e/bench_e2e.py --scenario C --skip-store

# 给 mock 加 20ms 响应延迟，观察上传队列是否堆积
uv run python tests/benchmark/e2e/bench_e2e.py --scenario C --latency-ms 20

# 多进程（DDP 场景）
uv run python tests/benchmark/multi_process/bench_multi_process.py --ranks 8 --steps 100 --keys 1000
uv run python tests/benchmark/multi_process/bench_multi_process.py --ranks 16 --steps 100 --keys 1000 --sleep 0.2
```

`make bench` 等价于 `uv run pytest tests/benchmark`，只跑带 `test_*` 入口的历史用例；
`bench_e2e.py` 和 `bench_multi_process.py` 是命令行脚本，不参与 pytest 收集。

## 场景

| 场景 | 负载 | 说明 |
|---|---|---|
| A | 100,000 步 × 9 指标 = 90 万条 | 高频少指标（强化学习、小步迭代） |
| B | 7 步 × 100,000 指标 = 70 万条 | 低频大批量；只报 init/finish/内存/体积，不报延迟分位 |
| C | 1,000 步 × 100 指标 = 10 万条 | 常规规模，支持 `--skip-store` 对照 |

`--profile` 固定 `core.record_batch` 与 `core.record_interval`，两档结果不要放在同一张表里比：

| 档位 | batch | interval | 用途 |
|---|---|---|---|
| `default`（默认） | 10,000 | 1.5 | 接近真实使用，会上传中触发队列堆积 |
| `isolated` | 10,000 | 3,600 | 只在 finish 时统一上传，便于切分生产/收尾两段耗时 |

`record_batch` 会影响多个环节：Python `BackgroundConsumer` 的小批固定是 `min(100, batch)`、
gRPC 批是 `min(2000, batch)`，Go 侧队列容量是 `8 × batch`。两个 core 的发送时机也不一样
（Python 缓冲非空就发，Go 要攒满 batch 或等满一个 interval）。

## 校验项（任一不通过则本轮数据不可用）

`bench_e2e.py` 每个 case 都会检查：

| 检查 | 含义 |
|---|---|
| 实际使用的 core 与预期一致 | Go 二进制缺失时 SDK 会静默退回 Python core，必须挡住 |
| 子进程未连真实后端 | 回读 `Settings`，`api_host` / `api_key` / `on_online_error` 必须是注入的值 |
| Go core 进程数 / Python 无 core 进程 | 确认进程模型与预期一致 |
| mock 收到标量条数 | 应等于 `steps × keys`（按 record 类型分别计数） |
| 本地文件标量条数 | `DataStoreReader` 回读后应等于 `steps × keys`（`--skip-store` 时跳过） |
| 上传全程未被暂停 | `upload.paused` 不能出现 true；出现通常说明 mock 的健康检查没做对 |
| core 日志无报错 / 子进程无异常堆栈 | 排除失败但被吞掉的场景 |

多进程脚本额外检查「共享模式下只有 1 个 core 进程」「每个 rank 各启动 1 个」「各 rank 用的
core 进程互不相同」，依据来自子进程自己上报的 pid（父进程轮询可能错过短命进程）。

## 测量方式

- `run.log()` 延迟用 `perf_counter` 逐次采样，报 mean/P50/P95/P99（微秒）。它反映的是主线程
  是否被后台消费拖慢，而不是 core 的绝对吞吐。
- 内存为 `ps -o rss=` 每 100ms 采样，Python 进程和 Go core 分别记录再相加；相加会把共享页
  重复计算，只能当作上界。CPU 取 `ps -o time=` 的累计时间差，macOS 上精度是 1 秒。
- 结果写入 `results/*.json`，含 git sha、core 二进制 sha256、`GOGC` 和实际使用的参数。

## 注意事项

1. 每组开跑前会先清掉残留的 core 进程（`detach=True` 默认要空闲 300s 才退出）。残留实例会
   带着上一轮的指标状态和内存继续服务，会污染多进程结果。
2. Mock Cloud 的 `/api/health`、`/api/house/health` 和未注册路径都必须返回 200 或 4xx；
   返回 5xx 或断连会让 Go core 认为站点不可达并暂停上传。
3. `--interpreter` 可以指定另一个解释器（例如装了 `swanlab==0.10.6` 的 venv），用来对比上一版正式发布。
4. `--keep` 保留子进程工作目录和 mock 日志，便于排查失败 case。

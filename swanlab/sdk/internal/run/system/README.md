# System

SwanLab 系统信息采集模块，负责在实验运行期间采集三类信息：

| 模块 | 采集方式 | 内容 |
|------|----------|------|
| [metadata](./metadata) | 启动时采集一次（快照） | CPU/GPU 型号、OS、Python 版本、pip 包、git 信息等 |
| [monitor](./monitor) | 周期性采集（动态监控） | CPU 使用率、GPU 占用率、内存、磁盘 IO 等性能指标 |
| [console_proxy](./console_proxy) | 事件驱动（终端输出） | stdout / stderr 日志 |

---

## 核心约束：metadata 先于 monitor

`monitor` 模块中的**内置采集函数**不能在程序启动时静态确定，原因是：

- 是否启用 GPU 监控，取决于系统中是否存在 GPU
- 使用 NVML 还是 ROCm，取决于驱动类型
- 磁盘 IO 监控的挂载点，也依赖运行时环境

因此，内置 Monitor 是 **`metadata` 采集结果（`SystemInfo`）的函数**，必须在 metadata 采集完成后才能构建。这决定了整个模块的初始化顺序。

---

## 初始化两阶段

`SystemManager`（`system/__init__.py`）是唯一知道这一顺序的地方，其余模块彼此不感知：

```
Phase 1  metadata.collect()
              │
              ▼
         SystemInfo          ← 不可变快照，描述当前系统环境
              │
Phase 2  hardware.build_monitors(SystemInfo)
              │
              ▼
         内置 Monitor 列表    ← 按环境条件构建，缺失硬件对应的 Monitor 自动跳过
              │
         + monitors.get_all() ← 用户通过 merge_monitors() 注入的 Monitor（无条件追加）
              │
              ▼
         SystemManager._monitors
              │
              ▼
         BackgroundConsumer  ← 周期性调用所有 Monitor.fn()
```

### `hardware/` 模块的接口约定

每个硬件子模块暴露一个 `build(info: SystemInfo) -> Monitor | None` 工厂函数：

```python
# monitor/hardware/gpu.py
def build(info: SystemInfo) -> Monitor | None:
    if info.gpu_count == 0:
        return None  # 无 GPU，跳过
    return Monitor(id="gpu", fn=GPUCollector(info).collect)
```

`monitor/hardware/__init__.py` 汇总所有子模块，过滤掉 `None`：

```python
def build_monitors(info: SystemInfo) -> list[Monitor]:
    candidates = [gpu.build(info), cpu.build(info), disk.build(info), ...]
    return [m for m in candidates if m is not None]
```

---

## 模块依赖关系

```
metadata/  ──produces──▶  SystemInfo  ──consumed by──▶  monitor/hardware/

monitors.py (用户注入 store)  ──────────────────────────────────┤
                                                               ▼
                                                       SystemManager
```

三个子模块彼此解耦：

- `metadata/` 只管采集环境信息，**不知道 monitor 的存在**
- `monitor/hardware/` 只依赖 `SystemInfo` 这个数据类型，**不依赖 metadata 的内部实现**
- `monitors.py` 是独立的全局 store，**与 SystemInfo 完全无关**

---

## 用户注入接口

### `swanlab.merge_monitors(monitors)`

在 `swanlab.init()` 之前调用，向 `monitors.py` 的全局 store 注册自定义采集函数。
id 重复时覆盖已有项。这些 Monitor 在 Phase 2 末尾**无条件**追加到内置 Monitor 列表后，不经过系统环境检测。

```python
swanlab.merge_monitors([
    Monitor(id="npu_util", fn=lambda: {"util": npu.utilization()}),
])
```

---

## Settings 对应关系

`system/` 的行为受 `Settings` 中以下字段控制：

```
Settings.metadata   ← 控制内置 metadata 各采集项的开关（hardware/runtime/requirements/conda/git）
Settings.monitor    ← 控制周期性监控的开关及参数（enable/interval/disk_io_dir）
Settings.console    ← 控制终端日志代理的行为（proxy_type/max_log_length）
```

Settings 只管"采不采"，`system/` 只管"怎么采"，两者通过 `SystemManager.initialize()` 衔接。
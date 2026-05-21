"""
@author: cunyue
@description: CLI ping 模块，连通性与环境诊断命令

预期输出样式：

╭─ SwanLab Diagnostics ─────────────────────────────────────╮
│ Status: OK                                                │
│ Endpoint: https://api.swanlab.cn                          │
╰───────────────────────────────────────────────────────────╯

SDK
  ✓ Version        0.7.0
  ✓ Python         3.11.8
  ✓ Mode           cloud

Cloud
  ✓ Endpoint       https://api.swanlab.cn
  ✓ Status         OK
  ✓ Latency        42 ms
  ✓ Login          OK
  ✓ Username       cunyue

System
  ✓ OS             macOS 15.5 arm64
  ✓ Python         3.11.8
  ✓ Executable     /path/to/python

Hardware
  ✓ Apple Silicon  Apple M4 Pro
  ✓ CPU Cores      14
  ✓ Unified Memory 48 GB

  Accelerators
    - Not detected

硬件展示说明：
  - Accelerators 不要展示成单个 GPU 字段。probe_python 中加速器是 list[AcceleratorSnapshot]，
    需要支持多卡、多厂商异构场景。
  - 有加速器时按 vendor 分组，再列出每个 device，单卡也沿用同一结构，避免输出形状随设备数量变化。

异构加速器示例：
  Accelerators
    ✓ NVIDIA        2 devices, driver 550.54, CUDA 12.4
      [0] RTX 4090  24 GB
      [1] RTX 4090  24 GB
    ✓ Ascend        8 devices, CANN 8.0
      [0] Ascend 910B 64 GB

状态标识：
  ✓ 正常
  ! 警告
  ✗ 异常
  - 未检测到 / 未配置 / 不适用
"""

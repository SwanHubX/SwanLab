"""
@author: cunyue
@file: __init__.py
@time: 2026/4/16 21:12
@description: 上下文组件，存储swanlab run 运行生命周期所需要的外置组件，他们包括：

1. core: 核心组件，负责处理record的持久化和后端交互
2. system: 系统组件，负责系统环境信息采集、监控指标上报
3. callback: 回调组件，负责处理用户自定义的回调函数
4. config: 配置组件，负责处理用户上传的自定义配置
"""

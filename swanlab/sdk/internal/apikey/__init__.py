"""
@author: cunyue
@file: __init__.py
@time: 2026/3/6 12:53
@description: SwanLab API Key 管理，主要与本地存储相关
我们将环境变量、本地配置文件等来源在本模块抹平，方便业务调用，这里的优先级从低到高依次为：
1. 本地配置文件
2. 环境变量
3. context中内容
"""

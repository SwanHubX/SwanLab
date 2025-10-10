"""
@author: cunyue
@file: v1.py
@time: 2025/6/20 15:16
@description: swanlab data transfer protocol version 1
基于protobuf的数据模型将自动实现数据从python到go的转换，由于我们的数据在最开始基于swanlab.toolkit中的模型定义，因此我们需要在这里实现从toolkit模型到protobuf模型的转换
仅需几个函数实现字段映射即可
"""

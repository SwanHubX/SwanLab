<p align="center">
  <img alt="SwanLab Library" src="readme_files/swanlab-logo.svg" width="120" height="120">
</p>
<h1 align="center"><a href="https://github.com/SwanHubX/SwanLab/tree/main">SwanLab</a></h1>

<p align="center">
跟踪与可视化你的机器学习全流程
</p>

<p align="center">
  <a href="https://github.com/SwanHubX/SwanLab/stargazers"><img src="https://img.shields.io/github/stars/SwanHubX/SwanLab?style=social" alt= /></a>
  <a href="https://github.com/SwanHubX/SwanLab/blob/main/LICENSE"><img src="https://img.shields.io/github/license/SwanHubX/SwanLab.svg?color=brightgreen" alt="license"></a>
  <a href="https://github.com/SwanHubX/SwanLab/commits/main"><img src="https://img.shields.io/github/last-commit/SwanHubX/SwanLab" alt="license"></a>
  <a href="https://pypi.python.org/pypi/swanlab"><img src="https://img.shields.io/pypi/v/swanlab?color=orange" alt= /></a>
  <a href="https://pepy.tech/project/swanlab"><img alt="pypi Download" src="https://static.pepy.tech/badge/swanlab/month"></a>
  <a href="https://github.com/SwanHubX/SwanLab/discussions"><img alt="Github Discussion" src="https://img.shields.io/badge/discussions-GitHub-333333?logo=github"></a>
</p>

<p align="center">
👋 加入我们的<a href="https://geektechstudio.feishu.cn/wiki/NIZ9wp5LRiSqQykizbGcVzUKnic">微信</a>
</p>

<p align="center">
<a href="README.md">English</a> | 中文
</p>

## 目录

- [核心功能](#核心功能)
- [更新日志](#更新日志)
- [案例](#案例)
- [安装](#安装)
- [快速开始](#快速开始)
- [更多技巧](#更多技巧)
- [协议](#协议)

## 核心功能

- **🧪 表格视图**: 对比关键指标，更快获得洞见

<div align="center">
  <img src="readme_files/experiments-gridView.gif" width="600">
</div>

- **📊 图表视图**: 可视化你的机器学习训练全过程

<div align="center">
  <img src="readme_files/charts-1.gif" width="600">
</div>

- **⚡️ 跟踪机器学习流程**: 自动保存超参数，配置，度量指标，终端日志记录，环境信息

<div align="center">
  <img src="readme_files/track-machine-learning-pipeline.gif" width="600">
</div>


- **🥔 离线支持**: SwanLab可以完全离线运行，无需任何对互联网的访问。例如，它可以在您的本地计算机上、企业防火墙后面或数据中心中运行。

> [!NOTE]
> 请查看SwanLab的<a href="https://geektechstudio.feishu.cn/wiki/MwXmw9yDeiZWyQkPnNgcixwWnwu">文档</a>和<a href="https://github.com/SwanHubX/SwanLab-examples">示例</a>，以获取有关SwanLab的完整介绍。

<br>

## 更新日志

[24/01/23] 🚨 我们使用SQLite数据库和Peewee库替代了之前的基础配置信息读写方案（[#114](https://github.com/SwanHubX/SwanLab/issues/114)），这是个极大有利于项目未来的改动，但缺陷是不兼容旧版本（swanlab<=v0.1.4）的日志数据文件，所以如需可视化旧版本产生的日志文件请安装<=0.1.4的版本。与此同时，我们增加了支持导出实验列表为CSV，新的环境记录项`Run path`和`logdir`，增加了快捷复制的交互，以及新的API `swanlab.config`。

[24/01/14] 🔥 我们发布了一个新的UI界面，以及跟踪更多环境信息，包括Command、git提交/分支、和机器内存。此外，我们还添加了一个`logdir` API，允许开发人员设置日志文件的目录。

[24/01/07] ✨ 我们支持在仪表板上删除实验和编辑实验信息。

[完整更新日志](https://github.com/SwanHubX/SwanLab/releases)

<br>

## 案例

通过以下用例学习如何更有效地使用SwanLab：

| 案例 | 介绍 | 
| ------- | ------- |
| [Hello World](https://github.com/SwanHubX/SwanLab-examples/tree/main/Hello_World) | 简单入门 |
| [MNIST](https://github.com/SwanHubX/SwanLab-examples/tree/main/MNIST) | 基于神经网络的MNIST手写体识别（使用pytorch、swanlab库） |
| [图像分类](https://github.com/SwanHubX/SwanLab-examples/blob/main/Resnet50) | ResNet50猫狗分类（使用pytorch、swanlab、Gradio库） [图文教程](https://zhuanlan.zhihu.com/p/676430630) |
| [文本生成](https://github.com/SwanHubX/SwanLab-examples/blob/main/Word_language_model) | 基于自然语言模型的文本生成 (RNN/LSTM/GRU/Transformer) |

<br>

## 安装

### pip安装

环境要求：Python 3.8+

使用[pip](https://pip.pypa.io/en/stable/)将安装我们稳定发布的版本，安装命令如下所示：

```bash
pip install -U swanlab
```

### 源码安装

如果您等不及发布，想体验最新的代码与特性，那么必须[从源代码安装此库](https://geektechstudio.feishu.cn/wiki/DvxSweHUKiAe8yksci3cMflbnwh#SMXHdJ1c1o4jzTxcDticHcwvnHd)。

<br>

## 快速开始

1. 首先，使用[pip](https://pip.pypa.io/en/stable/)安装SwanLab SDK:

```bash
pip install -U swanlab
```

2. 然后，使用下面的示例代码片段作为模板，将SwanLab集成到您的Python脚本中:
```Python
import swanlab

# Start a SwanLab Run with swanlab.init
swanlab.init(
  # save model inputs and hyperparameters in a swanlab.config object
  config={'learning_rate': 0.01},
)

# Model training code here...

# Log metrics over time for visualizing performance with swanlab.log
for epoch in range(1, 20):
    swanlab.log({"loss": loss})
```

3. 第三步，开启一个SwanLab仪表板: 
```bash
$ swanlab watch
```

就是这样！打开 http://127.0.0.1:5092 ，查看你的第一个SwanLab实验的仪表板。

<br>

## 更多技巧

🏄‍♀️ 实验看板相关：

<details>
<summary>设置日志文件保存路径，并运行基于它的实验看板</summary>

设置日志文件保存路径，比如`./logs`: 

```python
import swanlab 

swanlab.init(
  logdir="./logs"
)
```

运行基于`./logs`的看实验看板：

```bash
$ swanlab watch --logdir ./logs
```
</details>


<details>
<summary>设置实验看板的主机和端口</summary>

```bash
$ swanlab watch --host 0.0.0.0 --port 8080
```
</details>

<details>
<summary>远程访问实验看板</summary>

- 链接：[在远程服务器上进行训练时，访问SwanLab实验看板](https://zhuanlan.zhihu.com/p/677224865)

</details>

⚙️ 其他：

<details>
<summary>argparse与swanlab.config</summary>

swanlab.config支持直接传入argparse.Namespace类型的变量，如:
```python
import argparse
import swanlab

parser = argparse.ArgumentParser()
···
args = parser.parse_args()


swanlab.init(
    config=args
)
```
</details>

<br>

## 协议

此项目当前的许可证协议是 [Apache 2.0 License](https://github.com/SwanHubX/SwanLab/blob/main/LICENSE).

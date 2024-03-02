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
  <img alt="SwanLab Head Image" src="readme_files/swanlab-head-image.png" width="800">
</p>


<p align="center">
👀 查看<a href="https://geektechstudio.feishu.cn/wiki/MwXmw9yDeiZWyQkPnNgcixwWnwu">官方文档</a> | 👋 加入我们的<a href="https://geektechstudio.feishu.cn/wiki/NIZ9wp5LRiSqQykizbGcVzUKnic">微信交流群</a>
</p>


<p align="center">
<a href="README.md">English</a> | 中文
</p>

## 目录

- [更新日志](#更新日志)
- [核心功能](#核心功能)
- [安装](#安装)
- [快速开始](#快速开始)
- [使用教程](#使用教程)
- [案例](#案例)
- [协议](#协议)

<br>

## 更新日志

升级到最新版本: `pip install -U swanlab`。

[24/03/01] 🚀 依旧是超大杯的更新！我们支持了[文本图表](https://geektechstudio.feishu.cn/wiki/T0L7wYfzGiZUCKkxfehcFwYAnIh)以适配NLP、LLM、Agent等场景任务的需求; 我们对折线图的UI、图例、渲染速度做了大量优化，并提高了Logs的渲染性能，200k行的终端打印信息查看也不卡顿。（v0.2.1）

[24/02/08] 🔥 超大更新! 我们支持了[图像图表](https://geektechstudio.feishu.cn/wiki/LZFxwTuegiXxPGkhXcpcBUEXnHb)、[音频图表](https://geektechstudio.feishu.cn/wiki/SU6mwcVNbixMf1k95KbcZHDCnJe)、多实验图表以及一系列全面的优化和改进！可通过 `pip install -U swanlab` 升级到最新版本体验新特性。（v0.2.0）

[24/01/25] 😄 我们发布了新的Config/Summary表格组件，支持参数搜索。此外我们还使用了新的字体和配色。（v0.1.6）

[完整更新日志](https://github.com/SwanHubX/SwanLab/releases)

<br>

## 核心功能

- **📊 训练可视化**: 可视化你的机器学习训练全过程

<div align="center">
  <img src="readme_files/charts-1.gif" width="600">
</div>

- **🚀 多媒体图表**: 记录训练中的图像/音频/视频/文本/3D模型...

<div align="center">
  <img src="readme_files/mutilmedia-chart.gif" width="600">
</div>


- **🧪 表格视图**: 对比关键指标，更快获得洞见

<div align="center">
  <img src="readme_files/experiments-table.png" width="600">
</div>


- **⚡️ 自动保存环境信息**: 自动保存超参数，配置，指标，终端日志记录，pip环境信息等

- **🥔 离线支持**: SwanLab可以完全离线运行，无需任何对互联网的访问。例如，它可以在您的本地计算机上、企业防火墙后面或数据中心中运行。

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
```python
import swanlab

# 初始化swanlab
swanlab.init(
  config={'epochs': 20, 'learning_rate': 0.01},  # 通过config参数保存输入或超参数
  logdir="./logs",  # 指定日志文件的保存路径
)

# 把模型训练的代码放到这里...
...

# 使用swanlab.log记录指标变化的数据
for epoch in range(1, swanlab.config.epoch):
    swanlab.log({"loss": loss})
```

例如, 我们写1个模拟实验脚本:
```python
import swanlab
import random

offset = random.random() / 5

run = swanlab.init(
    experiment_name="Example",
    description="这是一个机器学习模拟实验",
    config={
        "learning_rate": 0.01,
        "epochs": 20,
    },
    logdir="./logs"
)

# 模拟机器学习训练过程
for epoch in range(2, run.config.epochs):
    acc = 1 - 2**-epoch - random.random() / epoch - offset
    loss = 2**-epoch + random.random() / epoch + offset
    swanlab.log({"loss": loss, "accuracy": acc})
```

3. 最后，打开终端，使用下面的指令，开启一个SwanLab仪表板: 
```bash
$ swanlab watch -l ./logs
```

运行完成后，SwanLab会给你1个URL链接（默认是http://127.0.0.1:5092），查看链接，即可在浏览器看到你的第一个实验可视化结果。

<div align="center">
  <img src="readme_files/get-started.png" width="600">
</div>

<br>

## 使用教程

**入门教程**
- [安装](#安装)
- [快速上手](https://geektechstudio.feishu.cn/wiki/UInBw9eaziv17IkwfrOcHCZ1nbc)
- [启动实验看板](https://geektechstudio.feishu.cn/wiki/YsEfwC79viJL2nk5TgPcAOUhn5U)  

**Python API**
- [init](https://geektechstudio.feishu.cn/wiki/H7Wbwt91LiCJtnkpHOzcar4TnCc)
- [log](https://geektechstudio.feishu.cn/wiki/RmjXwjmgUi5zGCkBPsTc5ygQn4g)
- [config](https://geektechstudio.feishu.cn/wiki/HkTOwxLkHiUC84kJNrlcohyGnuh)
- [Image - 图像图表](https://geektechstudio.feishu.cn/wiki/LZFxwTuegiXxPGkhXcpcBUEXnHb)
- [Audio - 音频图表](https://geektechstudio.feishu.cn/wiki/SU6mwcVNbixMf1k95KbcZHDCnJe)
- [Text - 文本图表](https://geektechstudio.feishu.cn/wiki/T0L7wYfzGiZUCKkxfehcFwYAnIh)

**CLI API**
- [watch - 开启实验看板](https://geektechstudio.feishu.cn/wiki/Q6I5wdyr9iRYkdkZ2gYcHQkxnCU)

**技巧**
- [远程访问实验看板](https://geektechstudio.feishu.cn/wiki/Icesw6coTidDsPkN960c0lNtnCb)
- [将argparse传入swanlab.config](https://geektechstudio.feishu.cn/wiki/CT1Xwo6ehimNH5kz7y9csTGkn0e)

<br>

## 案例

通过以下用例学习如何更有效地使用SwanLab：

| 案例 | 介绍 | 
| ------- | ------- |
| [Hello World](https://github.com/SwanHubX/SwanLab-examples/tree/main/Hello_World) | 简单入门 |
| [MNIST](https://github.com/SwanHubX/SwanLab-examples/tree/main/MNIST) | 基于神经网络的MNIST手写体识别（使用pytorch、swanlab库） |
| [图像分类](https://github.com/SwanHubX/SwanLab-examples/blob/main/Resnet50) | ResNet50猫狗分类（使用pytorch、swanlab、Gradio库） [图文教程](https://zhuanlan.zhihu.com/p/676430630) |
| [文本生成](https://github.com/SwanHubX/SwanLab-examples/blob/main/Word_language_model) | 基于自然语言模型的文本生成 (RNN/LSTM/GRU/Transformer) |
| [微调UIE](https://github.com/SwanHubX/SwanLab-examples/tree/main/UIE) | 如何使用个人数据来微调UIE模型并通过swanlab监控训练过程 |

<br>

## 协议

此项目当前的许可证协议是 [Apache 2.0 License](https://github.com/SwanHubX/SwanLab/blob/main/LICENSE).

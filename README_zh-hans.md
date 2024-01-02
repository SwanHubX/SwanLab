<p align="center">
  <img alt="SwanLab Library" src="readme_files/swanlab-logo.svg" width="120" height="120">
</p>

<h1 align="center"><a href="https://github.com/SwanHubX/SwanLab/tree/main">SwanLab</a></h1>

<p align="center">SwanLab是一个强大的开源机器学习训练管理工具，供研究人员跟踪记录自己的训练。通过使用SwanLab，研究人员可以积累训练经验并发现新的Idea。</p>

<p align="center">
  <b><a href="README.md">English</a> | 简体中文</b>
</p>
<p align="center">
  <a href="https://pypi.python.org/pypi/swanlab"><img src="https://img.shields.io/pypi/v/swanlab?color=blue" alt= /></a>
  <a href="https://pepy.tech/project/swanlab"><img alt="pypi Download" src="https://static.pepy.tech/badge/swanlab/month"></a>
  <a href="https://github.com/SwanHubX/SwanLab/discussions"><img alt="Github Discussion" src="https://img.shields.io/badge/discussions-GitHub-333333?logo=github"></a> 
  <a href="https://geektechstudio.feishu.cn/wiki/space/7310593325374013444?ccm_open_type=lark_wiki_spaceLink&open_tab_from=wiki_home"><img alt="Website" src="https://img.shields.io/badge/website-online-green"></a>
  <a href="https://github.com/SwanHubX/SwanLab/blob/main/LICENSE"><img src="https://img.shields.io/github/license/SwanHubX/SwanLab.svg?color=brightgreen" alt="license"></a>
</p>

<img alt="hello_world_main2" src="readme_files/hello_world_main2.gif" width=1535>

<br>

## ✨ 特色功能

1. ⚽ **实时指标记录：** 几行代码，即可实时记录你的训练指标
2. 🧪 **多实验对比**：支持多实验指标对比
3. 🤖 **机器学习支持：** 支持 PyTorch、TensorFlow、Transformers、mmdetection 等主流训练框架
4. 📝 **环境记录**：支持自动记录 Logging、报错、系统硬件、Python 环境等等一系列环境信息
5. 🖥 **端云均支持：** 即支持本地管理训练，也支持同步到公有云（即将）

<br>

## 🔥 使用案例

我们提供了一些案例代码和文章，来帮助你更好地理解与掌握 SwanLab：

- [Hello World](https://github.com/SwanHubX/SwanLab-examples/tree/main/Hello_World)
- [MNIST 手写体识别](https://github.com/SwanHubX/SwanLab-examples/tree/main/plain_net_mnist)
- [ResNet50 猫狗分类](https://github.com/SwanHubX/SwanLab-examples/tree/main/resnet50_cats_vs_dogs)

<br>

## ⚡️ 快速上手

Hi，无论你是开发人员还是日常用户，这篇快速上手教程都将帮助你入门并且向你展示如何使用 SwanLab：

- 记录训练配置信息
- 记录关键指标
- 可视化实验

<br>

### 🎯 第 1 步：安装 SwanLab

```bash
$ pip install -U swanlab
```

<br>

### 👋 第 2 步：Hello World

如果我们抽象机器学习的训练过程，其本质就是**配置参数**、再**循环**的过程，而我们关注的是中间的**指标**。

下面的 Python 代码模拟了这一过程：

```Python
import swanlab

# 初始化swanlab
swanlab.init()

for epoch in range(1, 20):
    print("epoch", epoch)
    # 跟踪指标epoch
    swanlab.log({"epoch": epoch})
```

其中，`swanlab.init`是必需的，作用是初始化实例以及配置参数；`swanlab.log`的作用是负责记录数据，接收的数据类型是 1 个字典（dict）。

运行上面的代码，你会看到下面的输出结果：

```Bash
[SwanLab-INFO]:        Run data will be saved locally in path/swanlog/majestic-hemlock-1
[SwanLab-INFO]:        Experiment_name: majestic-hemlock-1
[SwanLab-INFO]:        Run `swanlab watch` to view SwanLab Experiment Dashboard
epoch 1
epoch 2
epoch 3
epoch 4
epoch 5
epoch 6
epoch 7
epoch 8
epoch 9
[SwanLab-INFO]:        train successfully
```

并且根目录下会出现 1 个`swanlog`文件夹，里面是 SwanLab 自动生成的文件，记录了一系列实验数据。

<br>

### 🧪 第 3 步：开启实验看板

现在来查看我们使用 SwanLab 记录的指在每个循环步骤中的情况。运行命令`swanlab watch`：

```Bash
$ swanlab watch

[SwanLab-INFO]:        SwanLab Experiment Dashboard ready in 375ms
                       ➜  Local:   http://127.0.0.1:5092
```

访问http://127.0.0.1:5092 ，打开实验看板，访问刚刚运行的实验：

<img src="readme_files/hello_world_main1.gif" width=1535>

<br>

### 🚀 进阶一下

在这一节，让我们写 1 个进阶的训练脚本来模拟真实的机器学习训练。

首先初始化 swanlab，这次设置了实验名称、介绍和配置：

```Python
swanlab.init(
    # 设置实验名称
    experiment_name="sample_experiment"
    # 设置实验介绍
    description="This is a sample experiment for machine learning training.",
    # 记录跟踪的超参数和运行元数据
    config={
        "learning_rate": lr,
        "epochs": epochs,
    },
)
```

组合成 1 个完整的训练脚本，使用`swanlab.log`API 追踪损失值`loss`和准确率`accuracy`：

```Python
import swanlab
import time
import random

lr = 0.01
epochs = 20
offset = random.random() / 5

swanlab.init(
    # 设置实验名称
    experiment_name="sample_experiment",
    # 设置实验介绍
    description="This is a sample experiment for machine learning training.",
    # 记录跟踪的超参数和运行元数据
    config={
        "learning_rate": lr,
        "epochs": epochs,
    },
)

# 模拟机器学习训练过程
for epoch in range(2, epochs):
    acc = 1 - 2**-epoch - random.random() / epoch - offset
    loss = 2**-epoch + random.random() / epoch + offset
    print(f"epoch={epoch}, accuracy={acc}, loss={loss}")
    # 记录loss和acc
    swanlab.log({"loss": loss, "accuracy": acc})
    time.sleep(1)
```

同样的，运行`swanlab watch`启动实验看板：

<img alt="hello_world_main2" src="readme_files/hello_world_main2.gif" width=1535>

<br>

## 🌱 了解更多

- [官方文档](https://geektechstudio.feishu.cn/wiki/space/7310593325374013444?ccm_open_type=lark_wiki_spaceLink&open_tab_from=wiki_home)：完整的 APl 文档与引导。
- [案例仓库](https://github.com/SwanHubX/SwanLab-examples)：官方代码案例

<br>

## 💬 社区

加入 SwanLab 社区，分享您的想法、建议或问题，并与其他用户和贡献者交流。

WeChat 与 Github 社区入口：

[![PyPI - Downloads](https://img.shields.io/badge/wechat-online-green)](https://geektechstudio.feishu.cn/wiki/NIZ9wp5LRiSqQykizbGcVzUKnic)[![Discuss on GitHub](https://img.shields.io/badge/discussions-GitHub-333333?logo=github)](https://github.com/SwanHubX/SwanLab/discussions)

<br>

## 版权说明

SwanLab 是 [Apache 2.0 License](https://github.com/SwanHubX/SwanLab/blob/main/LICENSE)。

<p align="center">
  <img alt="SwanLab Library" src="readme_files/swanlab-logo-light.svg#gh-light-mode-only" width="352" height="59">
  <img alt="SwanLab Library" src="readme_files/swanlab-logo-dark.svg#gh-dark-mode-only" width="352" height="59">
  <br/>
  <br/>
</p>

<p align="center">
  <a href="https://pypi.python.org/pypi/swanlab"><img src="https://img.shields.io/pypi/v/swanlab" alt= /></a>
  <a href="https://pepy.tech/project/swanlab"><img alt="pypi Download" src="https://static.pepy.tech/badge/swanlab/month"></a>
  <a href="https://geektechstudio.feishu.cn/wiki/space/7310593325374013444?ccm_open_type=lark_wiki_spaceLink&open_tab_from=wiki_home"><img alt="Website" src="https://img.shields.io/badge/website-online-blue"></a>
  <a href="https://github.com/SwanHubX/SwanLab/blob/main/LICENSE"><img src="https://img.shields.io/github/license/SwanHubX/SwanLab.svg" alt="license"></a>
  <a href="https://github.com/SwanHubX/SwanLab/releases"><img alt="GitHub release" src="https://img.shields.io/github/release/SwanHubX/SwanLab.svg"></a>
</p>

<p align="center">
  <b><a href="https://github.com/SwanHubX/SwanLab/blob/main/README.md">English</a> | 简体中文</b>
</p>

SwanLab 是[SwanHub](https://swanhub.co)开源社区发布的新一代机器学习实验管理与可视化工具，旨在让机器学习训练有效地协作起
来。

SwanLab 提供简洁的 API，轻松实现机器学习指标跟踪与配置记录。同时，SwanLab 还提供了一个可视化看板，以最直观的方式**监看、
分析和对比**你的训练。

有关 SwanLab 功能的具体指南，请参
阅[用户指南](https://geektechstudio.feishu.cn/wiki/space/7310593325374013444?ccm_open_type=lark_wiki_spaceLink&open_tab_from=wiki_home)。

目前，SwanLab 正在快速迭代，并将持续添加新功能。

## 安装

此存储库在 Python 3.8+上进行了测试。

SwanLab 可以使用 pip 安装，如下所示:

```bash
pip install swanlab
```

## 快速开始

让我们模拟一个简单的机器学习训练过程，使用`swanlab.init`来初始化实验并记录配置信息，并使用`swanlab.log`跟踪关键指标（在
这个例子中是 `loss` 和 `acc`）：

```python
import swanlab
import time
import random

lr = 0.01
epochs = 20
offset = random.random() / 5

# Initialize the experiment and record configuration information
swanlab.init(
    description="This is a sample experiment for machine learning training.",
    config={
        "learning_rate": lr,
        "epochs": epochs,
    },
)

# Simulate a machine learning training process
for epoch in range(2, epochs):
    acc = 1 - 2**-epoch - random.random() / epoch - offset
    loss = 2**-epoch + random.random() / epoch + offset
    print(f"epoch={epoch}, accuracy={acc}, loss={loss}")
    # Track key metrics
    swanlab.log({"loss": loss, "accuracy": acc})
    time.sleep(1)
```

在程序运行过程中，目录下会生成一个`swanlog`文件夹，记录了你的训练数据。

如果要可视化你的实验，那么打开终端，进入根目录（不必进入`swanlog`文件夹），运行如下命令：

```bash
swanlab watch
```

看见如下输出则表示实验看板运行成功：

```console
[SwanLab-INFO]:        SwanLab Experiment Dashboard ready in 375ms

                        ➜  Local:   http://127.0.0.1:5092
```

此时访问`http://127.0.0.1:5092`，即可进入实验看板以浏览你的实验结果：

<img alt="swanlab-dashboard-1" src="readme_files/swanlab-dashborad-1.png" width="800">

## 了解更多

| Section                                                                                                                           | Details                 |
| --------------------------------------------------------------------------------------------------------------------------------- | :---------------------- |
| [文档](https://geektechstudio.feishu.cn/wiki/space/7310593325374013444?ccm_open_type=lark_wiki_spaceLink&open_tab_from=wiki_home) | 完整的 APl 文档与教程。 |
| [案例](https://github.com/SwanHubX/SwanLab-examples)                                                                              | 官方代码案例            |

## License

[Apache 2.0 License](https://github.com/SwanHubX/SwanLab/blob/main/LICENSE)

## Star History

[![Star History Chart](https://api.star-history.com/svg?repos=swanhubx/swanlab&type=Date)](https://star-history.com/#swanhubx/swanlab&Date)

## 贡献者

[![swanlab contributors](https://contrib.rocks/image?repo=swanhubx/swanlab&max=2000)](https://github.com/SwanHubX/SwanLab/graphs/contributors)

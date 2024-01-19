<p align="center">
  <img alt="SwanLab Library" src="readme_files/swanlab-logo.svg" width="120" height="120">
</p>
<h1 align="center"><a href="https://github.com/SwanHubX/SwanLab/tree/main">SwanLab</a></h1>

<p align="center">
Track and visualize all the pieces of your machine learning pipeline
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
👋 Join our <a href="https://geektechstudio.feishu.cn/wiki/NIZ9wp5LRiSqQykizbGcVzUKnic">WeChat</a>
</p>

<p align="center">
English | <a href="README_zh-hans.md">中文</a>
</p>

## Table of Contents

- [Key Function](#key-function)
- [Changelog](#changelog)
- [Use Case](#use-case)
- [Getting Started](#getting-started)
- [More Tips](#more-tips)
- [LICENSE](#license)

## Key Function

- **🧪 Experiments GridView**: compare your key metrics for inspiration faster

<div align="center">
  <img src="readme_files/experiments-gridView.gif" width="600">
</div>

- **📊 Charts**: visualize your entire training process

<div align="center">
  <img src="readme_files/charts-1.gif" width="600">
</div>

- **⚡️ Track machine-learning pipeline**: Hyperparameters, Config, Metric, Terminal logging, Environment Information auto save

<div align="center">
  <img src="readme_files/track-machine-learning-pipeline.gif" width="600">
</div>


- **🥔 Offline Support**: SwanLab can run entirely offile, ithout requiring any access to the Internet. For instance, this may be on your local machine, behind a corporate firewall, or in a datacenter

> [!NOTE]
> See the SwanLab <a href="https://geektechstudio.feishu.cn/wiki/MwXmw9yDeiZWyQkPnNgcixwWnwu">Documentaion</a> and <a href="https://github.com/SwanHubX/SwanLab-examples">Examples</a> for a full description of the SwanLab.

<br>

## Changelog

[24/01/14] 🔥 We supported a new UI, tracking additional environment information, including command, git commit/branch and memory. Additionally, we've added a `logdir` API, allowing developers to set the directory for log files.

[24/01/07] ✨ We supported delete experiment and edit experiment inforamation on Dashboard.

[24/01/01] We fixed some bugs to make SwanLab more stable.

[Full Changelog](https://github.com/SwanHubX/SwanLab/releases)

<br>

## Use Case

Learn how to use SwanLab more effectively by following these use cases:

| Code Cases | Description | 
| ------- | ------- |
| [Hello World](https://github.com/SwanHubX/SwanLab-examples/tree/main/Hello_World) | Getting Started |
| [MNIST](https://github.com/SwanHubX/SwanLab-examples/tree/main/MNIST) | Handwriting recognition based on a plain net and MNIST dataset with pytroch, swanlab. |
| [Image Classification](https://github.com/SwanHubX/SwanLab-examples/blob/main/Resnet50) | Cat and dog classification based on ResNet50 with pytorch, swanlab and gradio. [Tutorial](https://zhuanlan.zhihu.com/p/676430630). |
| [Text Generation](https://github.com/SwanHubX/SwanLab-examples/blob/main/Word_language_model) | Text generation based on Word_language_model (RNN/LSTM/GRU/Transformer) |

<br>

## Getting Started

1. First, install the SwanLab SDK with [pip](https://pip.pypa.io/en/stable/):

```bash
pip install -U swanlab
```

2. Second, Use the example code snippet below as a template to integrate SwanLab to your Python script:
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

3. Third, Run a Dashboard: 
```bash
$ swanlab watch
```

That's it! Open http://127.0.0.1:5092 to view a dashboard of your first SwanLab Experiment.

<br>

## More Tips

- Set a log directory save path and run the Dashboard using it:
```python
import swanlab 

swanlab.init(
  logdir="./logs"
)
```

```bash
$ swanlab watch --logdir ./logs_path
```

- Set the Host and Port for the Dashboard: 
```bash
$ swanlab watch --host 0.0.0.0 --port 8080
```
- Use Argparse init swanlab.config: 
```python
import argparse
import swanlab

parser = argparse.ArgumentParser()
···
args = parser.parse_args()


swanlab.init(
    config=vars(args)
)
```

- [Remotely access Dashboard](https://zhuanlan.zhihu.com/p/677224865): Access the SwanLab Dashboard While Training on a Remote Server.

<br>

## LICENSE

This project is currently licensed under [Apache 2.0 License](https://github.com/SwanHubX/SwanLab/blob/main/LICENSE).

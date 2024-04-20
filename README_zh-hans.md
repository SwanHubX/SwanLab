![Frame 2545](/Users/zeyilin/Downloads/Frame 2545.png)

<p align="center">
<a href="https://geektechstudio.feishu.cn/wiki/MwXmw9yDeiZWyQkPnNgcixwWnwu">SwanLab在线版</a> · <a href="https://geektechstudio.feishu.cn/wiki/MwXmw9yDeiZWyQkPnNgcixwWnwu">文档</a> · <a href="https://geektechstudio.feishu.cn/wiki/NIZ9wp5LRiSqQykizbGcVzUKnic">微信</a>
</p>

<p align="center">
  <a href="https://swanlab.cn" target="_blank">
        <img alt="Static Badge" src="https://img.shields.io/badge/Product-F04438"></a>
  <a href="https://github.com/SwanHubX/SwanLab/stargazers"><img src="https://img.shields.io/github/stars/SwanHubX/SwanLab?style=social" alt= /></a>
  <a href="https://github.com/SwanHubX/SwanLab/blob/main/LICENSE"><img src="https://img.shields.io/github/license/SwanHubX/SwanLab.svg?color=brightgreen" alt="license"></a>
  <a href="https://github.com/SwanHubX/SwanLab/commits/main"><img src="https://img.shields.io/github/last-commit/SwanHubX/SwanLab" alt="license"></a>
  <a href="https://pypi.python.org/pypi/swanlab"><img src="https://img.shields.io/pypi/v/swanlab?color=orange" alt= /></a>
  <a href="https://pepy.tech/project/swanlab"><img alt="pypi Download" src="https://static.pepy.tech/badge/swanlab/month"></a>
  <a href="https://github.com/SwanHubX/SwanLab/discussions"><img alt="Github Discussion" src="https://img.shields.io/badge/discussions-GitHub-333333?logo=github"></a>
  <a href="https://github.com/swanhubx/swanlab/issues"><img alt="issues" src="https://img.shields.io/github/issues/swanhubx/swanlab"></a>
</p>


<div align="center">
  <a href="./README.md"><img alt="英文文档" src="https://img.shields.io/badge/英文-d9d9d9"></a>
  <a href="./README_zh-hans.md"><img alt="中文文档" src="https://img.shields.io/badge/简体中文-d9d9d9"></a>
</div>



SwanLab是一款开源的AI实验管理平台。其简洁的API和UI结合了超参数跟踪、指标记录、多人协同训练、实验链接分享、实时消息通知等功能，让您可以快速跟踪ML实验、可视化过程、分享给同伴，旨在加速AI研发团队100倍的研发效率。

我们相信...

SwanLab提供了一套UI，用于探索和比较AI训练实验。以下是其核心功能列表：

1. 123
2. 234
3. 345
4. 



| 在你的机器学习流程中记录元数据 | 通过UI可视化和比较元数据 |
| ------------------------------ | ------------------------ |
|                                |                          |
| **有效地进行机器学习训练**     | **组织你的实验**         |
|                                |                          |



# 🎬 Demo







# 🏁 快速开始

请按照以下步骤开始使用 SwanLab。

## 1. 在你的环境中安装SwanLab

```bash
pip install -U swanlab
```



## 2. 登录SwanLab

...



ps：如果你想要离线使用，[自托管教程](#🏁 快速开始-自托管)



## 3. 将SwanLab与你的代码集成

```python
import swanlab

# 初始化一个新的swanlab实验
swanlab.init(
  # 记录超参数
  config={
    'batch_size': 32,
    'learning_rate': 0.01
  },  
  # 指定日志文件的保存路径
  logdir="./logs",  
)

# 记录指标
for i in range(10):
    swanlab.log({"loss": loss})
    swanlab.log({"acc": loss})
```



### 4. 在训练的同时启动SwanLab UI

打开终端，使用下面的指令，开启一个SwanLab仪表板: 

```bash
swanlab watch -l ./logs
```

运行完成后，SwanLab会给你1个URL链接（默认是http://127.0.0.1:5092），查看链接，即可在浏览器看到你的第一个实验可视化结果。

<div align="center">
  <img src="readme_files/get-started.png" width="600">
</div>
<br>



# 🏁 快速开始-自托管







# ⚡️更新日志

[24/04/...]  ...

[24/03/19] 🔧 我们修复了一些问题，并为即将到来的云端版本做准备。（v0.2.4）

[24/03/12] 👽 我们发布了折线图平滑功能，支持3种不同的平滑算法; 我们支持了在多实验图表中对比图像、音频图表; 同时改进了`swanlab.Image`，支持tensor作为输入。（v0.2.3）

[24/03/03] 🔧 我们修复了v0.2.1的一些问题，以及支持了通过按键切换多媒体图表的内容。（v0.2.2）

[24/03/01] 🚀 依旧是超大杯的更新！我们支持了[文本图表](https://geektechstudio.feishu.cn/wiki/T0L7wYfzGiZUCKkxfehcFwYAnIh)以适配NLP、LLM、Agent等场景任务的需求; 我们对折线图的UI、图例、渲染速度做了大量优化，并提高了Logs的渲染性能，200k行的终端打印信息查看也不卡顿。（v0.2.1）

[完整更新日志](https://github.com/SwanHubX/SwanLab/releases)

<br>



# 🆚与熟悉的工具的比较

SwanLab的设计立足在巨人的肩膀，感谢Tensorboard、Weights and Biases、MLFlow等工具在ML工具生态的巨大贡献！  

**Tensorboard vs SwanLab**

- **☁️支持在线使用**：通过SwanLab可以方便地将训练实验在云端在线同步与保存，便于远程查看训练进展、管理历史项目、分享实验链接、发送实时消息通知、多端看实验等。而Tensorboard是一个离线的实验跟踪工具。
- **👥多人协作**：在进行多人、跨团队的机器学习协作时，通过SwanLab可以轻松管理多人的训练项目、分享实验链接、跨空间交流讨论。而Tensorboard主要为个人设计，难以进行多人协作和分享实验。
- **💻持久、集中的仪表板**：无论你在何处训练模型，无论是在本地计算机上、在实验室集群还是在公有云的GPU实例中，你的结果都会记录到同一个集中式仪表板中。而使用TensorBoard需要花费时间从不同的机器复制和管理 TFEvent文件。
- **💪更强大的表格**：通过SwanLab表格可以查看、搜索、过滤来自不同实验的结果，可以轻松查看数千个模型版本并找到适合不同任务的最佳性能模型。 TensorBoard 不适用于大型项目。



**Weights and Biases vs SwanLab**

云端 vs 自托管+云端

- Weights and Biases 是一个必须联网使用的闭源MLOps平台
- SwanLab 不仅支持联网使用，也支持开源、免费、自托管的版本



# 👥 社区

## 社区与支持

- [GitHub Issues](https://github.com/SwanHubX/SwanLab/issues)：使用SwanLab时遇到的错误和问题
- [电子邮件支持](zeyi.lin@swanhub.co)：反馈关于使用SwanLab的问题
- <a href="https://geektechstudio.feishu.cn/wiki/NIZ9wp5LRiSqQykizbGcVzUKnic">微信交流群</a>：交流使用SwanLab的问题、分享最新的AI技术



## SwanLab README徽章

如果你喜欢在工作中使用 SwanLab，请将 SwanLab 徽章添加到你的README中：

[![swanlab](https://img.shields.io/badge/powered%20by-SwanLab-%23#b2d3bb)](https://github.com/swanhubx/swanlab)

```
[![swanlab](https://img.shields.io/badge/powered%20by-SwanLab-%23#b2d3bb)](https://github.com/swanhubx/swanlab)
```



## 在论文中引用SwanLab

如果您发现 SwanLab 对您的研究之旅有帮助，如果您能认可 SwanLab 的贡献，我们将非常开心：

```bibtex
@software{Zeyilin_SwanLab_2023,
  author = {Zeyi Lin, Shaohong Chen, Kang Li, Qiushan Jiang, Zirui Cai,  Kaifang Ji and {The SwanLab team}},
  license = {Apache-2.0},
  title = {{SwanLab}},
  url = {https://github.com/aimhubio/aim},
  year = {2023}
}
```



## 为SwanLab做出贡献

考虑为 SwanLab 做出贡献吗？首先，请花点时间阅读 CONTRIBUTING.md 指南。

同时，我们非常欢迎通过社交媒体、活动和会议的分享来支持SwanLab。

**贡献者们**

<a href="https://github.com/swanhubx/swanlab/graphs/contributors">
  <img src="https://contrib.rocks/image?repo=swanhubx/swanlab" />
</a>

(上图由 [contrib.rocks](https://contrib.rocks) 制作)

<br>



## 协议

本仓库遵循 [Apache 2.0 License](https://github.com/SwanHubX/SwanLab/blob/main/LICENSE) 开源协议

<div align="center">

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="readme_files/swanlab-logo-type2-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="readme_files/swanlab-logo-type2-light.svg">
  <img alt="SwanLab" src="readme_files/swanlab-logo-type2-light.svg" width="300" height="130">
</picture>

一个专业、现代化设计的AI训练分析平台  
面向模型训练团队，与50+主流框架集成，与你的实验代码轻松结合

<a href="https://swanlab.cn">🔥SwanLab 在线版</a> · <a href="https://docs.swanlab.cn">📃 文档</a> · <a href="https://github.com/swanhubx/swanlab/issues">报告问题</a> · <a href="https://rcnpx636fedp.feishu.cn/share/base/form/shrcnwT9dTedqUchdpgzPzKvcmd">建议反馈</a> · <a href="https://docs.swanlab.cn/guide_cloud/general/changelog.html">更新日志</a> · <img height="16" width="16" src="https://raw.githubusercontent.com/SwanHubX/assets/main/community.svg" alt="swanlab community Logo" /> <a href="https://swanlab.cn/benchmarks">基线社区</a> · <a href="https://github.com/SwanHubX/SwanLab-Skill">SKILL</a>

[![][release-shield]][release-link]
[![][dockerhub-shield]][dockerhub-link]
[![][github-stars-shield]][github-stars-link]
[![][github-issues-shield]][github-issues-shield-link]
[![][github-contributors-shield]][github-contributors-link]
[![][license-shield]][license-shield-link]  
[![][tracking-swanlab-shield]][tracking-swanlab-shield-link]
[![][last-commit-shield]][last-commit-shield-link]
[![][pypi-version-shield]][pypi-version-shield-link]
[![][wechat-shield]][wechat-shield-link]
[![][pypi-downloads-shield]][pypi-downloads-shield-link]
[![][colab-shield]][colab-shield-link]

![](readme_files/swanlab-overview.png)

中文 / [English](README_EN.md) / [日本語](README_JP.md) / [Русский](README_RU.md)

👋 加入我们的[微信群](https://docs.swanlab.cn/guide_cloud/community/online-support.html)

<a href="https://hellogithub.com/repository/b442a9fa270e4ccb8847c9ee3445e41b" target="_blank"><img src="https://abroad.hellogithub.com/v1/widgets/recommend.svg?rid=b442a9fa270e4ccb8847c9ee3445e41b&claim_uid=Oh5UaGjfrblg0yZ" alt="Featured｜HelloGitHub" style="width: 250px; height: 54px;" width="250" height="54" /></a>

</div>

<br/>

## 目录

- [👋🏻 什么是SwanLab](#-什么是swanlab)
- [📃 在线演示](#-在线演示)
- [🏁 快速开始](#-快速开始)
- [💻 自托管](#-自托管)
- [🔥 实战案例](#-实战案例)
- [🎮 硬件记录](#-硬件记录)
- [🚗 框架集成](#-框架集成)
- [🔌 插件与API](#-插件与api)
- [🆚 与熟悉的工具的比较](#-与熟悉的工具的比较)
- [👥 社区](#-社区)
- [📃 协议](#-协议)

<br/>

## 👋🏻 什么是SwanLab

SwanLab 是一款**AI训练分析与指标观测平台**，面向模型训练团队，提供训练可视化、自动日志记录、超参数记录、实验对比、多人协同等功能，帮助团队快速发现训练问题，加速模型迭代。

在SwanLab上，研究者能基于直观的可视化图表发现训练问题，对比多个实验找到研究灵感，并通过**在线网页**的分享与基于组织的**多人协同训练**，打破团队沟通的壁垒，提高组织训练效率。

https://github.com/user-attachments/assets/7965fec4-c8b0-4956-803d-dbf177b44f54

以下是其核心特性列表：

**1. 📊 实验指标与超参数跟踪**: 极简的代码嵌入您的机器学习 pipeline，跟踪记录训练关键指标

- ☁️ 支持**云端**使用（类似Weights & Biases），随时随地查看训练进展。[手机看实验的方法](https://docs.swanlab.cn/guide_cloud/general/app.html)

- 📝 支持**超参数记录**、**指标总结**、**表格分析**

- 🌸 **可视化训练过程**: 通过UI界面对实验跟踪数据进行可视化，可以让训练师直观地看到实验每一步的结果，分析指标走势，判断哪些变化导致了模型效果的提升，从而整体性地提升模型迭代效率。

- **支持的元数据类型**：标量指标、图像、音频、文本、视频、3D点云、生物化学分子、Echarts自定义图表...

![swanlab-table](readme_files/molecule.gif)

- **支持的图表类型**：折线图、媒体图（图像、音频、文本、视频）、3D点云、生物化学分子、柱状图、散点图、箱线图、热力图、饼状图、雷达图、[自定义图表](https://docs.swanlab.cn/guide_cloud/experiment_track/log-custom-chart.html)...

[![swanlab-echarts](readme_files/echarts.png)](https://swanlab.cn/@ZeyiLin/swanlab-echarts-demo/charts)

- **LLM生成内容可视化组件**：为大语言模型训练场景打造的文本内容可视化图表，支持Markdown渲染

![text-chart](readme_files/text-chart.gif)

- **后台自动记录**：日志logging、硬件环境、Git 仓库、Python 环境、Python 库列表、项目运行目录

- **断点续训记录**：支持在训练完成/中断后，补充新的指标数据到同个实验中

**2. ⚡️ 全面的框架集成**: PyTorch、🤗HuggingFace Transformers、PyTorch Lightning、🦙LLaMA Factory、MMDetection、Ultralytics、PaddleDetection、LightGBM、XGBoost、Keras、TensorBoard、Weights&Biases、Swift、XTuner、Stable-Baselines3、Hydra 在内的 **30+** 框架

![](readme_files/integrations.png)

**3. 💻 硬件监控**: 支持实时记录与监控CPU、NPU（**昇腾Ascend**）、GPU（**英伟达Nvidia**）、AMD（**AMD ROCm**）、MLU（**寒武纪Cambricon**）、XPU（**昆仑芯Kunlunxin**）、DCU（**海光DCU**）、MetaX GPU（**沐曦XPU**）、Moore Threads GPU（**摩尔线程**）、Iluvatar GPU（**天数智芯**）、内存的系统级硬件指标

**4. 📦 实验管理**: 通过专为训练场景设计的集中式仪表板，通过整体视图速览全局，快速管理多个项目与实验

**5. 🆚 比较结果**: 通过在线表格与对比图表比较不同实验的超参数和结果，挖掘迭代灵感

![](readme_files/swanlab-table.png)

**6. 👥 在线协作**: 您可以与团队进行协作式训练，支持将实验实时同步在一个项目下，您可以在线查看团队的训练记录，基于结果发表看法与建议

**7. ✉️ 分享结果**: 复制和发送持久的 URL 来共享每个实验，方便地发送给伙伴，或嵌入到在线笔记中

**8. 💻 支持自托管**: 支持离线环境使用，自托管的社区版同样可以查看仪表盘与管理实验，[使用攻略](#-自托管)

**9. 🔌 插件拓展**: 支持通过插件拓展SwanLab的使用场景，比如 [飞书通知](https://docs.swanlab.cn/plugin/notification-lark.html)、[Slack通知](https://docs.swanlab.cn/plugin/notification-slack.html)、[CSV记录器](https://docs.swanlab.cn/plugin/writer-csv.html)等

> \[!IMPORTANT]
>
> **收藏项目**，你将从 GitHub 上无延迟地接收所有发布通知～ ⭐️

![star-us](readme_files/star-us.png)

<br>

## 📃 在线演示

来看看 SwanLab 的在线演示：

|              [ResNet50 猫狗分类][demo-cats-dogs]               |                  [Yolov8-COCO128 目标检测][demo-yolo]                   |
| :------------------------------------------------------------: | :---------------------------------------------------------------------: |
|          [![][demo-cats-dogs-image]][demo-cats-dogs]           |                    [![][demo-yolo-image]][demo-yolo]                    |
| 跟踪一个简单的 ResNet50 模型在猫狗数据集上训练的图像分类任务。 | 使用 Yolov8 在 COCO128 数据集上进行目标检测任务，跟踪训练超参数和指标。 |

|             [Qwen2 指令微调][demo-qwen2-sft]              |               [LSTM Google 股票预测][demo-google-stock]                |
| :-------------------------------------------------------: | :--------------------------------------------------------------------: |
|        [![][demo-qwen2-sft-image]][demo-qwen2-sft]        |           [![][demo-google-stock-image]][demo-google-stock]            |
| 跟踪 Qwen2 大语言模型的指令微调训练，完成简单的指令遵循。 | 使用简单的 LSTM 模型在 Google 股价数据集上训练，实现对未来股价的预测。 |

|         [ResNeXt101 音频分类][demo-audio-classification]          |          [Qwen2-VL COCO数据集微调][demo-qwen2-vl]          |
| :---------------------------------------------------------------: | :--------------------------------------------------------: |
| [![][demo-audio-classification-image]][demo-audio-classification] |         [![][demo-qwen2-vl-image]][demo-qwen2-vl]          |
|         从ResNet到ResNeXt在音频分类任务上的渐进式实验过程         | 基于Qwen2-VL多模态大模型，在COCO2014数据集上进行Lora微调。 |

|  [EasyR1 多模态LLM RL训练][demo-easyr1-rl]  |    [Qwen2.5-0.5B GRPO训练][demo-qwen2-grpo]     |
| :-----------------------------------------: | :---------------------------------------------: |
| [![][demo-easyr1-rl-image]][demo-easyr1-rl] |  [![][demo-qwen2-grpo-image]][demo-qwen2-grpo]  |
|     使用EasyR1框架进行多模态LLM RL训练      | 基于Qwen2.5-0.5B模型在GSM8k数据集上进行GRPO训练 |

[更多案例](https://docs.swanlab.cn/examples/mnist.html)

<br>

## 🏁 快速开始

### 1.安装

```bash
pip install swanlab
```

<details><summary>源码安装</summary>

如果你想体验最新的特性，可以使用源码安装。

```bash
# 方式一
git clone https://github.com/SwanHubX/SwanLab.git
pip install -e .

# 方式二
pip install git+https://github.com/SwanHubX/SwanLab.git
```

</details>

<details><summary>离线看板拓展安装</summary>

[离线看板文档](https://docs.swanlab.cn/guide_cloud/offline-board.html)

```bash
pip install 'swanlab[dashboard]'
```

</details>

### 2.登录并获取 API Key

1. 免费[注册账号](https://swanlab.cn)

2. 登录账号，在用户设置 > [API Key](https://swanlab.cn/settings) 里复制您的 API Key

3. 打开终端，输入：

```bash
swanlab login
```

出现提示时，输入您的 API Key，按下回车，完成登陆。

### 3.将 SwanLab 与你的代码集成

```python
import swanlab

# 初始化一个新的swanlab实验
swanlab.init(
    project="my-first-ml",
    config={"learning-rate": 0.003},
)

# 记录指标
for i in range(10):
    swanlab.log({"loss": i, "acc": i})
```

大功告成！前往[SwanLab](https://swanlab.cn)查看你的第一个 SwanLab 实验。

<br>

## 💻 自托管

自托管社区版支持离线查看 SwanLab 仪表盘。

![swanlab-kubernetes](./readme_files/swanlab-kubernetes.png)

详细部署文档见：

- [🔗 Kubernetes部署文档](https://docs.swanlab.cn/self_host/kubernetes/deploy.html)
- [🔗 Docker部署文档](https://docs.swanlab.cn/self_host/docker/deploy.html)
- [🔗 从Docker到Kubernetes迁移文档](https://docs.swanlab.cn/self_host/docker/migration-docker-kubernetes.html)

<br>

## 🔥 实战案例

**使用SwanLab的优秀教程开源项目：**

- [happy-llm](https://github.com/datawhalechina/happy-llm)：从零开始的大语言模型原理与实践教程 ![GitHub Repo stars](https://img.shields.io/github/stars/datawhalechina/happy-llm)
- [self-llm](https://github.com/datawhalechina/self-llm)：《开源大模型食用指南》针对中国宝宝量身打造的基于Linux环境快速微调（全参数/Lora）、部署国内外开源大模型（LLM）/多模态大模型（MLLM）教程 ![GitHub Repo stars](https://img.shields.io/github/stars/datawhalechina/self-llm)
- [Minimind](https://github.com/jingyaogong/minimind)：🚀🚀 「大模型」2小时完全从0训练26M的小参数GPT！![GitHub Repo stars](https://img.shields.io/github/stars/jingyaogong/minimind)
- [unlock-deepseek](https://github.com/datawhalechina/unlock-deepseek)：DeepSeek 系列工作解读、扩展和复现 ![GitHub Repo stars](https://img.shields.io/github/stars/datawhalechina/unlock-deepseek)
- [Qwen3-SmVL](https://github.com/ShaohonChen/Qwen3-SmVL): 将SmolVLM2的视觉头与Qwen3-0.6B模型进行了拼接微调 ![GitHub Repo stars](https://img.shields.io/github/stars/ShaohonChen/Qwen3-SmVL)
- [OPPO/Agent_Foundation_Models](https://github.com/OPPO-PersonalAI/Agent_Foundation_Models): 通过多Agent蒸馏和Agent RL的端到端Agent基础模型。 ![GitHub Repo stars](https://img.shields.io/github/stars/OPPO-PersonalAI/Agent_Foundation_Models)
- [Tree-GRPO](https://github.com/AMAP-ML/Tree-GRPO): [ICLR 2026] 树搜索在LLM Agent RL中的应用 ![GitHub Repo stars](https://img.shields.io/github/stars/AMAP-ML/Tree-GRPO)
- [FIPO](https://github.com/qwenpilot/FIPO): This code implements the algorithm of FIPO, a value-free RL recipe for eliciting deeper reasoning from a clean base model. ![GitHub Repo stars](https://img.shields.io/github/stars/qwenpilot/FIPO)
- [llm-agent-rl-lab](https://github.com/KMnO4-zx/llm-agent-rl-lab): Reproducing and studying RL algorithms for LLM agents, including PPO, GRPO, GSPO, DAPO, OPD and beyond. ![GitHub Repo stars](https://img.shields.io/github/stars/KMnO4-zx/llm-agent-rl-lab)

**使用SwanLab的优秀论文：**

- [Accelerating AlphaZero Training for Connect6 via Behavioral Cloning Initialization and Prior Action Masking in Resource-Constrained Environments](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=6730709)
- [SeePhys Pro: Diagnosing Modality Transfer and Blind-Training Effects in Multimodal RLVR for Physics Reasoning](https://arxiv.org/abs/2605.09266)
- [MemFactory: Unified Inference & Training Framework for Agent Memory](https://arxiv.org/abs/2603.29493)
- [ILL: A Lightweight Large Language Model for Legal and Courtroom Assistance](https://www.mdpi.com/2079-9292/15/11/2401)
- [Retrieval, Reward, and Training Protocols: What Matters in Training Search Agents?](https://arxiv.org/pdf/2605.27881)
- [MolAct: An Agentic RL Framework for Molecular Editing and Property Optimization](https://arxiv.org/abs/2512.20135)
- [CQLLM: A Framework for Generating CodeQL Security Vulnerability Detection Code Based on Large Language Model](https://www.preprints.org/manuscript/202510.1458)
- [Animation Needs Attention: A Holistic Approach to Slides Animation Comprehension with Visual-Language Models](https://arxiv.org/abs/2507.03916)
- [Efficient Model Fine-Tuning with LoRA for Biomedical Named Entity Recognition](https://ieeexplore.ieee.org/abstract/document/11082049/)
- [SpectrumWorld: Artificial Intelligence Foundation for Spectroscopy](https://arxiv.org/abs/2508.01188)
- [CodeBoost: Boosting Code LLMs by Squeezing Knowledge from Code Snippets with RL](https://arxiv.org/pdf/2508.05242)
- [A Joint Classification Method for Traditional Chinese Medicine Diseases and Syndromes Based on BertChinese-RCNNATTN](https://dl.acm.org/doi/full/10.1145/3759972.3759979)
- [A Gradient-Norm-Aware Optimizer for Symmetry-Preserving and Stable Deep Learning](https://www.mdpi.com/2073-8994/18/2/374)
- [TrafficBM: A Dual-Modality Pre-Training Framework for Network Traffic Classification](https://ieeexplore.ieee.org/document/11304658)
- [Using Human Cumulative Prospect Theory to Understand Large Language Models Decision-Making](https://ieeexplore.ieee.org/document/11300712)
- [BiPACE: Bisimulation-Guided Policy Optimization with Action Counterfactual Estimation for LLM Agents](https://arxiv.org/abs/2606.25556)

**教程文章：**

- [MNIST手写体识别](https://docs.swanlab.cn/examples/mnist.html)
- [FashionMNIST服装分类](https://docs.swanlab.cn/examples/fashionmnist.html)
- [Cifar10图像分类](https://docs.swanlab.cn/examples/cifar10.html)
- [Resnet猫狗分类](https://docs.swanlab.cn/examples/cats_dogs_classification.html)
- [Yolo目标检测](https://docs.swanlab.cn/examples/yolo.html)
- [UNet医学影像分割](https://docs.swanlab.cn/examples/unet-medical-segmentation.html)
- [音频分类](https://docs.swanlab.cn/examples/audio_classification.html)
- [DQN强化学习-推车倒立摆](https://docs.swanlab.cn/examples/dqn_cartpole.html)
- [LSTM Google股票预测](https://docs.swanlab.cn/examples/audio_classification.html)
- [BERT文本分类](https://docs.swanlab.cn/examples/bert.html)
- [Stable Diffusion文生图微调](https://docs.swanlab.cn/examples/stable_diffusion.html)
- [LLM预训练](https://docs.swanlab.cn/examples/pretrain_llm.html)
- [GLM4指令微调](https://docs.swanlab.cn/examples/glm4-instruct.html)
- [Qwen下游任务训练](https://docs.swanlab.cn/examples/qwen_finetune.html)
- [NER命名实体识别](https://docs.swanlab.cn/examples/ner.html)
- [Qwen3医学模型微调](https://docs.swanlab.cn/examples/qwen3-medical.html)
- [Qwen2-VL多模态大模型微调实战](https://docs.swanlab.cn/examples/qwen_vl_coco.html)
- [GRPO大模型强化学习](https://docs.swanlab.cn/examples/qwen_grpo.html)
- [Qwen3-SmVL-0.6B多模态模型训练](https://docs.swanlab.cn/examples/qwen3_smolvlm_muxi.html)
- [LeRobot 具身智能入门](https://docs.swanlab.cn/examples/robot/lerobot-guide.html)
- [GLM-4.5-Air-LoRA 及 SwanLab 可视化记录](https://github.com/datawhalechina/self-llm/blob/master/models/GLM-4.5-Air/03-GLM-4.5-Air-Lora%20%E5%8F%8A%20Swanlab%20%E5%8F%AF%E8%A7%86%E5%8C%96%E5%BE%AE%E8%B0%83.md)
- [RAG怎么做？SwanLab文档助手方案开源了](https://docs.swanlab.cn/course/prompt_engineering_course/11-swanlab_rag/1.swanlab-rag.html)

🌟如果你有想收录的教程，欢迎提交PR！

<br>

## 🎮 硬件记录

SwanLab会对AI训练过程中所使用的**硬件信息**和**资源使用情况**进行记录，下面是支持情况表格：

| 硬件        | 信息记录 | 资源监控 | 脚本                                                                                                                                           |
| ----------- | -------- | -------- | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| 英伟达GPU   | ✅       | ✅       | [nvidia.py](https://github.com/SwanHubX/SwanLab/blob/main/swanlab/sdk/internal/probe_python/hardware_vendor/accelerator/nvidia.py)             |
| AMD GPU     | ✅       | ✅       | [amd.py](https://github.com/SwanHubX/SwanLab/blob/main/swanlab/sdk/internal/probe_python/hardware_vendor/accelerator/amd.py)                   |
| 昇腾NPU     | ✅       | ✅       | [ascend.py](https://github.com/SwanHubX/SwanLab/blob/main/swanlab/sdk/internal/probe_python/hardware_vendor/accelerator/huawei.py)             |
| 寒武纪MLU   | ✅       | ✅       | [cambricon.py](https://github.com/SwanHubX/SwanLab/blob/main/swanlab/sdk/internal/probe_python/hardware_vendor/accelerator/cambricon.py)       |
| 昆仑芯XPU   | ✅       | ✅       | [kunlunxin.py](https://github.com/SwanHubX/SwanLab/blob/main/swanlab/sdk/internal/probe_python/hardware_vendor/accelerator/kunlunxin.py)       |
| 摩尔线程GPU | ✅       | ✅       | [moorethreads.py](https://github.com/SwanHubX/SwanLab/blob/main/swanlab/sdk/internal/probe_python/hardware_vendor/accelerator/moorethreads.py) |
| 沐曦GPU     | ✅       | ✅       | [metax.py](https://github.com/SwanHubX/SwanLab/blob/main/swanlab/sdk/internal/probe_python/hardware_vendor/accelerator/metax.py)               |
| 天数智芯GPU | ✅       | ✅       | [iluvatar.py](https://github.com/SwanHubX/SwanLab/blob/main/swanlab/sdk/internal/probe_python/hardware_vendor/accelerator/iluvatar.py)         |
| 海光DCU     | ✅       | ✅       | [hygon.py](https://github.com/SwanHubX/SwanLab/blob/main/swanlab/sdk/internal/probe_python/hardware_vendor/accelerator/hygon.py)               |
| CPU         | ✅       | ✅       | [cpu.py](https://github.com/SwanHubX/SwanLab/blob/main/swanlab/sdk/internal/probe_python/hardware_vendor/cpu.py)                               |
| 内存        | ✅       | ✅       | [memory.py](https://github.com/SwanHubX/SwanLab/blob/main/swanlab/sdk/internal/probe_python/hardware_vendor/memory.py)                         |

如果你希望记录其他硬件，欢迎提交Issue与PR！

<br>

## 🚗 框架集成

将你最喜欢的框架与 SwanLab 结合使用！  
下面是我们已集成的框架列表，欢迎提交 [Issue](https://github.com/swanhubx/swanlab/issues) 来反馈你想要集成的框架。

**基础框架**

- [PyTorch](https://docs.swanlab.cn/guide_cloud/integration/integration-pytorch.html)
- [MindSpore](https://docs.swanlab.cn/guide_cloud/integration/integration-ascend.html)
- [Keras](https://docs.swanlab.cn/guide_cloud/integration/integration-keras.html)

**LLM训练框架**

- [HuggingFace Transformers](https://docs.swanlab.cn/guide_cloud/integration/integration-huggingface-transformers.html)
- [LLaMA Factory](https://docs.swanlab.cn/guide_cloud/integration/integration-llama-factory.html)
- [MS-Swift](https://docs.swanlab.cn/guide_cloud/integration/integration-swift.html)
- [Unsloth](https://docs.swanlab.cn/guide_cloud/integration/integration-unsloth.html)
- [MLX-LM](https://docs.swanlab.cn/guide_cloud/integration/integration-mlx-lm.html)
- [Torchtune](https://docs.swanlab.cn/guide_cloud/integration/integration-pytorch-torchtune.html)
- [PaddleNLP](https://docs.swanlab.cn/guide_cloud/integration/integration-paddlenlp.html)
- [Sentence Transformers](https://docs.swanlab.cn/guide_cloud/integration/integration-sentence-transformers.html)
- [XTuner](https://docs.swanlab.cn/guide_cloud/integration/integration-xtuner.html)
- [OpenMind](https://modelers.cn/docs/zh/openmind-library/1.0.0/basic_tutorial/finetune/finetune_pt.html#%E8%AE%AD%E7%BB%83%E7%9B%91%E6%8E%A7)

**LLM强化学习框架**

- [veRL](https://docs.swanlab.cn/guide_cloud/integration/integration-verl.html)
- [HuggingFace trl](https://docs.swanlab.cn/guide_cloud/integration/integration-huggingface-trl.html)
- [NVIDIA-NeMo RL](https://docs.swanlab.cn/guide_cloud/integration/integration-nvidia-nemo-rl.html)
- [EasyR1](https://docs.swanlab.cn/guide_cloud/integration/integration-easyr1.html)
- [AReaL](https://docs.swanlab.cn/guide_cloud/integration/integration-areal.html)
- [ROLL](https://docs.swanlab.cn/guide_cloud/integration/integration-roll.html)

**机器人框架**

- [RLinf](https://docs.swanlab.cn/guide_cloud/integration/integration-rlinf.html)

**文生图/视频训练框架**

- [DiffSynth Studio](https://docs.swanlab.cn/guide_cloud/integration/integration-diffsynth-studio.html)

**深度学习框架**

- [PyTorch Lightning](https://docs.swanlab.cn/guide_cloud/integration/integration-pytorch-lightning.html)
- [MMEngine](https://docs.swanlab.cn/guide_cloud/integration/integration-mmengine.html)
- [FastAI](https://docs.swanlab.cn/guide_cloud/integration/integration-fastai.html)

**计算机视觉**

- [Ultralytics](https://docs.swanlab.cn/guide_cloud/integration/integration-ultralytics.html)
- [MMDetection](https://docs.swanlab.cn/guide_cloud/integration/integration-mmdetection.html)
- [MMSegmentation](https://docs.swanlab.cn/guide_cloud/integration/integration-mmsegmentation.html)
- [PaddleDetection](https://docs.swanlab.cn/guide_cloud/integration/integration-paddledetection.html)
- [PaddleYOLO](https://docs.swanlab.cn/guide_cloud/integration/integration-paddleyolo.html)

**机器学习框架**

- [LightGBM](https://docs.swanlab.cn/guide_cloud/integration/integration-lightgbm.html)
- [XGBoost](https://docs.swanlab.cn/guide_cloud/integration/integration-xgboost.html)
- [CatBoost](https://docs.swanlab.cn/guide_cloud/integration/integration-catboost.html)

**评估框架**

- [EvalScope](https://docs.swanlab.cn/guide_cloud/integration/integration-evalscope.html)

**传统强化学习框架**

- [Stable Baseline3](https://docs.swanlab.cn/guide_cloud/integration/integration-sb3.html)

**其他框架：**

- [Tensorboard](https://docs.swanlab.cn/guide_cloud/integration/integration-tensorboard.html)
- [Weights&Biases](https://docs.swanlab.cn/guide_cloud/integration/integration-wandb.html)
- [MLFlow](https://docs.swanlab.cn/guide_cloud/integration/integration-mlflow.html)
- [HuggingFace Accelerate](https://docs.swanlab.cn/guide_cloud/integration/integration-huggingface-accelerate.html)
- [Ray](https://docs.swanlab.cn/guide_cloud/integration/integration-ray.html)
- [Hydra](https://docs.swanlab.cn/guide_cloud/integration/integration-hydra.html)
- [Omegaconf](https://docs.swanlab.cn/guide_cloud/integration/integration-omegaconf.html)
- [SpecForge](https://docs.swanlab.cn/guide_cloud/integration/integration-specforge.html)

[更多集成](https://docs.swanlab.cn/guide_cloud/integration/integration-any-library.html)

<br>

## 🔌 插件与API

欢迎通过插件来拓展SwanLab的功能，增强你的实验管理体验！

- [自定义你的插件](https://docs.swanlab.cn/plugin/custom-plugin.html)
- [邮件通知](https://docs.swanlab.cn/plugin/notification-email.html)
- [飞书通知](https://docs.swanlab.cn/plugin/notification-lark.html)
- [钉钉通知](https://docs.swanlab.cn/plugin/notification-dingtalk.html)
- [企业微信通知](https://docs.swanlab.cn/plugin/notification-wxwork.html)
- [Discord通知](https://docs.swanlab.cn/plugin/notification-discord.html)
- [Slack通知](https://docs.swanlab.cn/plugin/notification-slack.html)
- [Bark通知](https://docs.swanlab.cn/plugin/notification-bark.html)
- [CSV记录器](https://docs.swanlab.cn/plugin/writer-csv.html)

开放接口：

- [OpenAPI](https://docs.swanlab.cn/api/py-openapi.html)

<br>

## 🆚 与熟悉的工具的比较

### Tensorboard vs SwanLab

- **☁️ 支持在线使用**：
  通过 SwanLab 可以方便地将训练实验在云端在线同步与保存，便于远程查看训练进展、管理历史项目、分享实验链接、发送实时消息通知、多端看实验等。而 Tensorboard 是一个离线的实验跟踪工具。

- **👥 多人协作**：
  在进行多人、跨团队的机器学习协作时，通过 SwanLab 可以轻松管理多人的训练项目、分享实验链接、跨空间交流讨论。而 Tensorboard 主要为个人设计，难以进行多人协作和分享实验。

- **💻 持久、集中的仪表板**：
  无论你在何处训练模型，无论是在本地计算机上、在实验室集群还是在公有云的 GPU 实例中，你的结果都会记录到同一个集中式仪表板中。而使用 TensorBoard 需要花费时间从不同的机器复制和管理
  TFEvent 文件。

- **💪 更强大的表格**：
  通过 SwanLab 表格可以查看、搜索、过滤来自不同实验的结果，可以轻松查看数千个模型版本并找到适合不同任务的最佳性能模型。
  TensorBoard 不适用于大型项目。

### Weights and Biases vs SwanLab

- Weights and Biases 是一个必须联网使用的闭源 MLOps 平台

- SwanLab 不仅支持联网使用，也支持开源、免费、自托管的版本

<br>

## 👥 社区

### 周边仓库

- [self-hosted](https://github.com/swanhubx/self-hosted)：私有化部署脚本仓库
- [SwanLab-Docs](https://github.com/swanhubx/swanlab-docs)：官方文档仓库
- [SwanLab-Dashboard](https://github.com/swanhubx/swanlab-dashboard)：离线看板仓库，存放了由`swanlab watch`打开的轻量离线看板的web代码

### 社区与支持

- [GitHub Issues](https://github.com/SwanHubX/SwanLab/issues)：使用 SwanLab 时遇到的错误和问题
- [电子邮件支持](mailto:contact@swanlab.cn)：反馈关于使用 SwanLab 的问题
- <a href="https://docs.swanlab.cn/guide_cloud/community/online-support.html">微信交流群</a>：交流使用 SwanLab 的问题、分享最新的 AI 技术

### SwanLab README 徽章

如果你喜欢在工作中使用 SwanLab，请将 SwanLab 徽章添加到你的 README 中：

[![][tracking-swanlab-shield]][tracking-swanlab-shield-link]、[![][visualize-swanlab-shield]][visualize-swanlab-shield-link]

```
[![](https://raw.githubusercontent.com/SwanHubX/assets/main/badge2.svg)](your experiment url)
[![](https://raw.githubusercontent.com/SwanHubX/assets/main/badge1.svg)](your experiment url)
```

更多设计素材：[assets](https://github.com/SwanHubX/assets)

### 在论文中引用 SwanLab

如果您发现 SwanLab 对您的研究之旅有帮助，请考虑以下列格式引用：

```bibtex
@software{Zeyilin_SwanLab_2023,
  author = {Zeyi Lin, Shaohong Chen, Kang Li, Qiushan Jiang, Zirui Cai,  Kaifang Ji and {The SwanLab team}},
  doi = {10.5281/zenodo.11100550},
  license = {Apache-2.0},
  title = {{SwanLab}},
  url = {https://github.com/swanhubx/swanlab},
  year = {2023}
}
```

### 为 SwanLab 做出贡献

考虑为 SwanLab 做出贡献吗？首先，请花点时间阅读 [贡献指南](CONTRIBUTING.md)。

同时，我们非常欢迎通过社交媒体、活动和会议的分享来支持 SwanLab，衷心感谢！

<br>

**Contributors**

<a href="https://github.com/swanhubx/swanlab/graphs/contributors">
  <img src="https://contrib.rocks/image?repo=swanhubx/swanlab" />
</a>

<br>

<img src="./readme_files/swanlab-and-user.png" width="50%" />

## 📃 协议

本仓库遵循 [Apache 2.0 License](https://github.com/SwanHubX/SwanLab/blob/main/LICENSE) 开源协议

## Star History

<a href="https://www.star-history.com/?repos=swanhubx%2Fswanlab&type=date&legend=top-left">
 <picture>
   <source media="(prefers-color-scheme: dark)" srcset="https://api.star-history.com/chart?repos=swanhubx/swanlab&type=date&theme=dark&legend=top-left&sealed_token=77bt2XJC-p2PM1tRv6wxVLEw_RaNVx2MVWL-FFg2yjKRstx-vq6M_Y0yszVSAAYT4CvkjYMlgSOPU-mvyq_FfRWs7hh2-VFJiwQMtvyr1-7i0F6VgeQwQw" />
   <source media="(prefers-color-scheme: light)" srcset="https://api.star-history.com/chart?repos=swanhubx/swanlab&type=date&legend=top-left&sealed_token=77bt2XJC-p2PM1tRv6wxVLEw_RaNVx2MVWL-FFg2yjKRstx-vq6M_Y0yszVSAAYT4CvkjYMlgSOPU-mvyq_FfRWs7hh2-VFJiwQMtvyr1-7i0F6VgeQwQw" />
   <img alt="Star History Chart" src="https://api.star-history.com/chart?repos=swanhubx/swanlab&type=date&legend=top-left&sealed_token=77bt2XJC-p2PM1tRv6wxVLEw_RaNVx2MVWL-FFg2yjKRstx-vq6M_Y0yszVSAAYT4CvkjYMlgSOPU-mvyq_FfRWs7hh2-VFJiwQMtvyr1-7i0F6VgeQwQw" />
 </picture>
</a>

<!-- link -->

[release-shield]: https://img.shields.io/github/v/release/swanhubx/swanlab?color=369eff&labelColor=black&logo=github&style=flat-square
[release-link]: https://github.com/swanhubx/swanlab/releases
[license-shield]: https://img.shields.io/badge/license-apache%202.0-white?labelColor=black&style=flat-square
[license-shield-link]: https://github.com/SwanHubX/SwanLab/blob/main/LICENSE
[last-commit-shield]: https://img.shields.io/github/last-commit/swanhubx/swanlab?color=c4f042&labelColor=black&style=flat-square
[last-commit-shield-link]: https://github.com/swanhubx/swanlab/commits/main
[pypi-version-shield]: https://img.shields.io/pypi/v/swanlab?color=orange&labelColor=black&style=flat-square
[pypi-version-shield-link]: https://pypi.org/project/swanlab/
[pypi-downloads-shield]: https://static.pepy.tech/badge/swanlab?labelColor=black&style=flat-square
[pypi-downloads-shield-link]: https://pepy.tech/project/swanlab
[swanlab-cloud-shield]: https://img.shields.io/badge/Product-SwanLab云端版-636a3f?labelColor=black&style=flat-square
[swanlab-cloud-shield-link]: https://swanlab.cn/
[wechat-shield]: https://img.shields.io/badge/WeChat-微信-4cb55e?labelColor=black&style=flat-square
[wechat-shield-link]: https://docs.swanlab.cn/guide_cloud/community/online-support.html
[colab-shield]: https://colab.research.google.com/assets/colab-badge.svg
[colab-shield-link]: https://colab.research.google.com/drive/1RWsrY_1bS8ECzaHvYtLb_1eBkkdzekR3?usp=sharing
[github-stars-shield]: https://img.shields.io/github/stars/swanhubx/swanlab?labelColor&style=flat-square&color=ffcb47
[github-stars-link]: https://github.com/swanhubx/swanlab
[github-issues-shield]: https://img.shields.io/github/issues/swanhubx/swanlab?labelColor=black&style=flat-square&color=ff80eb
[github-issues-shield-link]: https://github.com/swanhubx/swanlab/issues
[github-contributors-shield]: https://img.shields.io/github/contributors/swanhubx/swanlab?color=c4f042&labelColor=black&style=flat-square
[github-contributors-link]: https://github.com/swanhubx/swanlab/graphs/contributors
[demo-cats-dogs]: https://swanlab.cn/@ZeyiLin/Cats_Dogs_Classification/runs/jzo93k112f15pmx14vtxf/chart
[demo-cats-dogs-image]: readme_files/example-catsdogs.png
[demo-yolo]: https://swanlab.cn/@ZeyiLin/ultratest/runs/yux7vclmsmmsar9ear7u5/chart
[demo-yolo-image]: readme_files/example-yolo.png
[demo-qwen2-sft]: https://swanlab.cn/@ZeyiLin/Qwen2-fintune/runs/cfg5f8dzkp6vouxzaxlx6/chart
[demo-qwen2-sft-image]: readme_files/example-qwen2.png
[demo-google-stock]: https://swanlab.cn/@ZeyiLin/Google-Stock-Prediction/charts
[demo-google-stock-image]: readme_files/example-lstm.png
[demo-audio-classification]: https://swanlab.cn/@ZeyiLin/PyTorch_Audio_Classification/charts
[demo-audio-classification-image]: readme_files/example-audio-classification.png
[demo-qwen2-vl]: https://swanlab.cn/@ZeyiLin/Qwen2-VL-finetune/runs/pkgest5xhdn3ukpdy6kv5/chart
[demo-qwen2-vl-image]: readme_files/example-qwen2-vl.jpg
[demo-easyr1-rl]: https://swanlab.cn/@Kedreamix/easy_r1/runs/wzezd8q36bb6dlza6wtpc/chart
[demo-easyr1-rl-image]: readme_files/example-easyr1-rl.png
[demo-qwen2-grpo]: https://swanlab.cn/@kmno4/Qwen-R1/runs/t0zr3ak5r7188mjbjgdsc/chart
[demo-qwen2-grpo-image]: readme_files/example-qwen2-grpo.png
[tracking-swanlab-shield-link]: https://swanlab.cn
[tracking-swanlab-shield]: https://raw.githubusercontent.com/SwanHubX/assets/main/badge2.svg
[visualize-swanlab-shield-link]: https://swanlab.cn
[visualize-swanlab-shield]: https://raw.githubusercontent.com/SwanHubX/assets/main/badge1.svg
[dockerhub-shield]: https://img.shields.io/docker/v/swanlab/swanlab-next?color=369eff&label=docker&labelColor=black&logoColor=white&style=flat-square
[dockerhub-link]: https://hub.docker.com/r/swanlab/swanlab-next/tags

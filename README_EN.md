<div align="center">

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="readme_files/swanlab-logo-type2-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="readme_files/swanlab-logo-type2-light.svg">
  <img alt="SwanLab" src="readme_files/swanlab-logo-type2-light.svg" width="300" height="130">
</picture>

A Modern, Purpose-Built AI Training Analytics Platform  
Built for ML Teams — Integrated with 50+ Mainstream Frameworks and Designed to Drop Effortlessly into Your Training Code

<a href="https://swanlab.cn">🔥SwanLab Cloud</a> · <a href="https://docs.swanlab.cn/en/">📃 Documentation</a> · <a href="https://github.com/swanhubx/swanlab/issues">Report Issues</a> · <a href="https://rcnpx636fedp.feishu.cn/share/base/form/shrcnwT9dTedqUchdpgzPzKvcmd">Feedback</a> · <a href="https://docs.swanlab.cn/en/guide_cloud/general/changelog.html">Changelog</a> · <img height="16" width="16" src="https://raw.githubusercontent.com/SwanHubX/assets/main/community.svg" alt="swanlab community Logo" /> <a href="https://swanlab.cn/benchmarks">Benchmarks</a> · <a href="https://github.com/SwanHubX/SwanLab-Skill">SKILL</a>

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

[中文](README.md) / English / [日本語](README_JP.md) / [Русский](README_RU.md)

👋 Join our [WeChat Group](https://docs.swanlab.cn/en/guide_cloud/community/online-support.html)

<a href="https://hellogithub.com/repository/b442a9fa270e4ccb8847c9ee3445e41b" target="_blank"><img src="https://abroad.hellogithub.com/v1/widgets/recommend.svg?rid=b442a9fa270e4ccb8847c9ee3445e41b&claim_uid=Oh5UaGjfrblg0yZ" alt="Featured｜HelloGitHub" style="width: 250px; height: 54px;" width="250" height="54" /></a>

</div>

<br/>

## Table of Contents

- [👋🏻 What is SwanLab](#-what-is-swanlab)
- [📃 Online Demos](#-online-demos)
- [🏁 Quick Start](#-quick-start)
- [💻 Self-Hosting](#-self-hosting)
- [🔥 Case Studies](#-case-studies)
- [🎮 Hardware Monitoring](#-hardware-monitoring)
- [🚗 Framework Integrations](#-framework-integrations)
- [🔌 Plugins & API](#-plugins--api)
- [🆚 Comparison with Alternatives](#-comparison-with-alternatives)
- [👥 Community](#-community)
- [📃 License](#-license)

<br/>

## 👋🏻 What is SwanLab

SwanLab is an **AI training analytics and metrics observability platform** built for model training teams. It provides **training visualization, automated logging, hyperparameter tracking, experiment comparison, and multi-user collaboration** to help teams diagnose training anomalies early and accelerate model iteration.

With SwanLab, researchers can spot training bottlenecks through intuitive interactive charts, compare multiple runs side by side to uncover insights, and eliminate team communication silos via **web-based sharing** and **organization-level collaborative workspaces**.

https://github.com/user-attachments/assets/7965fec4-c8b0-4956-803d-dbf177b44f54

Key features include:

**1. 📊 Experiment Metrics & Hyperparameter Tracking**: Instrument your ML pipeline with minimal code to track and log critical training metrics.

- ☁️ **Cloud-Ready** (similar to Weights & Biases): Monitor training progress anytime, anywhere. [View experiments on mobile](https://docs.swanlab.cn/en/guide_cloud/general/app.html).

- 📝 Supports **hyperparameter logging**, **metric summaries**, and **tabular analysis**.

- 🌸 **Training Visualization**: Visualize experiment tracking data through an interactive UI to inspect step-by-step dynamics, analyze metric trajectories, and pinpoint which architectural or hyperparameter changes drive model improvements.

- **Supported Metadata Types**: Scalar metrics, images, audio, text, video, 3D point clouds, molecular structures, custom ECharts...

![swanlab-table](readme_files/molecule.gif)

- **Supported Chart Types**: Line charts, media panels (images, audio, text, video), 3D point clouds, molecular viewers, bar charts, scatter plots, box plots, heatmaps, pie charts, radar charts, [custom charts](https://docs.swanlab.cn/en/guide_cloud/experiment_track/log-custom-chart.html)...

[![swanlab-echarts](readme_files/echarts.png)](https://swanlab.cn/@ZeyiLin/swanlab-echarts-demo/charts)

- **LLM Generation Visualizer**: Purpose-built text visualization component for LLM training workflows with native Markdown rendering.

![text-chart](readme_files/text-chart.gif)

- **Automated Background Logging**: Console logs (`stdout`/`stderr`), hardware specs, Git repository state, Python runtime environment, installed packages, and working directory.

- **Checkpoint & Run Resumption**: Seamlessly append new metrics to an existing run after training interruption or completion.

**2. ⚡️ Comprehensive Framework Integrations**: Out-of-the-box integration with **30+** frameworks including PyTorch, 🤗HuggingFace Transformers, PyTorch Lightning, 🦙LLaMA Factory, MMDetection, Ultralytics, PaddleDetection, LightGBM, XGBoost, Keras, TensorBoard, Weights & Biases, MS-Swift, XTuner, Stable-Baselines3, and Hydra.

![](readme_files/integrations.png)

**3. 💻 System & Accelerator Monitoring**: Real-time tracking of system-level hardware metrics across CPU, NPU (**Huawei Ascend**), GPU (**NVIDIA**), AMD (**AMD ROCm**), MLU (**Cambricon**), XPU (**Kunlunxin**), DCU (**Hygon**), MetaX GPU (**MetaX**), Moore Threads GPU (**Moore Threads**), Iluvatar GPU (**Iluvatar CoreX**), and system memory.

**4. 📦 Centralized Experiment Management**: Gain a bird's-eye view of all projects and runs through a centralized dashboard tailored for ML workloads.

**5. 🆚 Side-by-Side Run Comparison**: Compare hyperparameters and metrics across runs using interactive tables and overlay plots to guide model iteration.

![](readme_files/swanlab-table.png)

**6. 👥 Team Collaboration**: Train collaboratively by syncing runs in real time under shared team projects, reviewing teammates' experiments, and discussing findings online.

**7. ✉️ Effortless Sharing**: Share persistent URLs for individual runs with collaborators or embed live charts directly into online docs and notebooks.

**8. 💻 Self-Hosting Support**: Deploy in air-gapped or on-premises environments with the self-hosted Community Edition to view dashboards and manage experiments locally. See the [Self-Hosting Guide](#-self-hosting).

**9. 🔌 Extensible Plugin Ecosystem**: Extend SwanLab workflows via plugins such as [Lark Notifications](https://docs.swanlab.cn/en/plugin/notification-lark.html), [Slack Notifications](https://docs.swanlab.cn/en/plugin/notification-slack.html), [CSV Logger](https://docs.swanlab.cn/en/plugin/writer-csv.html), and more.

> \[!IMPORTANT]
>
> **Star this repository** to receive instant notifications for all new releases on GitHub! ⭐️

![star-us](readme_files/star-us.png)

<br>

## 📃 Online Demos

Explore live interactive demos on SwanLab:

|             [ResNet50 Cats vs. Dogs Classification][demo-cats-dogs]              |                            [YOLOv8 COCO128 Object Detection][demo-yolo]                            |
| :------------------------------------------------------------------------------: | :------------------------------------------------------------------------------------------------: |
|                   [![][demo-cats-dogs-image]][demo-cats-dogs]                    |                                 [![][demo-yolo-image]][demo-yolo]                                  |
| Track a ResNet50 image classification training run on the Cats vs. Dogs dataset. | Object detection on COCO128 with YOLOv8, tracking training hyperparameters and evaluation metrics. |

|                      [Qwen2 Instruction Fine-Tuning][demo-qwen2-sft]                       |                  [LSTM Google Stock Price Prediction][demo-google-stock]                   |
| :----------------------------------------------------------------------------------------: | :----------------------------------------------------------------------------------------: |
|                        [![][demo-qwen2-sft-image]][demo-qwen2-sft]                         |                     [![][demo-google-stock-image]][demo-google-stock]                      |
| Track supervised instruction fine-tuning (SFT) of the Qwen2 LLM for instruction following. | Train an LSTM time-series model on historical Google stock data to forecast future prices. |

|              [ResNeXt101 Audio Classification][demo-audio-classification]              |                  [Qwen2-VL Fine-Tuning on COCO][demo-qwen2-vl]                  |
| :------------------------------------------------------------------------------------: | :-----------------------------------------------------------------------------: |
|           [![][demo-audio-classification-image]][demo-audio-classification]            |                    [![][demo-qwen2-vl-image]][demo-qwen2-vl]                    |
| Iterative model exploration from ResNet to ResNeXt on audio classification benchmarks. | LoRA fine-tuning of the Qwen2-VL vision-language model on the COCO2014 dataset. |

|        [EasyR1 Multimodal LLM RL Training][demo-easyr1-rl]        |            [Qwen2.5-0.5B GRPO Training][demo-qwen2-grpo]            |
| :---------------------------------------------------------------: | :-----------------------------------------------------------------: |
|            [![][demo-easyr1-rl-image]][demo-easyr1-rl]            |            [![][demo-qwen2-grpo-image]][demo-qwen2-grpo]            |
| Multimodal LLM reinforcement learning using the EasyR1 framework. | GRPO reinforcement learning with Qwen2.5-0.5B on the GSM8K dataset. |

[More Examples](https://docs.swanlab.cn/en/examples/mnist.html)

<br>

## 🏁 Quick Start

### 1. Installation

```bash
pip install swanlab
```

<details><summary>Install from Source</summary>

To try out the latest bleeding-edge features, install directly from source:

```bash
# Option 1: Editable install
git clone https://github.com/SwanHubX/SwanLab.git
pip install -e .

# Option 2: Install via Git URL
pip install git+https://github.com/SwanHubX/SwanLab.git
```

</details>

<details><summary>Install Offline Dashboard Extension</summary>

[Offline Dashboard Documentation](https://docs.swanlab.cn/en/guide_cloud/offline-board.html)

```bash
pip install 'swanlab[dashboard]'
```

</details>

### 2. Log In & Configure API Key

1. Sign up for a free account at [SwanLab](https://swanlab.cn).

2. Log in and copy your API Key from **User Settings** > [API Key](https://swanlab.cn/settings).

3. Open your terminal and run:

```bash
swanlab login
```

Paste your API Key when prompted and press Enter to complete authentication.

### 3. Integrate SwanLab into Your Training Script

```python
import swanlab

# Initialize a new SwanLab run
swanlab.init(
    project="my-first-ml",
    config={"learning-rate": 0.003},
)

# Log metrics during training
for i in range(10):
    swanlab.log({"loss": i, "acc": i})
```

That's it! Head over to [SwanLab](https://swanlab.cn) to inspect your first experiment run.

<br>

## 💻 Self-Hosting

The self-hosted Community Edition supports running the SwanLab dashboard in offline or private environments.

![swanlab-kubernetes](./readme_files/swanlab-kubernetes.png)

Deployment guides:

- [🔗 Kubernetes Deployment Guide](https://docs.swanlab.cn/en/self_host/kubernetes/deploy.html)
- [🔗 Docker Deployment Guide](https://docs.swanlab.cn/en/self_host/docker/deploy.html)
- [🔗 Migrating from Docker to Kubernetes](https://docs.swanlab.cn/en/self_host/docker/migration-docker-kubernetes.html)

<br>

## 🔥 Case Studies

**Featured Open-Source Projects & Tutorials Using SwanLab:**

- [happy-llm](https://github.com/datawhalechina/happy-llm): Hands-on tutorial on building Large Language Models from scratch. ![GitHub Repo stars](https://img.shields.io/github/stars/datawhalechina/happy-llm)
- [self-llm](https://github.com/datawhalechina/self-llm): A practical cookbook for fast fine-tuning (full-parameter / LoRA) and deployment of open-source LLMs and MLLMs on Linux. ![GitHub Repo stars](https://img.shields.io/github/stars/datawhalechina/self-llm)
- [Minimind](https://github.com/jingyaogong/minimind): 🚀🚀 Train a 26M-parameter GPT from scratch in just 2 hours! ![GitHub Repo stars](https://img.shields.io/github/stars/jingyaogong/minimind)
- [unlock-deepseek](https://github.com/datawhalechina/unlock-deepseek): Deep-dive analysis, extensions, and reproductions of the DeepSeek model family. ![GitHub Repo stars](https://img.shields.io/github/stars/datawhalechina/unlock-deepseek)
- [Qwen3-SmVL](https://github.com/ShaohonChen/Qwen3-SmVL): Multimodal grafting and fine-tuning combining the SmolVLM2 vision encoder with Qwen3-0.6B. ![GitHub Repo stars](https://img.shields.io/github/stars/ShaohonChen/Qwen3-SmVL)
- [OPPO/Agent_Foundation_Models](https://github.com/OPPO-PersonalAI/Agent_Foundation_Models): End-to-end agent foundation models trained via multi-agent distillation and agentic RL. ![GitHub Repo stars](https://img.shields.io/github/stars/OPPO-PersonalAI/Agent_Foundation_Models)
- [Tree-GRPO](https://github.com/AMAP-ML/Tree-GRPO): [ICLR 2026] Tree Search for LLM Agent Reinforcement Learning. ![GitHub Repo stars](https://img.shields.io/github/stars/AMAP-ML/Tree-GRPO)
- [FIPO](https://github.com/qwenpilot/FIPO): Official implementation of FIPO, a value-free RL recipe for eliciting deeper reasoning from a clean base model. ![GitHub Repo stars](https://img.shields.io/github/stars/qwenpilot/FIPO)
- [llm-agent-rl-lab](https://github.com/KMnO4-zx/llm-agent-rl-lab): Reproducing and studying RL algorithms for LLM agents, including PPO, GRPO, GSPO, DAPO, OPD, and beyond. ![GitHub Repo stars](https://img.shields.io/github/stars/KMnO4-zx/llm-agent-rl-lab)

**Research Papers Powered by SwanLab:**

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

**Step-by-Step Tutorials:**

- [MNIST Handwritten Digit Recognition](https://docs.swanlab.cn/en/examples/mnist.html)
- [FashionMNIST Apparel Classification](https://docs.swanlab.cn/en/examples/fashionmnist.html)
- [CIFAR-10 Image Classification](https://docs.swanlab.cn/en/examples/cifar10.html)
- [ResNet Cats vs. Dogs Classification](https://docs.swanlab.cn/en/examples/cats_dogs_classification.html)
- [YOLO Object Detection](https://docs.swanlab.cn/en/examples/yolo.html)
- [U-Net Medical Image Segmentation](https://docs.swanlab.cn/en/examples/unet-medical-segmentation.html)
- [Audio Classification](https://docs.swanlab.cn/en/examples/audio_classification.html)
- [DQN Reinforcement Learning — CartPole](https://docs.swanlab.cn/en/examples/dqn_cartpole.html)
- [LSTM Google Stock Price Prediction](https://docs.swanlab.cn/en/examples/audio_classification.html)
- [BERT Text Classification](https://docs.swanlab.cn/en/examples/bert.html)
- [Stable Diffusion Text-to-Image Fine-Tuning](https://docs.swanlab.cn/en/examples/stable_diffusion.html)
- [LLM Pre-training from Scratch](https://docs.swanlab.cn/en/examples/pretrain_llm.html)
- [GLM-4 Instruction Fine-Tuning](https://docs.swanlab.cn/en/examples/glm4-instruct.html)
- [Qwen Downstream Task Fine-Tuning](https://docs.swanlab.cn/en/examples/qwen_finetune.html)
- [Named Entity Recognition (NER)](https://docs.swanlab.cn/en/examples/ner.html)
- [Qwen3 Medical Domain Fine-Tuning](https://docs.swanlab.cn/en/examples/qwen3-medical.html)
- [Qwen2-VL Vision-Language Model Fine-Tuning](https://docs.swanlab.cn/en/examples/qwen_vl_coco.html)
- [GRPO Reinforcement Learning for LLMs](https://docs.swanlab.cn/en/examples/qwen_grpo.html)
- [Training the Qwen3-SmVL-0.6B Multimodal Model](https://docs.swanlab.cn/en/examples/qwen3_smolvlm_muxi.html)
- [Getting Started with Embodied AI Using LeRobot](https://docs.swanlab.cn/en/examples/robot/lerobot-guide.html)
- [GLM-4.5-Air LoRA Fine-Tuning & SwanLab Visualization](https://github.com/datawhalechina/self-llm/blob/master/models/GLM-4.5-Air/03-GLM-4.5-Air-Lora%20%E5%8F%8A%20Swanlab%20%E5%8F%AF%E8%A7%86%E5%8C%96%E5%BE%AE%E8%B0%83.md)
- [Building a RAG Pipeline: Open-Sourcing the SwanLab Docs Assistant](https://docs.swanlab.cn/course/prompt_engineering_course/11-swanlab_rag/1.swanlab-rag.html)

🌟 Have a tutorial you'd like to showcase here? PRs are welcome!

<br>

## 🎮 Hardware Monitoring

SwanLab automatically captures **hardware specifications** and **real-time resource utilization** during training. Supported hardware platforms include:

| Hardware          | Spec Logging | Resource Telemetry | Probe Script                                                                                                                                   |
| ----------------- | ------------ | ------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| NVIDIA GPU        | ✅           | ✅                 | [nvidia.py](https://github.com/SwanHubX/SwanLab/blob/main/swanlab/sdk/internal/probe_python/hardware_vendor/accelerator/nvidia.py)             |
| AMD GPU           | ✅           | ✅                 | [amd.py](https://github.com/SwanHubX/SwanLab/blob/main/swanlab/sdk/internal/probe_python/hardware_vendor/accelerator/amd.py)                   |
| Huawei Ascend NPU | ✅           | ✅                 | [ascend.py](https://github.com/SwanHubX/SwanLab/blob/main/swanlab/sdk/internal/probe_python/hardware_vendor/accelerator/huawei.py)             |
| Cambricon MLU     | ✅           | ✅                 | [cambricon.py](https://github.com/SwanHubX/SwanLab/blob/main/swanlab/sdk/internal/probe_python/hardware_vendor/accelerator/cambricon.py)       |
| Kunlunxin XPU     | ✅           | ✅                 | [kunlunxin.py](https://github.com/SwanHubX/SwanLab/blob/main/swanlab/sdk/internal/probe_python/hardware_vendor/accelerator/kunlunxin.py)       |
| Moore Threads GPU | ✅           | ✅                 | [moorethreads.py](https://github.com/SwanHubX/SwanLab/blob/main/swanlab/sdk/internal/probe_python/hardware_vendor/accelerator/moorethreads.py) |
| MetaX GPU         | ✅           | ✅                 | [metax.py](https://github.com/SwanHubX/SwanLab/blob/main/swanlab/sdk/internal/probe_python/hardware_vendor/accelerator/metax.py)               |
| Iluvatar GPU      | ✅           | ✅                 | [iluvatar.py](https://github.com/SwanHubX/SwanLab/blob/main/swanlab/sdk/internal/probe_python/hardware_vendor/accelerator/iluvatar.py)         |
| Hygon DCU         | ✅           | ✅                 | [hygon.py](https://github.com/SwanHubX/SwanLab/blob/main/swanlab/sdk/internal/probe_python/hardware_vendor/accelerator/hygon.py)               |
| CPU               | ✅           | ✅                 | [cpu.py](https://github.com/SwanHubX/SwanLab/blob/main/swanlab/sdk/internal/probe_python/hardware_vendor/cpu.py)                               |
| System Memory     | ✅           | ✅                 | [memory.py](https://github.com/SwanHubX/SwanLab/blob/main/swanlab/sdk/internal/probe_python/hardware_vendor/memory.py)                         |

Want to add support for additional hardware accelerators? Feel free to open an Issue or submit a PR!

<br>

## 🚗 Framework Integrations

Pair SwanLab with your favorite ML stack!  
Below are the frameworks supported out of the box. Feel free to open an [Issue](https://github.com/swanhubx/swanlab/issues) to request additional integrations.

**Foundational Frameworks**

- [PyTorch](https://docs.swanlab.cn/en/guide_cloud/integration/integration-pytorch.html)
- [MindSpore](https://docs.swanlab.cn/en/guide_cloud/integration/integration-ascend.html)
- [Keras](https://docs.swanlab.cn/en/guide_cloud/integration/integration-keras.html)

**LLM Training Frameworks**

- [HuggingFace Transformers](https://docs.swanlab.cn/en/guide_cloud/integration/integration-huggingface-transformers.html)
- [LLaMA Factory](https://docs.swanlab.cn/en/guide_cloud/integration/integration-llama-factory.html)
- [MS-Swift](https://docs.swanlab.cn/en/guide_cloud/integration/integration-swift.html)
- [Unsloth](https://docs.swanlab.cn/en/guide_cloud/integration/integration-unsloth.html)
- [MLX-LM](https://docs.swanlab.cn/en/guide_cloud/integration/integration-mlx-lm.html)
- [torchtune](https://docs.swanlab.cn/en/guide_cloud/integration/integration-pytorch-torchtune.html)
- [PaddleNLP](https://docs.swanlab.cn/en/guide_cloud/integration/integration-paddlenlp.html)
- [Sentence Transformers](https://docs.swanlab.cn/en/guide_cloud/integration/integration-sentence-transformers.html)
- [XTuner](https://docs.swanlab.cn/en/guide_cloud/integration/integration-xtuner.html)
- [OpenMind](https://modelers.cn/docs/en/openmind-library/1.0.0/basic_tutorial/finetune/finetune_pt.html#install-swanlab)

**LLM Reinforcement Learning Frameworks**

- [veRL](https://docs.swanlab.cn/en/guide_cloud/integration/integration-verl.html)
- [HuggingFace TRL](https://docs.swanlab.cn/en/guide_cloud/integration/integration-huggingface-trl.html)
- [NVIDIA NeMo RL](https://docs.swanlab.cn/en/guide_cloud/integration/integration-nvidia-nemo-rl.html)
- [EasyR1](https://docs.swanlab.cn/en/guide_cloud/integration/integration-easyr1.html)
- [AReaL](https://docs.swanlab.cn/en/guide_cloud/integration/integration-areal.html)
- [ROLL](https://docs.swanlab.cn/en/guide_cloud/integration/integration-roll.html)

**Robotics & Embodied AI**

- [RLinf](https://docs.swanlab.cn/en/guide_cloud/integration/integration-rlinf.html)

**Text-to-Image / Video Generation**

- [DiffSynth Studio](https://docs.swanlab.cn/en/guide_cloud/integration/integration-diffsynth-studio.html)

**Deep Learning Training Engines**

- [PyTorch Lightning](https://docs.swanlab.cn/en/guide_cloud/integration/integration-pytorch-lightning.html)
- [MMEngine](https://docs.swanlab.cn/en/guide_cloud/integration/integration-mmengine.html)
- [fastai](https://docs.swanlab.cn/en/guide_cloud/integration/integration-fastai.html)

**Computer Vision**

- [Ultralytics](https://docs.swanlab.cn/en/guide_cloud/integration/integration-ultralytics.html)
- [MMDetection](https://docs.swanlab.cn/en/guide_cloud/integration/integration-mmdetection.html)
- [MMSegmentation](https://docs.swanlab.cn/en/guide_cloud/integration/integration-mmsegmentation.html)
- [PaddleDetection](https://docs.swanlab.cn/en/guide_cloud/integration/integration-paddledetection.html)
- [PaddleYOLO](https://docs.swanlab.cn/en/guide_cloud/integration/integration-paddleyolo.html)

**Classical Machine Learning / Gradient Boosting**

- [LightGBM](https://docs.swanlab.cn/en/guide_cloud/integration/integration-lightgbm.html)
- [XGBoost](https://docs.swanlab.cn/en/guide_cloud/integration/integration-xgboost.html)
- [CatBoost](https://docs.swanlab.cn/en/guide_cloud/integration/integration-catboost.html)

**Model Evaluation**

- [EvalScope](https://docs.swanlab.cn/en/guide_cloud/integration/integration-evalscope.html)

**Reinforcement Learning**

- [Stable-Baselines3](https://docs.swanlab.cn/en/guide_cloud/integration/integration-sb3.html)

**Other Tools & Orchestrators:**

- [TensorBoard](https://docs.swanlab.cn/en/guide_cloud/integration/integration-tensorboard.html)
- [Weights & Biases](https://docs.swanlab.cn/en/guide_cloud/integration/integration-wandb.html)
- [MLflow](https://docs.swanlab.cn/en/guide_cloud/integration/integration-mlflow.html)
- [HuggingFace Accelerate](https://docs.swanlab.cn/en/guide_cloud/integration/integration-huggingface-accelerate.html)
- [Ray](https://docs.swanlab.cn/en/guide_cloud/integration/integration-ray.html)
- [Hydra](https://docs.swanlab.cn/en/guide_cloud/integration/integration-hydra.html)
- [OmegaConf](https://docs.swanlab.cn/en/guide_cloud/integration/integration-omegaconf.html)
- [SpecForge](https://docs.swanlab.cn/en/guide_cloud/integration/integration-specforge.html)

[More Integrations](https://docs.swanlab.cn/en/guide_cloud/integration/integration-any-library.html)

<br>

## 🔌 Plugins & API

Extend SwanLab's capabilities and streamline your experiment workflows with plugins:

- [Build a Custom Plugin](https://docs.swanlab.cn/en/plugin/custom-plugin.html)
- [Email Notifications](https://docs.swanlab.cn/en/plugin/notification-email.html)
- [Lark Notifications](https://docs.swanlab.cn/en/plugin/notification-lark.html)
- [DingTalk Notifications](https://docs.swanlab.cn/en/plugin/notification-dingtalk.html)
- [WeCom (WXWork) Notifications](https://docs.swanlab.cn/en/plugin/notification-wxwork.html)
- [Discord Notifications](https://docs.swanlab.cn/en/plugin/notification-discord.html)
- [Slack Notifications](https://docs.swanlab.cn/en/plugin/notification-slack.html)
- [Bark Notifications](https://docs.swanlab.cn/en/plugin/notification-bark.html)
- [CSV Logger](https://docs.swanlab.cn/en/plugin/writer-csv.html)

Programmatic API:

- [OpenAPI Reference](https://docs.swanlab.cn/en/api/py-openapi.html)

<br>

## 🆚 Comparison with Alternatives

### TensorBoard vs. SwanLab

- **☁️ Cloud-Native Sync & Remote Access**:
  SwanLab syncs and stores training runs in the cloud, making it effortless to monitor progress remotely, manage historical projects, share experiment links, receive real-time alerts, and inspect runs across devices. TensorBoard is strictly a local, offline visualization viewer.

- **👥 Built-In Team Collaboration**:
  For multi-user, cross-team ML workflows, SwanLab provides shared workspaces to organize team projects, share permanent links, and collaborate on findings. TensorBoard is designed for single-user local inspection and lacks native collaboration features.

- **💻 Persistent, Centralized Dashboard**:
  Whether you train on a local workstation, an on-prem cluster, or public cloud GPU instances, SwanLab aggregates all runs into a single centralized dashboard. With TensorBoard, you have to manually copy and manage `tfevents` files across machines.

- **💪 Scalable Tabular Analytics**:
  SwanLab's interactive tables let you view, search, sort, and filter metrics and hyperparameters across thousands of model runs to quickly identify top-performing checkpoints. TensorBoard struggles to scale to large numbers of runs.

### Weights & Biases vs. SwanLab

- **Weights & Biases** is a proprietary, closed-source MLOps platform that requires cloud connectivity or paid enterprise licensing.

- **SwanLab** offers both a full-featured cloud service and an **open-source, free, self-hosted edition** suitable for air-gapped and private infrastructure.

<br>

## 👥 Community

### Related Repositories

- [self-hosted](https://github.com/swanhubx/self-hosted): Deployment scripts and configurations for self-hosted SwanLab.
- [SwanLab-Docs](https://github.com/swanhubx/swanlab-docs): Official documentation source repository.
- [SwanLab-Dashboard](https://github.com/swanhubx/swanlab-dashboard): Web frontend for the lightweight offline dashboard launched via `swanlab watch`.

### Community & Support

- [GitHub Issues](https://github.com/SwanHubX/SwanLab/issues): Report bugs or request features.
- [Email Support](mailto:contact@swanlab.cn): Reach out directly with feedback or questions.
- <a href="https://docs.swanlab.cn/en/guide_cloud/community/online-support.html">WeChat Community</a>: Discuss SwanLab usage and exchange the latest AI/ML insights.

### SwanLab README Badges

If you use SwanLab in your projects, showcase it by adding a SwanLab badge to your README:

[![][tracking-swanlab-shield]][tracking-swanlab-shield-link] [![][visualize-swanlab-shield]][visualize-swanlab-shield-link]

```markdown
[![](https://raw.githubusercontent.com/SwanHubX/assets/main/badge2.svg)](your experiment url)
[![](https://raw.githubusercontent.com/SwanHubX/assets/main/badge1.svg)](your experiment url)
```

More brand assets: [assets](https://github.com/SwanHubX/assets)

### Citing SwanLab

If SwanLab has been helpful in your research, please consider citing our work:

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

### Contributing to SwanLab

Interested in contributing? Please start by reading our [Contributing Guide](CONTRIBUTING.md).

We also deeply appreciate shout-outs on social media, community events, and conferences. Thank you for your support!

<br>

**Contributors**

<a href="https://github.com/swanhubx/swanlab/graphs/contributors">
  <img src="https://contrib.rocks/image?repo=swanhubx/swanlab" />
</a>

<br>

<img src="./readme_files/swanlab-and-user.png" width="50%" />

## 📃 License

This repository is licensed under the [Apache 2.0 License](https://github.com/SwanHubX/SwanLab/blob/main/LICENSE).

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
[wechat-shield-link]: https://docs.swanlab.cn/en/guide_cloud/community/online-support.html
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

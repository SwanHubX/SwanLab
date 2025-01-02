[![Overview](readme_files/swanlab-overview-new.png)](https://swanlab.cn/)

<div align="center">

<a href="https://swanlab.cn">🔥SwanLab Cloud</a> · <a href="https://docs.swanlab.cn">📃Document</a> · <a href="https://github.com/swanhubx/swanlab/issues">Report Issue</a> · <a href="https://geektechstudio.feishu.cn/share/base/form/shrcnyBlK8OMD0eweoFcc2SvWKc">Feedback</a> · <a href="https://docs.swanlab.cn/zh/guide_cloud/general/changelog.html">Changelog</a>

[![license][license-shield]][license-shield-link]
[![last-commit][last-commit-shield]][last-commit-shield-link]
[![pypi-version][pypi-version-shield]][pypi-version-shield-link]
[![pypi-downloads][pypi-downloads-shield]][pypi-downloads-shield-link]
[![issues][issues-shield]][issues-shield-link]
<br>
[![swanlab-cloud][swanlab-cloud-shield]][swanlab-cloud-shield-link]
[![wechat][wechat-shield]][wechat-shield-link]
[![colab][colab-shield]][colab-shield-link]

[中文](README.md) / English / [日本語](README_JP.md) / [Русский](README_RU.md)

👋 Join our [WeChat group](https://docs.swanlab.cn/zh/guide_cloud/community/online-support.html)

</div>

## TOC

- [🌟 Recent Updates](#-recent-updates)
- [👋🏻 Introduction](#-Introduction)
- [📃 Demo](#-demo)
- [🏁 Quick Start](#-quick-start)
- [💻 Self-hosted](#-self-hosted)
- [🚗 Integration](#-integration)
- [🆚 Comparison with Familiar Tools](#-comparison-with-familiar-tools)
- [👥 Community](#-community)
- [📃 License](#-license)

<br/>

## 🌟 Recent Updates

- 2025.01.01: Added **persistent smoothing** for line charts and drag-to-resize functionality for line charts, enhancing the chart browsing experience.

- 2024.12.22: We have completed integration with [LLaMA Factory](https://github.com/hiyouga/LLaMA-Factory), now you can use SwanLab in LLaMA Factory to **track and visualize large model fine-tuning experiments**, [Usage Guide](https://github.com/hiyouga/LLaMA-Factory?tab=readme-ov-file#use-swanlab-logger).

- 2024.12.15: **Hardware Monitoring (0.4.0)** feature launched, supporting system-level information recording and monitoring for CPU, NPU (Ascend), and GPU (Nvidia).

- 2024.12.06: Added integration with [LightGBM](https://docs.swanlab.cn/guide_cloud/integration/integration-lightgbm.html) and [XGBoost](https://docs.swanlab.cn/guide_cloud/integration/integration-xgboost.html); increased the line length limit for log recording.

- 2024.11.26: Environment tab - Hardware section now supports recognition of **Huawei Ascend NPU** and **Kunpeng CPU**; Cloud vendor section supports recognition of QingCloud **GritStone**.

- 2024.11.24: **Brand new UI launched** with redesigned design language, optimized interaction details, and dark mode support.

- 2024.11.11: [VSCode](https://marketplace.visualstudio.com/items?itemName=SwanLab.swanlab&ssr=false#overview) plugin launched; Added integration with [Keras](https://docs.swanlab.cn/guide_cloud/integration/integration-keras.html) framework; Added `run.public` method to support retrieving project name, experiment name, link and other information, [#732](https://github.com/swanhubx/swanlab/issues/732)

<br>

## 👋🏻 Introduction

SwanLab is an open-source, lightweight AI model training tracking and visualization tool, offering a platform for tracking, recording, comparing, and collaborating on experiments.

Designed for artificial intelligence researchers, SwanLab features a user-friendly Python API and an elegant UI interface, providing functionalities such as **training visualization, automatic logging, hyperparameter recording, experiment comparison, and multi-user collaboration**. With SwanLab, researchers can identify training issues through intuitive visual charts, find research inspiration by comparing multiple experiments, and break down team communication barriers through **online web sharing** and **organization-based collaborative training**, thereby enhancing organizational training efficiency.

Here is a list of its core features:

**1. 📊 Experiment Metrics & Hyperparameter Tracking**: Minimal code integration into your machine learning pipeline to track and record key training metrics.

- Supports **cloud** usage (similar to Weights & Biases), allowing you to view training progress anytime, anywhere. [How to view experiments on mobile](https://docs.swanlab.cn/guide_cloud/general/app.html)
- Supports **hyperparameter recording** and tabular display
- **Supported metadata types**: Scalar metrics, images, audio, text, ...
- **Supported chart types**: Line charts, media charts (images, audio, text), ...
- **Automatic background logging**: Logging, hardware environment, Git repository, Python environment, Python library list, project run directory

![](readme_files/overview-2.png)

**2. ⚡️ Comprehensive Framework Integration**: Over **30+** frameworks including PyTorch, 🤗HuggingFace Transformers, PyTorch Lightning, 🦙LLaMA Factory, MMDetection, Ultralytics, PaddleDetection, LightGBM, XGBoost, Keras, Tensorboard, Weights&Biases, OpenAI, Swift, XTuner, Stable Baseline3, and Hydra.

**3. 💻 Hardware Monitoring**: Supports real-time recording and monitoring of system-level hardware metrics for CPU, NPU (Ascend), GPU (Nvidia), and memory.

**4. 📦 Experiment Management**: Through a centralized dashboard designed for training scenarios, quickly get an overview and manage multiple projects and experiments.

**4. 🆚 Result Comparison**: Compare hyperparameters and results of different experiments through online tables and comparison charts to uncover iteration insights.

**5. 👥 Online Collaboration**: Collaborate with your team on training, supporting real-time synchronization of experiments under one project. You can view your team's training records online and provide feedback and suggestions based on results.

**6. ✉️ Share Results**: Copy and send persistent URLs to share each experiment, easily send to partners, or embed in online notes.

**7. 💻 Self-Hosting Support**: Supports offline environment usage; the self-hosted community edition also allows viewing dashboards and managing experiments.

> \[!IMPORTANT]
>
> **Star Us**, You will receive all release notifications from GitHub without any delay ~ ⭐️

![star-us](readme_files/star-us.png)

<br>

## 📃 Demo

Check out SwanLab's online demo：

| [ResNet50 Cats vs Dogs][demo-cats-dogs] | [Yolov8-COCO128][demo-yolo] |
| :--------: | :--------: |
| [![][demo-cats-dogs-image]][demo-cats-dogs] | [![][demo-yolo-image]][demo-yolo] |
| Track the image classification task of training a simple ResNet50 model on the cats and dogs dataset. | Perform object detection tasks using Yolov8 on the COCO128 dataset, tracking training hyperparameters and metrics. |

| [Qwen2 Instruction Finetune][demo-qwen2-sft] | [LSTM Google Stock Prediction][demo-google-stock] |
| :--------: | :--------: |
| [![][demo-qwen2-sft-image]][demo-qwen2-sft] | [![][demo-google-stock-image]][demo-google-stock] |
| Track the instruction fine-tuning training of the Qwen2 large language model, completing simple instruction following. | Train a simple LSTM model on the Google stock price dataset to predict future stock prices. |

| [ResNeXt101 Audio Classification][demo-audio-classification] | [Qwen2-VL-COCO2014][demo-qwen2-vl] |
| :--------: | :--------: |
| [![][demo-audio-classification-image]][demo-audio-classification] | [![][demo-qwen2-vl-image]][demo-qwen2-vl] |
| Progressive Experimental Process from ResNet to ResNeXt in Audio Classification Tasks | Fine-tuning Qwen2-VL Multimodal Large Model on the COCO2014 Dataset Using LoRA |


[More examples](https://docs.swanlab.cn/zh/examples/mnist.html)

<br>

## 🏁 Quick Start

### 1.Installation

```bash
pip install swanlab
```

### 2.Log in and get the API Key

1. **Free [Sign Up](https://swanlab.cn)**

2. **Log in to your account**, go to User Settings > [API Key](https://swanlab.cn/settings) and copy your API Key.

3. **Open your terminal and enter**:

```bash
swanlab login
```

When prompted, enter your API Key and press Enter to complete the login.

### 3. Integrate SwanLab with Your Code

```python
import swanlab

# Create a new SwanLab experiment
swanlab.init(
    project="my-first-ml",
    config={'learning-rate': 0.003}
)

# Log metrics
for i in range(10):
    swanlab.log({"loss": i})
```

**All set!** Visit [SwanLab](https://swanlab.cn) to see your first SwanLab experiment.

![MNIST](readme_files/readme-mnist.png)

<br>

## 💻 Self-hosted

The community edition supports offline viewing of SwanLab dashboards.

### Offline Experiment Tracking

Set the parameters `logir` and `mode` in swanlab.init to track experiments offline:

```python
...

swanlab.init(
    logdir='./logs',
    mode='local',
)

...
```

- The parameter `mode` is set to `local`, which disables synchronizing the experiment to the cloud.

- The setting of the parameter `logdir` is optional, and it specifies the location for saving SwanLab log files (by
  default saved in the `swanlog` folder).

- Log files will be created and updated during tracking of experiments, and launching offline dashboards will also be
  based on these log files.

Other parts are completely consistent with cloud usage.

### Open Offline Board

Open the terminal and use the following command to open a SwanLab dashboard:

```bash
swanlab watch ./logs
```

After the operation is completed, SwanLab will provide you with a local URL link (default
is [http://127.0.0.1:5092](http://127.0.0.1:5092)).

Visit this link to view the experiment offline in the browser dashboard.

<br>

## 🚗 Integration

Combine your favorite framework with
SwanLab, [More Integration](https://docs.swanlab.cn/zh/guide_cloud/integration/integration-pytorch-lightning.html).

<details>
  <summary>
    <strong>⚡️ PyTorch Lightning</strong>
  </summary>
  <br>

Create an instance using `SwanLabLogger` and pass it into the `logger` parameter of `Trainer` to enable SwanLab to
record training metrics.

```python
from swanlab.integration.pytorch_lightning import SwanLabLogger
import importlib.util
import os
import pytorch_lightning as pl
from torch import nn, optim, utils
from torchvision.datasets import MNIST
from torchvision.transforms import ToTensor

encoder = nn.Sequential(nn.Linear(28 * 28, 128), nn.ReLU(), nn.Linear(128, 3))
decoder = nn.Sequential(nn.Linear(3, 128), nn.ReLU(), nn.Linear(128, 28 * 28))


class LitAutoEncoder(pl.LightningModule):
    def __init__(self, encoder, decoder):
        super().__init__()
        self.encoder = encoder
        self.decoder = decoder

    def training_step(self, batch, batch_idx):
        # training_step defines the train loop.
        # it is independent of forward
        x, y = batch
        x = x.view(x.size(0), -1)
        z = self.encoder(x)
        x_hat = self.decoder(z)
        loss = nn.functional.mse_loss(x_hat, x)
        # Logging to SwanLab (if installed) by default
        self.log("train_loss", loss)
        return loss

    def test_step(self, batch, batch_idx):
        # test_step defines the test loop.
        # it is independent of forward
        x, y = batch
        x = x.view(x.size(0), -1)
        z = self.encoder(x)
        x_hat = self.decoder(z)
        loss = nn.functional.mse_loss(x_hat, x)
        # Logging to SwanLab (if installed) by default
        self.log("test_loss", loss)
        return loss

    def configure_optimizers(self):
        optimizer = optim.Adam(self.parameters(), lr=1e-3)
        return optimizer


# init the autoencoder
autoencoder = LitAutoEncoder(encoder, decoder)

# setup data
dataset = MNIST(os.getcwd(), train=True, download=True, transform=ToTensor())
train_dataset, val_dataset = utils.data.random_split(dataset, [55000, 5000])
test_dataset = MNIST(os.getcwd(), train=False, download=True, transform=ToTensor())

train_loader = utils.data.DataLoader(train_dataset)
val_loader = utils.data.DataLoader(val_dataset)
test_loader = utils.data.DataLoader(test_dataset)

swanlab_logger = SwanLabLogger(
    project="swanlab_example",
    experiment_name="example_experiment",
    cloud=False,
)

trainer = pl.Trainer(limit_train_batches=100, max_epochs=5, logger=swanlab_logger)

trainer.fit(model=autoencoder, train_dataloaders=train_loader, val_dataloaders=val_loader)
trainer.test(dataloaders=test_loader)

```

</details>

<details>
<summary>
  <strong> 🤗HuggingFace Transformers</strong>
</summary>

<br>

Create an instance using `SwanLabCallback` and pass it into the `callbacks` parameter of `Trainer` to enable SwanLab to
record training metrics.

```python
import evaluate
import numpy as np
import swanlab
from swanlab.integration.huggingface import SwanLabCallback
from datasets import load_dataset
from transformers import AutoModelForSequenceClassification, AutoTokenizer, Trainer, TrainingArguments


def tokenize_function(examples):
    return tokenizer(examples["text"], padding="max_length", truncation=True)


def compute_metrics(eval_pred):
    logits, labels = eval_pred
    predictions = np.argmax(logits, axis=-1)
    return metric.compute(predictions=predictions, references=labels)


dataset = load_dataset("yelp_review_full")

tokenizer = AutoTokenizer.from_pretrained("bert-base-cased")

tokenized_datasets = dataset.map(tokenize_function, batched=True)

small_train_dataset = tokenized_datasets["train"].shuffle(seed=42).select(range(1000))
small_eval_dataset = tokenized_datasets["test"].shuffle(seed=42).select(range(1000))

metric = evaluate.load("accuracy")

model = AutoModelForSequenceClassification.from_pretrained("bert-base-cased", num_labels=5)

training_args = TrainingArguments(
    output_dir="test_trainer",
    report_to="none",
    num_train_epochs=3,
    logging_steps=50,
)

swanlab_callback = SwanLabCallback(experiment_name="TransformersTest", cloud=False)

trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=small_train_dataset,
    eval_dataset=small_eval_dataset,
    compute_metrics=compute_metrics,
    callbacks=[swanlab_callback],
)

trainer.train()
```

</details>

<details>
<summary>
  <strong> MMEngine(MMDetection etc.)</strong>
</summary>
<br>

Integrate `SwanlabVisBackend` into MMEngine to enable automatic logging of training metrics by SwanLab.

Add the following code snippet to your MM config file to start training:

```python
custom_imports = dict(imports=["swanlab.integration.mmengine"], allow_failed_imports=False)

vis_backends = [
    dict(
        type="SwanlabVisBackend",
        save_dir="runs/swanlab",
        init_kwargs={
            "project": "swanlab-mmengine",
        },
    ),
]

visualizer = dict(
    type="Visualizer",
    vis_backends=vis_backends,
)
```

</details>

<details>
<summary>
  <strong> Ultralytics</strong>
</summary>
<br>

Integrating SwanLab into Ultralytics is very simple; you can use the `add_swanlab_callback` function:

```python
from ultralytics import YOLO
from swanlab.integration.ultralytics import add_swanlab_callback

model = YOLO("yolov8n.yaml")
model.load()

add_swanlab_callback(model)

model.train(
    data="./coco.yaml",
    epochs=50,
    imgsz=320,
)
```

</details>

<br>

## 🆚 Comparison with familiar tools

### Tensorboard vs SwanLab

- **☁️ Online Usage Support**:
  With SwanLab, training experiments can be conveniently synchronized and saved in the cloud, allowing for remote
  monitoring of training progress, managing historical projects, sharing experiment links, sending real-time
  notification messages, and viewing experiments across multiple devices. In contrast, TensorBoard is an offline
  experiment tracking tool.

- **👥 Collaborative Multi-user Environment**:
  SwanLab facilitates easy management of multi-person training projects and sharing of experiment links for
  collaborative machine learning across teams. It also enables cross-space communication and discussion. On the other
  hand, TensorBoard is primarily designed for individual use, making it difficult to collaborate and share experiments
  with multiple users.

- **💻 Persistent, Centralized Dashboard**:
  Regardless of where you are training your models, be it on a local computer, a lab cluster, or on public cloud GPU
  instances, your results are logged to the same centralized dashboard. Using TensorBoard, on the other hand, requires
  spending time copying and managing TFEvent files from different machines.

- **💪 More Powerful Tables**:
  SwanLab tables allow you to view, search, and filter results from various experiments, making it easy to review
  thousands of model versions to find the best-performing models for different tasks. TensorBoard is not well-suited for
  large-scale projects.

### Weights and Biases vs SwanLab

- Weights and Biases is an online-only, proprietary MLOps platform.

- Not only does SwanLab support online usage, but it also offers an open-source, free, and self-hosted version.

<br>

## 👥 Community

### Community and support

- [GitHub Issues](https://github.com/SwanHubX/SwanLab/issues)：Errors and issues encountered when using SwanLab
- [Email support](zeyi.lin@swanhub.co)：Feedback on issues with using SwanLab
- <a href="https://docs.swanlab.cn/guide_cloud/community/online-support.html">WeChat</a>：Discuss issues using SwanLab,
  share the latest AI technology.

### SwanLab README Badge

If you like to use SwanLab in your work, please add the SwanLab badge to your README:

[![swanlab](https://img.shields.io/badge/powered%20by-SwanLab-438440)](https://github.com/swanhubx/swanlab)

```
[![swanlab](https://img.shields.io/badge/powered%20by-SwanLab-438440)](https://github.com/swanhubx/swanlab)
```

### Citing SwanLab in the paper

If you find SwanLab helpful for your research journey, please consider citing in the following format:

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

### Contribute to SwanLab

Considering contributing to SwanLab? First, please take some time to read
the [Contribution Guidelines](CONTRIBUTING.md).

At the same time, we warmly welcome support for SwanLab through social media, events, and conference sharing. Thank you!

<br>

**Contributors**

<a href="https://github.com/swanhubx/swanlab/graphs/contributors">
  <img src="https://contrib.rocks/image?repo=swanhubx/swanlab" />
</a>

### Download Icon

[SwanLab-Icon-SVG](readme_files/swanlab-logo.svg)

<br>

## 📃 License

This repository follows the [Apache 2.0 License](https://github.com/SwanHubX/SwanLab/blob/main/LICENSE) open source
license.

[license-shield]: https://img.shields.io/github/license/SwanHubX/SwanLab.svg?color=brightgreen
[license-shield-link]: https://github.com/SwanHubX/SwanLab/blob/main/LICENSE
[last-commit-shield]: https://img.shields.io/github/last-commit/SwanHubX/SwanLab
[last-commit-shield-link]: https://github.com/SwanHubX/SwanLab/commits/main
[pypi-version-shield]: https://img.shields.io/pypi/v/swanlab?color=orange
[pypi-version-shield-link]: https://pypi.org/project/swanlab/
[pypi-downloads-shield]: https://static.pepy.tech/badge/swanlab
[pypi-downloads-shield-link]: https://pepy.tech/project/swanlab
[issues-shield]: https://img.shields.io/github/issues/swanhubx/swanlab
[issues-shield-link]: https://github.com/swanhubx/swanlab/issues
[swanlab-cloud-shield]: https://img.shields.io/badge/Product-SwanLab云端版-636a3f
[swanlab-cloud-shield-link]: https://swanlab.cn/
[wechat-shield]: https://img.shields.io/badge/WeChat-微信-4cb55e
[wechat-shield-link]: https://geektechstudio.feishu.cn/wiki/NIZ9wp5LRiSqQykizbGcVzUKnic
[colab-shield]: https://colab.research.google.com/assets/colab-badge.svg
[colab-shield-link]: https://colab.research.google.com/drive/1RWsrY_1bS8ECzaHvYtLb_1eBkkdzekR3?usp=sharing

[demo-cats-dogs]: https://swanlab.cn/@ZeyiLin/Cats_Dogs_Classification/runs/jzo93k112f15pmx14vtxf/chart
[demo-cats-dogs-image]: readme_files/example-catsdogs.png
[demo-yolo]: https://swanlab.cn/@ZeyiLin/ultratest/runs/yux7vclmsmmsar9ear7u5/chart
[demo-yolo-image]: readme_files/example-yolo.png
[demo-qwen2-sft]: https://swanlab.cn/@ZeyiLin/Qwen2-fintune/runs/cfg5f8dzkp6vouxzaxlx6/chart
[demo-qwen2-sft-image]: readme_files/example-qwen2.png
[demo-google-stock]:https://swanlab.cn/@ZeyiLin/Google-Stock-Prediction/charts
[demo-google-stock-image]: readme_files/example-lstm.png
[demo-audio-classification]:https://swanlab.cn/@ZeyiLin/Google-Stock-Prediction/charts
[demo-audio-classification-image]: readme_files/example-audio-classification.png
[demo-qwen2-vl]:https://swanlab.cn/@ZeyiLin/Google-Stock-Prediction/charts
[demo-qwen2-vl-image]: readme_files/example-qwen2-vl.jpg
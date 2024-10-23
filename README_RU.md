[![Обзор](readme_files/swanlab-overview-new.png)](https://swanlab.cn/)

<div align="center">

<a href="https://swanlab.cn">🔥Онлайн-версия SwanLab</a> · <a href="https://docs.swanlab.cn">📃 Документация</a> · <a href="https://github.com/swanhubx/swanlab/issues">Сообщить о проблеме</a> · <a href="https://geektechstudio.feishu.cn/share/base/form/shrcnyBlK8OMD0eweoFcc2SvWKc">Обратная связь</a> · <a href="https://docs.swanlab.cn/zh/guide_cloud/general/changelog.html">Журнал обновлений</a>

[![лицензия][license-shield]][license-shield-link]
[![последний-коммит][last-commit-shield]][last-commit-shield-link]
[![версия-pypi][pypi-version-shield]][pypi-version-shield-link]
[![загрузки-pypi][pypi-downloads-shield]][pypi-downloads-shield-link]
[![issues][issues-shield]][issues-shield-link]
<br>
[![swanlab-cloud][swanlab-cloud-shield]][swanlab-cloud-shield-link]
[![wechat][wechat-shield]][wechat-shield-link]
[![colab][colab-shield]][colab-shield-link]

[中文](README.md) / [English](README_EN.md) / [日本語](README_JP.md) / Русский

👋 Присоединяйтесь к нашему [чату в WeChat](https://geektechstudio.feishu.cn/wiki/NIZ9wp5LRiSqQykizbGcVzUKnic)

</div>

## Содержание

- [👋🏻 Что такое SwanLab](#-что-такое-swanlab)
- [📃 Онлайн демонстрация](#-онлайн-демонстрация)
- [🏁 Быстрый старт](#-быстрый-старт)
- [💻 Самостоятельное размещение](#-самостоятельное-размещение)
- [🚗 Интеграция с фреймворками](#-интеграция-с-фреймворками)
- [🆚 Сравнение с привычными инструментами](#-сравнение-с-привычными-инструментами)
- [👥 Сообщество](#-сообщество)
- [📃 Лицензия](#-лицензия)

<br/>

## 👋🏻 Что такое SwanLab

SwanLab — это инструмент для отслеживания экспериментов с ИИ с открытым исходным кодом и легким весом, предоставляющий платформу для отслеживания, сравнения и совместной работы над экспериментами.

SwanLab предлагает удобный API и привлекательный интерфейс, сочетающий такие возможности, как отслеживание гиперпараметров, запись метрик, онлайн-сотрудничество и совместное использование ссылок на эксперименты. Это позволяет быстро отслеживать эксперименты с ИИ, визуализировать процессы, записывать гиперпараметры и делиться ими с коллегами.

Ниже приведен список его основных особенностей:

**1. 📊 Отслеживание метрик экспериментов и гиперпараметров**: Минимальная инсталляция кода в ваш pipeleine машинного обучения позволяет отслеживать ключевые метрики

- Свободное ведение записей гиперпараметров и конфигурации экспериментов
- Поддерживаемые типы метаданных: скалярные метрики, изображения, аудио, текст и др.
- Поддерживаемые типы диаграмм: линейные графики, медиа-графики (изображения, аудио, текст) и др.
- Автоматическая запись: регистрирование в консоли, аппаратное обеспечение GPU, информация Git, список библиотек Python, директория кода

![](readme_files/overview-2.png)

**2. ⚡️ Весьма интегрировано с фреймворками**: PyTorch, Tensorflow, PyTorch Lightning, 🤗HuggingFace, Transformers, MMEngine, Ultralytics, fastai, Tensorboard, OpenAI, ZhipuAI, Hydra и др.

**3. 📦 Организация экспериментов**: Централизованная приборная панель для быстрой и простой организации нескольких проектов и экспериментов, чтобы взглянуть на общее состояние обучения

**4. 🆚 Сравнение результатов**: Сравнивайте гиперпараметры и результаты различных экспериментов с помощью онлайн табло и графиков, чтобы вдохновляться на итерации

**5. 👥 Онлайн сотрудничество**: Вы можете работать над проектом вместе с вашей командой, в реальном времени синхронизируя эксперименты в рамках проекта, просматривать журналы обучения вашей команды и высказывать свои мнения и предложения по результатам

**6. ✉️ Совместное использование результатов**: Копируйте и отправляйте постоянные URL-адреса, чтобы поделиться каждым экспериментом, легко отправляйте их коллегам или интегрируйте в онлайн-заметки

**7. 💻 Поддержка самостоятельного размещения**: Поддержка использования в режиме оффлайн, самостоятелный вариант развертывания также позволяет просматривать панель управления и управлять экспериментами

> \[!ВАЖНО]
>
> **Добавьте проект в избранное**, чтобы вы могли получать все уведомления о публикации на GitHub без задержек ～ ⭐️

![star-us](readme_files/star-us.png)

<br>

## 📃 Онлайн демонстрация

Посмотрите онлайн демонстрацию SwanLab:

|                    [ResNet50 Классификация по категориям "Кошки и собаки"](https://swanlab.cn/@ZeyiLin/Cats_Dogs_Classification/runs/jzo93k112f15pmx14vtxf/chart)                    |                [Yolov8-COCO128 Обнаружение объектов](https://swanlab.cn/@ZeyiLin/ultratest/runs/yux7vclmsmmsar9ear7u5/chart)                 |
| :----------------------------------------------------------------------------------------------------------------------------------------------: | :------------------------------------------------------------------------------------------------------------------------------: |
| <a href="https://swanlab.cn/@ZeyiLin/Cats_Dogs_Classification/runs/jzo93k112f15pmx14vtxf/chart"> <img src="readme_files/example-mnist.png"> </a> | <a href="https://swanlab.cn/@ZeyiLin/ultratest/runs/yux7vclmsmmsar9ear7u5/chart"> <img src="readme_files/example-yolo.png"> </a> |
|                                          Отслеживание задачи классификации изображений с использованием простой модели ResNet50 на датасете "Кошки и собаки".                                          |                             Обнаружение объектов с использованием Yolov8 на датасете COCO128, отслеживание гиперпараметров обучения и метрик.                              |

|                     [Тонкая настройка Qwen2](https://swanlab.cn/@ZeyiLin/Qwen2-fintune/runs/cfg5f8dzkp6vouxzaxlx6/chart)                      |                  [Прогнозирование акций Google с помощью LSTM](https://swanlab.cn/@ZeyiLin/Google-Stock-Prediction/charts)                  |
| :-----------------------------------------------------------------------------------------------------------------------------------: | :------------------------------------------------------------------------------------------------------------------: |
| <a href="https://swanlab.cn/@ZeyiLin/Qwen2-fintune/runs/cfg5f8dzkp6vouxzaxlx6/chart"> <img src="readme_files/example-qwen2.png"> </a> | <a href="https://swanlab.cn/@ZeyiLin/Google-Stock-Prediction/charts"> <img src="readme_files/example-lstm.png"> </a> |
|                                       Отслеживание тонкой настройки Qwen2 модели, выполняющей простые инструкции.                                       |                        Обучение модели LSTM на датасете цен акций Google для прогнозирования будущих значений цен.                        |

[Больше примеров](https://docs.swanlab.cn/zh/examples/mnist.html)

<br>

## 🏁 Быстрый старт

### 1. Установка

```bash
pip install swanlab
```

### 2. Вход и получение API ключа

1. Бесплатная [регистрация](https://swanlab.cn)

2. Войти, скопировать ваш API ключ из настроек пользователя > [API Key](https://swanlab.cn/settings)

3. Откройте терминал и введите:

```bash
swanlab login
```

Когда система запросит, введите ваш API ключ, нажмите Enter, чтобы завершить вход.

### 3. Интеграция SwanLab с вашим кодом

```python
import swanlab

# Инициализация нового эксперимента swanlab
swanlab.init(
    project="my-first-ml",
    config={'learning-rate': 0.003},
)

# Запись метрик
for i in range(10):
    swanlab.log({"loss": i, "acc": i})
```

Сделано! Перейдите на [SwanLab](https://swanlab.cn), чтобы увидеть ваш первый эксперимент инга SwanLab.

![MNIST](readme_files/readme-mnist.png)

<br>

## 💻 Самостоятельное размещение

Вариант самостоятельного размещения позволяет использовать панель инструментов SwanLab в оффлайне.

### Отслеживание экспериментов в режиме оффлайн

Установите параметры `logir` и `mode` в swanlab.init для оффлайн отслеживания:

```python
...

swanlab.init(
    logdir='./logs',
    mode='local',
)

...
```

- Параметр `mode` установить в `local`, чтобы отключить синхронизацию экспериментов с облаком

- Параметр `logdir` устанавливается по желанию, он определяет, где сохраняются журналы SwanLab (по умолчанию в папке `swanlog`)

  - Журналы создаются и обновляются в процессе отслеживания экспериментов, запуск оффлайн-доски также будет основан на этих журналах

Остальные части настроек идентичны облачному использованию.

### Запуск оффлайн-доски

Откройте терминал и выполните следующую команду для запуска панели управления SwanLab:

```bash
swanlab watch ./logs
```

После запуска, SwanLab предоставит вам локальную ссылку (по умолчанию [http://127.0.0.1:5092](http://127.0.0.1:5092))

Посетите эту ссылку, чтобы просматривать эксперименты через оффлайн доску.

<br>

## 🚗 Интеграция с фреймворками

Используйте фреймворки вместе с SwanLab, [узнайте больше об интеграции](https://docs.swanlab.cn/zh/guide_cloud/integration/integration-pytorch-lightning.html).

<details>
  <summary>
    <strong>⚡️ PyTorch Lightning</strong>
  </summary>
  <br>

Создайте `SwanLabLogger`, затем передайте его в параметр `logger` у `Trainer`, чтобы SwanLab записывал метрики обучения.

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
        # training_step определяет цикл обучения.
        # это не связано с forward
        x, y = batch
        x = x.view(x.size(0), -1)
        z = self.encoder(x)
        x_hat = self.decoder(z)
        loss = nn.functional.mse_loss(x_hat, x)
        # Логирование в SwanLab (если установлен) по умолчанию
        self.log("train_loss", loss)
        return loss

    def test_step(self, batch, batch_idx):
        # test_step определяет цикл тестирования.
        # это не связано с forward
        x, y = batch
        x = x.view(x.size(0), -1)
        z = self.encoder(x)
        x_hat = self.decoder(z)
        loss = nn.functional.mse_loss(x_hat, x)
        # Логирование в SwanLab (если установлен) по умолчанию
        self.log("test_loss", loss)
        return loss

    def configure_optimizers(self):
        optimizer = optim.Adam(self.parameters(), lr=1e-3)
        return optimizer


# инициализация автоэнкодера
autoencoder = LitAutoEncoder(encoder, decoder)

# настройка данных
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

Создайте `SwanLabCallback`, затем передайте его в параметр `callbacks` у `Trainer`, чтобы SwanLab записывал метрики обучения.

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
  <strong> MMEngine(MMDetection и др.)</strong>
</summary>
<br>

Интегрируйте `SwanlabVisBackend`, разработанный специально для MMEngine, чтобы автоматически записывать метрики обучения в SwanLab.

Добавьте следующий фрагмент кода в файл вашей конфигурации MM для начала обучения.

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

Интеграция SwanLab с Ultralytics очень проста; просто используйте функцию `add_swanlab_callback`:

```python
from ultralytics import YOLO
from swanlab.integration.ultralytics import add_swanlab_callback

model = YOLO("yolov8n.yaml")
model.load()

# Добавление обратного вызова swanlab
add_swanlab_callback(model)

model.train(
    data="./coco.yaml",
    epochs=50,
    imgsz=320,
)
```

</details>

<br>

## 🆚 Сравнение с привычными инструментами

### Tensorboard vs SwanLab

- **☁️ Поддержка онлайн использования**:
  С помощью SwanLab можно легко сохранять и синхронизировать эксперименты в облаке, что позволяет удаленно следить за прогрессом обучения, управлять историческими проектами, делиться ссылками на эксперименты, отправлять уведомления в реальном времени и просматривать эксперименты на разных устройствах. Tensorboard же является оффлайн-инструментом для отслеживания экспериментов.

- **👥 Мульти-пользовательская сотрудничество**:
  Во время многопользовательского, межкомандного сотрудничества в машинном обучении, с помощью SwanLab можно легко управлять проектами нескольких пользователей, делиться ссылками на эксперименты и общаться в разных условиях. Tensorboard в основном рассчитан на одиночное использование и плохо подходит для коллективной работы и совместного использования экспериментов.

- **💻 Постоянная, централизованная панель мониторинга**:
  Где вы ни обучали бы модели, будь то на локальной машине, в лабораторном кластере или в GPU-инстансах публичного облака, ваши результаты будут записаны в единую централизованную панель мониторинга. Используя TensorBoard, вам придется тратить время на копирование и управление TFEvent-файлами с разных машин.

- **💪 Более мощная таблица**:
  В SwanLab можно просматривать, искать и фильтровать результаты разных экспериментов с помощью таблицы. Это позволяет легко просматривать тысячи версий моделей и находить наилучшие для разных задач производительности модели. TensorBoard не подходит для крупных проектов.

### Weights and Biases vs SwanLab

- Weights and Biases — это проприетарная MLOps платформа, которую необходимо использовать только в режиме онлайн

- SwanLab поддерживает как онлайн использование, так и бесплатное, с открытым кодом, самостоятельное размещение

<br>

## 👥 Сообщество

### Сообщество и поддержка

- [GitHub Issues](https://github.com/SwanHubX/SwanLab/issues): ошибки и проблемы, возникающие при использовании SwanLab
- [Электронная почта поддержки](zeyi.lin@swanhub.co): обратная связь по вопросам использования SwanLab
- <a href="https://geektechstudio.feishu.cn/wiki/NIZ9wp5LRiSqQykizbGcVzUKnic">WeChat группа</a>: обсуждение вопросов использования SwanLab, обмен новинками в области технологий ИИ

### Бейдж для README SwanLab

Если вам нравится использовать SwanLab в своей работе, добавьте бейдж SwanLab в ваш README:

[![swanlab](https://img.shields.io/badge/powered%20by-SwanLab-438440)](https://github.com/swanhubx/swanlab)

```
[![swanlab](https://img.shields.io/badge/powered%20by-SwanLab-438440)](https://github.com/swanhubx/swanlab)
```

### Цитирование SwanLab в публикациях

Если SwanLab оказалась полезной для ваших научных изысканий, пожалуйста, цитируйте ее в следующем формате:

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

### Вклад в развитие SwanLab

Задумываетесь о вкладе в развитие SwanLab? Сначала потратьте немного времени, чтобы ознакомиться с [руководством по внесению вклада](CONTRIBUTING.md).

Мы также очень признательны за любую поддержку в виде распространения информации о SwanLab через социальные сети, участие в мероприятиях и конференциях, огромное спасибо!

### Загрузить Icon

[SwanLab-Icon-SVG](readme_files/swanlab-logo.svg)

<br>

**Участники**

<a href="https://github.com/swanhubx/swanlab/graphs/contributors">
  <img src="https://contrib.rocks/image?repo=swanhubx/swanlab" />
</a>

<br>

## 📃 Лицензия

Этот репозиторий лицензирован в соответствии с [лицензией Apache 2.0](https://github.com/SwanHubX/SwanLab/blob/main/LICENSE)

<!-- ссылки -->

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
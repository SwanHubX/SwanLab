<div align="center">

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="readme_files/swanlab-logo-single-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="readme_files/swanlab-logo-single.svg">
  <img alt="SwanLab" src="readme_files/swanlab-logo-single.svg" width="70" height="70">
</picture>

<h1>SwanLab</h1>

Открытый инструмент для отслеживания и визуализации обучения глубоких нейронных сетей с современным дизайном  
Поддерживает как облачное, так и оффлайн использование, совместим с 30+ популярными фреймворками, легко интегрируется с вашим кодом экспериментов

<a href="https://swanlab.cn">🔥SwanLab Online</a> · <a href="https://docs.swanlab.cn">📃 Документация</a> · <a href="https://github.com/swanhubx/swanlab/issues">Сообщить о проблеме</a> · <a href="https://geektechstudio.feishu.cn/share/base/form/shrcnyBlK8OMD0eweoFcc2SvWKc">Оставить отзыв</a> · <a href="https://docs.swanlab.cn/en/guide_cloud/general/changelog.html">История изменений</a>

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

[中文](README.md) / [English](README_EN.md) / [日本語](README_JP.md) / Русский

👋 Присоединяйтесь к нашей [группе в WeChat](https://docs.swanlab.cn/en/guide_cloud/community/online-support.html)

</div>


## Содержание

- [🌟 Последние обновления](#-последние-обновления)
- [👋🏻 Что такое SwanLab](#-что-такое-swanlab)
- [📃 Онлайн-демонстрация](#-онлайн-демонстрация)
- [🏁 Быстрый старт](#-быстрый-старт)
- [💻 Самостоятельное размещение](#-самостоятельное-размещение)
- [🚗 Интеграция с фреймворками](#-интеграция-с-фреймворками)
- [🔌 Плагины](#-Плагины)
- [🆚 Сравнение с известными инструментами](#-сравнение-с-известными-инструментами)
- [👥 Сообщество](#-сообщество)
- [📃 Лицензия](#-лицензия)

<br/>

## 🌟 Последние обновления

• 2025.03.30: Добавлена поддержка метода **swanlab.Settings**, обеспечивающего более точный контроль над поведением экспериментов; добавлена поддержка мониторинга оборудования **Cambricon MLU**; интегрированы уведомления через [Slack](https://docs.swanlab.cn/plugin/notification-slack.html) и [Discord](https://docs.swanlab.cn/plugin/notification-discord.html).

- 2025.03.21: 🎉🤗 HuggingFace Transformers официально интегрировал SwanLab (версия >=4.50.0), [#36433](https://github.com/huggingface/transformers/pull/36433). Добавлена поддержка графиков Object3D, теперь вы можете отслеживать и визуализировать 3D облака точек, [docs](https://docs.swanlab.cn/en/api/py-object3d.html). Мониторинг оборудования поддерживает запись памяти GPU (МБ), использования диска, а также отправленных и полученных данных сети.

- 2025.03.12: 🎉🎉 SwanLab**самостоятельное размещение** теперь доступно! [🔗 Документация](https://docs.swanlab.cn/en/guide_cloud/self_host/docker-deploy.html); SwanLab теперь поддерживает расширение плагинов, таких как [Email Notification](https://docs.swanlab.cn/en/plugin/notification-email.html) и [Lark Notification](https://docs.swanlab.cn/en/plugin/notification-lark.html).

- 2025.03.09: Добавлена интеграция с [MLFlow](https://docs.swanlab.cn/en/guide_cloud/integration/integration-mlflow.html), теперь вы можете использовать SwanLab для **отслеживания и визуализации экспериментов MLFlow**.

- 2025.03.06: Мы объединили усилия с [DiffSynth Studio](https://github.com/modelscope/diffsynth-studio), теперь вы можете использовать SwanLab в DiffSynth Studio для **отслеживания и визуализации экспериментов по тонкой настройке больших моделей** [инструкция по использованию](https://docs.swanlab.cn/en/guide_cloud/integration/integration-diffsynth-studio.html).

- 2025.03.04: Добавлена интеграция с [MLFlow](https://docs.swanlab.cn/en/guide_cloud/integration/integration-mlflow.html), теперь вы можете использовать SwanLab для **отслеживания и визуализации экспериментов MLFlow**.

- 2025.03.01: Добавлена функция **перемещения экспериментов**, теперь вы можете перемещать эксперименты между организациями и проектами.

- 2025.02.24: Мы объединили усилия с [EasyR1](https://github.com/hiyouga/EasyR1), теперь вы можете использовать SwanLab в EasyR1 для **отслеживания и визуализации экспериментов по тонкой настройке больших моделей** [инструкция по использованию](https://github.com/hiyouga/EasyR1?tab=readme-ov-file#merge-checkpoint-in-hugging-face-format).

- 2025.02.18: Мы объединили усилия с [Swift](https://github.com/modelscope/ms-swift), теперь вы можете использовать SwanLab в Swift's CLI/WebUI для **отслеживания и визуализации экспериментов по тонкой настройке больших моделей** [инструкция по использованию](https://docs.swanlab.cn/en/guide_cloud/integration/integration-swift.html).


<details><summary>Полный список изменений</summary>

- 2025.02.16: Добавлены функции перемещения групп и создания групп.

- 2025.02.09: Мы объединили усилия с [veRL](https://github.com/volcengine/verl), теперь вы можете использовать SwanLab в veRL для **отслеживания и визуализации экспериментов по тонкой настройке больших моделей** в LLaMA Factory, [инструкция по использованию](https://docs.swanlab.cn/en/guide_cloud/integration/integration-verl.html).

- 2025.02.05: `swanlab.log` поддерживает вложенные словари [#812](https://github.com/SwanHubX/SwanLab/pull/812), поддерживает параметры `name` и `notes`.

- 2025.01.22: Добавлена функция `sync_tensorboardX` и `sync_tensorboard_torch`, поддерживающая синхронизацию с этими двумя TensorBoard фреймворками.

- 2025.01.17: Добавлена функция `sync_wandb`, [документация](https://docs.swanlab.cn/en/guide_cloud/integration/integration-wandb.html), поддерживающая синхронизацию с отслеживанием экспериментов Weights & Biases; значительно улучшена производительность рендеринга логов.

- 2025.01.11: Облачная версия значительно оптимизировала производительность таблиц проектов и добавила поддержку таких функций, как перетаскивание, сортировка и фильтрация.

- 2025.01.01: Добавлено **сглаживание графиков** и возможность изменения размера графиков перетаскиванием, улучшено взаимодействие с графиками.

- 2024.12.22: Интеграция с [LLaMA Factory](https://github.com/hiyouga/LLaMA-Factory), теперь можно использовать SwanLab для **отслеживания и визуализации экспериментов по тонкой настройке больших моделей** в LLaMA Factory, [инструкция по использованию](https://github.com/hiyouga/LLaMA-Factory?tab=readme-ov-file#use-swanlab-logger).

- 2024.12.15: Добавлена функция **мониторинга оборудования (0.4.0)**, поддерживается запись и мониторинг системной информации для CPU, NPU (Ascend), GPU (Nvidia).

- 2024.12.06: Добавлена интеграция с [LightGBM](https://docs.swanlab.cn/en/guide_cloud/integration/integration-lightgbm.html) и [XGBoost](https://docs.swanlab.cn/en/guide_cloud/integration/integration-xgboost.html); увеличено ограничение на длину строки в логах.

- 2024.11.26: Вкладка "Окружение" — раздел "Оборудование" теперь поддерживает распознавание **Huawei Ascend NPU** и **Kunpeng CPU**; раздел "Облачные провайдеры" поддерживает распознавание **QingCloud**.

</details>

<br>

## 👋🏻 Что такое SwanLab

SwanLab — это открытый и легковесный инструмент для отслеживания и визуализации обучения моделей искусственного интеллекта, предоставляющий платформу для отслеживания, записи, сравнения и совместной работы над экспериментами.

SwanLab ориентирован на исследователей в области ИИ, предлагая удобный Python API и красивый интерфейс, а также функции **визуализации обучения, автоматической записи логов, записи гиперпараметров, сравнения экспериментов и совместной работы**. С помощью SwanLab исследователи могут обнаруживать проблемы в обучении на основе наглядных графиков, сравнивать несколько экспериментов для поиска идей и делиться результатами через **онлайн-страницы** и **совместное обучение в организациях**, что упрощает коммуникацию в команде и повышает эффективность обучения.

Основные функции:

**1. 📊 Отслеживание метрик и гиперпараметров**: Минималистичный код для встраивания в ваш ML pipeline, отслеживание ключевых метрик обучения.

- Поддержка **облачного** использования (аналогично Weights & Biases), возможность просмотра прогресса обучения в любое время. [Как смотреть эксперименты на телефоне](https://docs.swanlab.cn/en/guide_cloud/general/app.html)
- Поддержка **записи гиперпараметров** и их отображения в таблицах.
- **Поддерживаемые типы данных**: скалярные метрики, изображения, аудио, текст, 3D точки, ...
- **Поддерживаемые типы графиков**: линейные графики, медиа-графики (изображения, аудио, текст, 3D точки), ...
- **Автоматическая запись логов**: логирование, информация об оборудовании, Git-репозитории, окружение Python, список библиотек Python, рабочая директория проекта.

**2. ⚡️ Полная интеграция с фреймворками**: PyTorch, 🤗HuggingFace Transformers, PyTorch Lightning, 🦙LLaMA Factory, MMDetection, Ultralytics, PaddleDetetion, LightGBM, XGBoost, Keras, Tensorboard, Weights&Biases, OpenAI, Swift, XTuner, Stable Baseline3, Hydra и более **30+** фреймворков.

![](readme_files/integrations.png)

**3. 💻 Мониторинг оборудования**: Поддержка записи и мониторинга системных показателей CPU, NPU (Ascend), GPU (Nvidia), памяти.

**4. 📦 Управление экспериментами**: Централизованная панель управления, разработанная для сценариев обучения, позволяет быстро просматривать и управлять несколькими проектами и экспериментами.

**5. 🆚 Сравнение результатов**: Сравнение гиперпараметров и результатов разных экспериментов через онлайн-таблицы и графики, поиск идей для улучшения.

**6. 👥 Онлайн-сотрудничество**: Возможность совместного обучения с командой, синхронизация экспериментов в реальном времени в одном проекте, просмотр записей обучения команды и обсуждение результатов.

**7. ✉️ Поделиться результатами**: Копирование и отправка постоянных URL для обмена каждым экспериментом, удобная отправка партнерам или встраивание в онлайн-заметки.

**8. 💻 Поддержка самостоятельного размещения**: Поддержка оффлайн использования, локальная версия также позволяет просматривать панель управления и управлять экспериментами.

> \[!IMPORTANT]
>
> **Добавьте проект в избранное**, чтобы получать уведомления о всех новых выпусках без задержек～ ⭐️

![star-us](readme_files/star-us.png)

<br>

## 📃 Онлайн-демонстрация

Ознакомьтесь с онлайн-демонстрацией SwanLab:

| [Классификация кошек и собак с ResNet50][demo-cats-dogs] | [Обнаружение объектов с Yolov8-COCO128][demo-yolo] |
| :--------: | :--------: |
| [![][demo-cats-dogs-image]][demo-cats-dogs] | [![][demo-yolo-image]][demo-yolo] |
| Отслеживание простой модели ResNet50 для задачи классификации изображений на наборе данных кошек и собак. | Использование Yolov8 для задачи обнаружения объектов на наборе данных COCO128, отслеживание гиперпараметров и метрик обучения. |

| [Тонкая настройка Qwen2][demo-qwen2-sft] | [Прогнозирование акций Google с LSTM][demo-google-stock] |
| :--------: | :--------: |
| [![][demo-qwen2-sft-image]][demo-qwen2-sft] | [![][demo-google-stock-image]][demo-google-stock] |
| Отслеживание тонкой настройки большой языковой модели Qwen2 для выполнения простых инструкций. | Использование простой модели LSTM для прогнозирования будущих цен акций Google на наборе данных акций. |

| [Классификация аудио с ResNeXt101][demo-audio-classification] | [Тонкая настройка Qwen2-VL на наборе данных COCO][demo-qwen2-vl] |
| :--------: | :--------: |
| [![][demo-audio-classification-image]][demo-audio-classification] | [![][demo-qwen2-vl-image]][demo-qwen2-vl] |
| Постепенный процесс экспериментов от ResNet к ResNeXt в задаче классификации аудио. | Тонкая настройка мультимодальной модели Qwen2-VL на наборе данных COCO2014 с использованием Lora. |

[Больше примеров](https://docs.swanlab.cn/en/examples/mnist.html)

<br>

## 🏁 Быстрый старт

### 1. Установка

```bash
pip install swanlab
```

<details><summary>Установка из исходного кода</summary>

Если вы хотите испытать новейшие функции, вы можете установить программу из исходного кода.

```bash
# Method 1
git clone https://github.com/SwanHubX/SwanLab.git
pip install -e .

# Method 2
pip install git+https://github.com/SwanHubX/SwanLab.git
```

</details>

<details><summary>Установка расширения панели управления</summary>

[Документация по расширению панели управления](https://docs.swanlab.cn/en/guide_cloud/self_host/offline-board.html)

```bash
pip install 'swanlab[dashboard]'
```

</details>

### 2. Вход и получение API Key

1. Бесплатная [регистрация аккаунта](https://swanlab.cn)

2. Войдите в аккаунт, скопируйте ваш API Key в разделе пользовательских настроек > [API Key](https://swanlab.cn/settings)

3. Откройте терминал и введите:

```bash
swanlab login
```

При появлении запроса введите ваш API Key, нажмите Enter, чтобы завершить вход.

### 3. Интеграция SwanLab с вашим кодом

```python
import swanlab

# Инициализация нового эксперимента SwanLab
swanlab.init(
    project="my-first-ml",
    config={'learning-rate': 0.003},
)

# Запись метрик
for i in range(10):
    swanlab.log({"loss": i, "acc": i})
```

Готово! Перейдите на [SwanLab](https://swanlab.cn), чтобы увидеть ваш первый эксперимент.

<br>

## 💻 Самостоятельное размещение

Самостоятельная версия для сообщества поддерживает офлайн-просмотр панели управления SwanLab.

![swanlab-docker](./readme_files/swanlab-docker.png)

### 1. Развертывание самостоятельной версии с использованием Docker

Подробные инструкции см. в: [Документация](https://docs.swanlab.cn/en/guide_cloud/self_host/docker-deploy.html)

```bash
git clone https://github.com/SwanHubX/self-hosted.git
cd self-hosted/docker
```

Быстрая установка для Китая:

```bash
./install.sh
```

Скачивание и установка образа из DockerHub:

```bash
./install-dockerhub.sh
```

### 2. Направление экспериментов в самостоятельный сервис

Вход в самостоятельный сервис:

```bash
swanlab login --host http://localhost:8000
```

После входа вы можете записывать эксперименты в самостоятельный сервис.

<br>

## 🚗 Интеграция с фреймворками

Используйте ваш любимый фреймворк вместе с SwanLab!  
Ниже приведен список уже интегрированных фреймворков. Если вы хотите предложить интеграцию с другим фреймворком, создайте [Issue](https://github.com/swanhubx/swanlab/issues).

**Основные фреймворки**
- [PyTorch](https://docs.swanlab.cn/en/guide_cloud/integration/integration-pytorch.html)
- [MindSpore](https://docs.swanlab.cn/en/guide_cloud/integration/integration-ascend.html)
- [Keras](https://docs.swanlab.cn/en/guide_cloud/integration/integration-keras.html)

**Специализированные/фреймворки для тонкой настройки**
- [PyTorch Lightning](https://docs.swanlab.cn/en/guide_cloud/integration/integration-pytorch-lightning.html)
- [HuggingFace Transformers](https://docs.swanlab.cn/en/guide_cloud/integration/integration-huggingface-transformers.html)
- [LLaMA Factory](https://docs.swanlab.cn/en/guide_cloud/integration/integration-llama-factory.html)
- [Modelscope Swift](https://docs.swanlab.cn/en/guide_cloud/integration/integration-swift.html)
- [DiffSynth Studio](https://docs.swanlab.cn/en/guide_cloud/integration/integration-diffsynth-studio.html)
- [Sentence Transformers](https://docs.swanlab.cn/en/guide_cloud/integration/integration-sentence-transformers.html)
- [OpenMind](https://modelers.cn/docs/zh/openmind-library/1.0.0/basic_tutorial/finetune/finetune_pt.html#%E8%AE%AD%E7%BB%83%E7%9B%91%E6%8E%A7)
- [Torchtune](https://docs.swanlab.cn/en/guide_cloud/integration/integration-pytorch-torchtune.html)
- [XTuner](https://docs.swanlab.cn/en/guide_cloud/integration/integration-xtuner.html)
- [MMEngine](https://docs.swanlab.cn/en/guide_cloud/integration/integration-mmengine.html)
- [FastAI](https://docs.swanlab.cn/en/guide_cloud/integration/integration-fastai.html)
- [LightGBM](https://docs.swanlab.cn/en/guide_cloud/integration/integration-lightgbm.html)
- [XGBoost](https://docs.swanlab.cn/en/guide_cloud/integration/integration-xgboost.html)


**Компьютерное зрение**
- [Ultralytics](https://docs.swanlab.cn/en/guide_cloud/integration/integration-ultralytics.html)
- [MMDetection](https://docs.swanlab.cn/en/guide_cloud/integration/integration-mmdetection.html)
- [MMSegmentation](https://docs.swanlab.cn/en/guide_cloud/integration/integration-mmsegmentation.html)
- [PaddleDetection](https://docs.swanlab.cn/en/guide_cloud/integration/integration-paddledetection.html)
- [PaddleYOLO](https://docs.swanlab.cn/en/guide_cloud/integration/integration-paddleyolo.html)

**Обучение с подкреплением**
- [Stable Baseline3](https://docs.swanlab.cn/en/guide_cloud/integration/integration-sb3.html)
- [veRL](https://docs.swanlab.cn/en/guide_cloud/integration/integration-verl.html)
- [HuggingFace trl](https://docs.swanlab.cn/en/guide_cloud/integration/integration-huggingface-trl.html)
- [EasyR1](https://docs.swanlab.cn/en/guide_cloud/integration/integration-easyr1.html)

**Другие фреймворки:**
- [Tensorboard](https://docs.swanlab.cn/en/guide_cloud/integration/integration-tensorboard.html)
- [Weights&Biases](https://docs.swanlab.cn/en/guide_cloud/integration/integration-wandb.html)
- [MLFlow](https://docs.swanlab.cn/en/guide_cloud/integration/integration-mlflow.html)
- [HuggingFace Accelerate](https://docs.swanlab.cn/en/guide_cloud/integration/integration-huggingface-accelerate.html)
- [Unsloth](https://docs.swanlab.cn/en/guide_cloud/integration/integration-unsloth.html)
- [Hydra](https://docs.swanlab.cn/en/guide_cloud/integration/integration-hydra.html)
- [Omegaconf](https://docs.swanlab.cn/en/guide_cloud/integration/integration-omegaconf.html)
- [OpenAI](https://docs.swanlab.cn/en/guide_cloud/integration/integration-openai.html)
- [ZhipuAI](https://docs.swanlab.cn/en/guide_cloud/integration/integration-zhipuai.html)

[Больше интеграций](https://docs.swanlab.cn/en/guide_cloud/integration/integration-pytorch-lightning.html)

<br>

## 🆚 Сравнение с известными инструментами

### Tensorboard vs SwanLab

- **☁️ Поддержка онлайн-использования**:
  SwanLab позволяет легко синхронизировать и сохранять эксперименты в облаке, что удобно для удаленного просмотра прогресса обучения, управления историей проектов, обмена ссылками на эксперименты, отправки уведомлений и просмотра экспериментов на разных устройствах. Tensorboard — это оффлайн инструмент для отслеживания экспериментов.

- **👥 Совместная работа**:
  При совместной работе над проектами машинного обучения SwanLab упрощает управление проектами, обмен ссылками на эксперименты и обсуждение результатов. Tensorboard в основном предназначен для индивидуального использования и не поддерживает совместную работу.

- **💻 Постоянная и централизованная панель управления**:
  Независимо от того, где вы обучаете модель — на локальном компьютере, в лабораторном кластере или на облачном GPU, результаты будут записываться в одну централизованную панель управления. Tensorboard требует ручного копирования и управления файлами TFEvent с разных машин.

- **💪 Более мощные таблицы**:
  SwanLab позволяет просматривать, искать и фильтровать результаты из разных экспериментов, что упрощает поиск лучшей модели для различных задач. Tensorboard не подходит для крупных проектов.

### Weights and Biases vs SwanLab

- Weights and Biases — это закрытая MLOps платформа, требующая подключения к интернету.

- SwanLab поддерживает как онлайн, так и оффлайн использование, а также предоставляет открытую и бесплатную версию для самостоятельного размещения.

<br>

## 🔌 Плагины

Расширьте функциональность SwanLab и улучшите управление экспериментами с помощью плагинов!

- [Настройте свой плагин](https://docs.swanlab.cn/en/plugin/custom-plugin.html)
-  [Уведомления по электронной почте](https://docs.swanlab.cn/en/plugin/notification-email.html)
-  [Уведомления в Lark](https://docs.swanlab.cn/en/plugin/notification-lark.html)
-  [Уведомления в DingTalk](https://docs.swanlab.cn/en/plugin/notification-dingtalk.html)
-  [Уведомления в WXWork](https://docs.swanlab.cn/en/plugin/notification-wxwork.html)
-  [Уведомления в Discord](https://docs.swanlab.cn/en/plugin/notification-discord.html)
-  [Уведомления в Slack](https://docs.swanlab.cn/en/plugin/notification-slack.html)
-  [CSV-логгер](https://docs.swanlab.cn/en/plugin/writer-csv.html)

<br>

## 👥 Сообщество

### Периферийные репозитории

• [SwanLab-Docs](https://github.com/swanhubx/swanlab-docs): Официальный репозиторий документации.
• [SwanLab-Dashboard](https://github.com/swanhubx/swanlab-dashboard): Репозиторий оффлайн-дашборда, содержащий веб-код для легковесного оффлайн-дашборда, открываемого командой `swanlab watch`.
• [self-hosted](https://github.com/swanhubx/self-hosted): Репозиторий скриптов для частного развертывания.

### Сообщество и поддержка

- [GitHub Issues](https://github.com/SwanHubX/SwanLab/issues): Ошибки и проблемы при использовании SwanLab.
- [Электронная почта](zeyi.lin@swanhub.co): Отправка отзывов и вопросов по использованию SwanLab.
- <a href="https://docs.swanlab.cn/en/guide_cloud/community/online-support.html">Группа в WeChat</a>: Обсуждение вопросов по использованию SwanLab, обмен новыми технологиями в области ИИ.

### Значок SwanLab для README

Если вам нравится использовать SwanLab в вашей работе, добавьте значок SwanLab в ваш README:

[![][tracking-swanlab-shield]][tracking-swanlab-shield-link]、[![][visualize-swanlab-shield]][visualize-swanlab-shield-link]

```
[![](https://raw.githubusercontent.com/SwanHubX/assets/main/badge2.svg)](your experiment url)
[![](https://raw.githubusercontent.com/SwanHubX/assets/main/badge1.svg)](your experiment url)
```

Больше дизайнерских материалов：[assets](https://github.com/SwanHubX/assets)

### Цитирование SwanLab в научных работах

Если SwanLab помог вам в ваших исследованиях, рассмотрите возможность цитирования в следующем формате:

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

Хотите внести вклад в SwanLab? Сначала ознакомьтесь с [руководством по вкладу](CONTRIBUTING.md).

Мы также приветствуем поддержку через социальные сети, мероприятия и конференции. Спасибо!

<br>

**Участники**

<a href="https://github.com/swanhubx/swanlab/graphs/contributors">
  <img src="https://contrib.rocks/image?repo=swanhubx/swanlab" />
</a>

<br>

## 📃 Лицензия

Этот репозиторий распространяется под лицензией [Apache 2.0 License](https://github.com/SwanHubX/SwanLab/blob/main/LICENSE).

## История звезд

[![Star History Chart](https://api.star-history.com/svg?repos=swanhubx/swanlab&type=Date)](https://star-history.com/#swanhubx/swanlab&Date)

<!-- ссылки -->

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

[demo-google-stock]:https://swanlab.cn/@ZeyiLin/Google-Stock-Prediction/charts
[demo-google-stock-image]: readme_files/example-lstm.png

[demo-audio-classification]:https://swanlab.cn/@ZeyiLin/PyTorch_Audio_Classification/charts
[demo-audio-classification-image]: readme_files/example-audio-classification.png

[demo-qwen2-vl]:https://swanlab.cn/@ZeyiLin/Qwen2-VL-finetune/runs/pkgest5xhdn3ukpdy6kv5/chart
[demo-qwen2-vl-image]: readme_files/example-qwen2-vl.jpg

[tracking-swanlab-shield-link]:https://swanlab.cn
[tracking-swanlab-shield]: https://raw.githubusercontent.com/SwanHubX/assets/main/badge2.svg

[visualize-swanlab-shield-link]:https://swanlab.cn
[visualize-swanlab-shield]: https://raw.githubusercontent.com/SwanHubX/assets/main/badge1.svg

[dockerhub-shield]: https://img.shields.io/docker/v/swanlab/swanlab-next?color=369eff&label=docker&labelColor=black&logoColor=white&style=flat-square
[dockerhub-link]: https://hub.docker.com/r/swanlab/swanlab-next/tags
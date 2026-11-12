<div align="center">

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="readme_files/swanlab-logo-type2-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="readme_files/swanlab-logo-type2-light.svg">
  <img alt="SwanLab" src="readme_files/swanlab-logo-type2-light.svg" width="300" height="130">
</picture>

プロフェッショナルかつモダンな AI 学習分析プラットフォーム  
モデル開発チーム向けに設計 — 50 以上の主要フレームワークと連携し、実験コードにシームレスに統合

<a href="https://swanlab.cn">🔥SwanLab クラウド版</a> · <a href="https://docs.swanlab.cn/en/">📃 ドキュメント</a> · <a href="https://github.com/swanhubx/swanlab/issues">Issue 報告</a> · <a href="https://rcnpx636fedp.feishu.cn/share/base/form/shrcnwT9dTedqUchdpgzPzKvcmd">フィードバック</a> · <a href="https://docs.swanlab.cn/en/guide_cloud/general/changelog.html">更新履歴</a> · <img height="16" width="16" src="https://raw.githubusercontent.com/SwanHubX/assets/main/community.svg" alt="swanlab community Logo" /> <a href="https://swanlab.cn/benchmarks">ベンチマークコミュニティ</a> · <a href="https://github.com/SwanHubX/SwanLab-Skill">SKILL</a>

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

[中文](README.md) / [English](README_EN.md) / 日本語 / [Русский](README_RU.md)

👋 [WeChat コミュニティ](https://docs.swanlab.cn/en/guide_cloud/community/online-support.html)に参加する

<a href="https://hellogithub.com/repository/b442a9fa270e4ccb8847c9ee3445e41b" target="_blank"><img src="https://abroad.hellogithub.com/v1/widgets/recommend.svg?rid=b442a9fa270e4ccb8847c9ee3445e41b&claim_uid=Oh5UaGjfrblg0yZ" alt="Featured｜HelloGitHub" style="width: 250px; height: 54px;" width="250" height="54" /></a>

</div>

<br/>

## 目次

- [👋🏻 SwanLab とは](#-swanlab-とは)
- [📃 オンラインデモ](#-オンラインデモ)
- [🏁 クイックスタート](#-クイックスタート)
- [💻 セルフホスティング](#-セルフホスティング)
- [🔥 実践事例](#-実践事例)
- [🎮 ハードウェア監視](#-ハードウェア監視)
- [🚗 フレームワーク統合](#-フレームワーク統合)
- [🔌 プラグインと API](#-プラグインと-api)
- [🆚 既存ツールとの比較](#-既存ツールとの比較)
- [👥 コミュニティ](#-コミュニティ)
- [📃 ライセンス](#-ライセンス)

<br/>

## 👋🏻 SwanLab とは

SwanLab は、モデル開発チーム向けの **AI 学習分析・メトリクス可観測性（Observability）プラットフォーム**です。**学習プロセスの可視化、自動ロギング、ハイパーパラメータ管理、実験比較、チームコラボレーション**などの機能を備え、学習中の課題を素早く特定し、モデルのイテレーションを加速させます。

直感的なインタラクティブチャートで学習のボトルネックを発見し、複数の実験を比較して改善のヒントを得られるほか、**Web ブラウザでの共有**や組織単位での**共同実験管理**によってチーム内の情報共有を円滑にし、研究開発の生産性を高めます。

https://github.com/user-attachments/assets/7965fec4-c8b0-4956-803d-dbf177b44f54

主な機能は以下の通りです：

**1. 📊 実験メトリクスとハイパーパラメータのトラッキング**: 最小限のコード変更で機械学習パイプラインに組み込み、重要な学習指標を記録します。

- ☁️ **クラウド対応**（Weights & Biases と同様）：いつでもどこからでも学習の進捗を確認可能。[モバイル端末での確認方法](https://docs.swanlab.cn/en/guide_cloud/general/app.html)

- 📝 **ハイパーパラメータの記録**、**メトリクスサマリー**、**テーブル分析**をサポート。

- 🌸 **学習プロセスの可視化**: トラッキングデータを洗練された UI で可視化することで、各ステップの挙動やメトリクスの推移を直感的に把握でき、どの変更がモデルの精度向上に寄与したかを迅速に分析できます。

- **対応メタデータ型**: スカラー指標、画像、音声、テキスト、動画、3D 点群（Point Cloud）、分子構造、ECharts カスタムチャートなど...

![swanlab-table](readme_files/molecule.gif)

- **対応チャート種別**: 折れ線グラフ、メディアパネル（画像・音声・テキスト・動画）、3D 点群、分子ビューア、棒グラフ、散布図、箱ひげ図、ヒートマップ、円グラフ、レーダーチャート、[カスタムチャート](https://docs.swanlab.cn/en/guide_cloud/experiment_track/log-custom-chart.html)など...

[![swanlab-echarts](readme_files/echarts.png)](https://swanlab.cn/@ZeyiLin/swanlab-echarts-demo/charts)

- **LLM 生成テキスト可視化コンポーネント**: 大規模言語モデル（LLM）の学習向けに設計されたテキスト可視化チャートで、Markdown レンダリングにネイティブ対応しています。

![text-chart](readme_files/text-chart.gif)

- **バックグラウンド自動ロギング**: コンソールログ（`stdout`/`stderr`）、ハードウェア環境、Git リポジトリ情報、Python 実行環境、インストール済みパッケージ一覧、プロジェクト実行ディレクトリを自動収集。

- **チェックポイントからの再開記録（Resume Logging）**: 学習の完了後や中断後に、同一の実験ランへ新しいメトリクスを追記可能。

**2. ⚡️ 幅広いフレームワーク統合**: PyTorch、🤗HuggingFace Transformers、PyTorch Lightning、🦙LLaMA Factory、MMDetection、Ultralytics、PaddleDetection、LightGBM、XGBoost、Keras、TensorBoard、Weights & Biases、MS-Swift、XTuner、Stable-Baselines3、Hydra を含む **30 以上**のフレームワークと連携。

![](readme_files/integrations.png)

**3. 💻 ハードウェア・アクセラレータ監視**: CPU、NPU（**Huawei Ascend**）、GPU（**NVIDIA**）、AMD（**AMD ROCm**）、MLU（**Cambricon**）、XPU（**Kunlunxin**）、DCU（**Hygon**）、MetaX GPU（**MetaX**）、Moore Threads GPU（**Moore Threads**）、Iluvatar GPU（**Iluvatar CoreX**）、システムメモリのハードウェア指標をリアルタイムに記録・監視。

**4. 📦 実験管理**: 機械学習ワークロード向けに設計された集中型ダッシュボードにより、複数プロジェクトや実験の全体像を一目で把握・管理。

**5. 🆚 実験結果の比較**: インタラクティブなテーブルや比較チャートを通じて、複数実験のハイパーパラメータや評価指標を横断的に比較し、次のイテレーションの洞察を獲得。

![](readme_files/swanlab-table.png)

**6. 👥 オンラインコラボレーション**: 同一プロジェクト下でチームの実験をリアルタイムに同期し、メンバーの学習ログを閲覧しながら結果に基づいたフィードバックや議論が可能。

**7. ✉️ 簡単な結果共有**: 実験ごとに発行される永続 URL をコピーしてチームメンバーに共有したり、Notion や Feishu などのオンラインドキュメントへチャートを埋め込んだりできます。

**8. 💻 セルフホスティング対応**: オフライン（エアギャップ）環境やオンプレミス環境での運用に対応。セルフホスト可能なコミュニティ版でもダッシュボードの閲覧と実験管理が可能です（詳細は[セルフホスティング](#-セルフホスティング)を参照）。

**9. 🔌 プラグイン拡張**: [Lark 通知](https://docs.swanlab.cn/en/plugin/notification-lark.html)、[Slack 通知](https://docs.swanlab.cn/en/plugin/notification-slack.html)、[CSV ロガー](https://docs.swanlab.cn/en/plugin/writer-csv.html)などのプラグインによって SwanLab のワークフローを柔軟に拡張できます。

> \[!IMPORTANT]
>
> **リポジトリを Star** していただくと、GitHub から最新リリースの通知をリアルタイムで受け取れます！ ⭐️

![star-us](readme_files/star-us.png)

<br>

## 📃 オンラインデモ

SwanLab のインタラクティブデモをお試しください：

|               [ResNet50 犬・猫画像分類][demo-cats-dogs]                |                     [YOLOv8-COCO128 物体検出][demo-yolo]                     |
| :--------------------------------------------------------------------: | :--------------------------------------------------------------------------: |
|              [![][demo-cats-dogs-image]][demo-cats-dogs]               |                      [![][demo-yolo-image]][demo-yolo]                       |
| 犬・猫データセットを用いた ResNet50 画像分類モデルの学習トラッキング。 | YOLOv8 と COCO128 データセットによる物体検出のハイパーパラメータ・指標記録。 |

|           [Qwen2 指示ファインチューニング (SFT)][demo-qwen2-sft]            |                  [LSTM Google 株価予測][demo-google-stock]                  |
| :-------------------------------------------------------------------------: | :-------------------------------------------------------------------------: |
|                 [![][demo-qwen2-sft-image]][demo-qwen2-sft]                 |              [![][demo-google-stock-image]][demo-google-stock]              |
| Qwen2 大規模言語モデルの Instruction Tuning（指示追従）学習のトラッキング。 | Google の過去株価データセットで LSTM モデルを学習し、将来の株価推移を予測。 |

|           [ResNeXt101 音声分類][demo-audio-classification]            |               [Qwen2-VL COCO データセット微調整][demo-qwen2-vl]                |
| :-------------------------------------------------------------------: | :----------------------------------------------------------------------------: |
|   [![][demo-audio-classification-image]][demo-audio-classification]   |                   [![][demo-qwen2-vl-image]][demo-qwen2-vl]                    |
| 音声分類タスクにおける ResNet から ResNeXt への段階的な実験プロセス。 | マルチモーダルモデル Qwen2-VL に対する COCO2014 データセットでの LoRA 微調整。 |

|        [EasyR1 マルチモーダル LLM 強化学習][demo-easyr1-rl]        |         [Qwen2.5-0.5B GRPO 強化学習][demo-qwen2-grpo]          |
| :----------------------------------------------------------------: | :------------------------------------------------------------: |
|            [![][demo-easyr1-rl-image]][demo-easyr1-rl]             |         [![][demo-qwen2-grpo-image]][demo-qwen2-grpo]          |
| EasyR1 フレームワークを用いたマルチモーダル LLM の強化学習（RL）。 | GSM8K データセット上での Qwen2.5-0.5B モデルの GRPO 強化学習。 |

[その他の事例・チュートリアル](https://docs.swanlab.cn/en/examples/mnist.html)

<br>

## 🏁 クイックスタート

### 1. インストール

```bash
pip install swanlab
```

<details><summary>ソースコードからのインストール</summary>

最新の開発版機能を試したい場合は、ソースコードからインストールできます。

```bash
# 方法 1：Editable インストール
git clone https://github.com/SwanHubX/SwanLab.git
pip install -e .

# 方法 2：Git URL から直接インストール
pip install git+https://github.com/SwanHubX/SwanLab.git
```

</details>

<details><summary>オフラインダッシュボード拡張機能のインストール</summary>

[オフラインダッシュボードのドキュメント](https://docs.swanlab.cn/en/guide_cloud/offline-board.html)

```bash
pip install 'swanlab[dashboard]'
```

</details>

### 2. ログインと API キーの取得

1. [SwanLab](https://swanlab.cn) で無料アカウントを作成します。

2. ログイン後、**ユーザー設定** > [API Key](https://swanlab.cn/settings) から API キーをコピーします。

3. ターミナルを開き、以下のコマンドを実行します：

```bash
swanlab login
```

プロンプトが表示されたら API キーを貼り付けて Enter キーを押し、認証を完了します。

### 3. 学習コードへの統合

```python
import swanlab

# 新しい SwanLab 実験ランを初期化
swanlab.init(
    project="my-first-ml",
    config={"learning-rate": 0.003},
)

# メトリクスをロギング
for i in range(10):
    swanlab.log({"loss": i, "acc": i})
```

準備完了です！[SwanLab](https://swanlab.cn) にアクセスして、最初の実験結果を確認しましょう。

<br>

## 💻 セルフホスティング

セルフホスティング対応のコミュニティ版を利用すると、オフライン環境やプライベート環境で SwanLab ダッシュボードを運用できます。

![swanlab-kubernetes](./readme_files/swanlab-kubernetes.png)

デプロイメントガイド：

- [🔗 Kubernetes デプロイメントガイド](https://docs.swanlab.cn/en/self_host/kubernetes/deploy.html)
- [🔗 Docker デプロイメントガイド](https://docs.swanlab.cn/en/self_host/docker/deploy.html)
- [🔗 Docker から Kubernetes への移行ガイド](https://docs.swanlab.cn/en/self_host/docker/migration-docker-kubernetes.html)

<br>

## 🔥 実践事例

**SwanLab を活用した注目のオープンソースプロジェクト・チュートリアル：**

- [happy-llm](https://github.com/datawhalechina/happy-llm)：ゼロから学ぶ大規模言語モデル（LLM）の理論と実践チュートリアル ![GitHub Repo stars](https://img.shields.io/github/stars/datawhalechina/happy-llm)
- [self-llm](https://github.com/datawhalechina/self-llm)：Linux 環境におけるオープンソース LLM / MLLM の高速ファインチューニング（フルパラメータ / LoRA）およびデプロイ実践ガイド ![GitHub Repo stars](https://img.shields.io/github/stars/datawhalechina/self-llm)
- [Minimind](https://github.com/jingyaogong/minimind)：🚀🚀 わずか 2 時間でゼロから 26M パラメータの小型 GPT を学習！ ![GitHub Repo stars](https://img.shields.io/github/stars/jingyaogong/minimind)
- [unlock-deepseek](https://github.com/datawhalechina/unlock-deepseek)：DeepSeek シリーズの技術解説・拡張・再現プロジェクト ![GitHub Repo stars](https://img.shields.io/github/stars/datawhalechina/unlock-deepseek)
- [Qwen3-SmVL](https://github.com/ShaohonChen/Qwen3-SmVL): SmolVLM2 のビジョンエンコーダと Qwen3-0.6B を結合したマルチモーダル微調整 ![GitHub Repo stars](https://img.shields.io/github/stars/ShaohonChen/Qwen3-SmVL)
- [OPPO/Agent_Foundation_Models](https://github.com/OPPO-PersonalAI/Agent_Foundation_Models): マルチエージェント蒸留と Agent RL によるエンドツーエンドのエージェント基盤モデル ![GitHub Repo stars](https://img.shields.io/github/stars/OPPO-PersonalAI/Agent_Foundation_Models)
- [Tree-GRPO](https://github.com/AMAP-ML/Tree-GRPO): [ICLR 2026] LLM エージェント強化学習における木探索の応用 ![GitHub Repo stars](https://img.shields.io/github/stars/AMAP-ML/Tree-GRPO)
- [FIPO](https://github.com/qwenpilot/FIPO): ベースモデルから深い推論能力を引き出す Value-Free 強化学習アルゴリズム FIPO の公式実装 ![GitHub Repo stars](https://img.shields.io/github/stars/qwenpilot/FIPO)
- [llm-agent-rl-lab](https://github.com/KMnO4-zx/llm-agent-rl-lab): PPO、GRPO、GSPO、DAPO、OPD などを含む LLM エージェント向け RL アルゴリズムの再現と研究 ![GitHub Repo stars](https://img.shields.io/github/stars/KMnO4-zx/llm-agent-rl-lab)

**SwanLab を使用した研究論文：**

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

**チュートリアル記事：**

- [MNIST 手書き数字認識](https://docs.swanlab.cn/en/examples/mnist.html)
- [FashionMNIST 衣類画像分類](https://docs.swanlab.cn/en/examples/fashionmnist.html)
- [CIFAR-10 画像分類](https://docs.swanlab.cn/en/examples/cifar10.html)
- [ResNet 犬・猫画像分類](https://docs.swanlab.cn/en/examples/cats_dogs_classification.html)
- [YOLO 物体検出](https://docs.swanlab.cn/en/examples/yolo.html)
- [U-Net 医用画像セグメンテーション](https://docs.swanlab.cn/en/examples/unet-medical-segmentation.html)
- [音声分類](https://docs.swanlab.cn/en/examples/audio_classification.html)
- [DQN 強化学習 — CartPole（倒立振子）](https://docs.swanlab.cn/en/examples/dqn_cartpole.html)
- [LSTM Google 株価予測](https://docs.swanlab.cn/en/examples/audio_classification.html)
- [BERT テキスト分類](https://docs.swanlab.cn/en/examples/bert.html)
- [Stable Diffusion 画像生成ファインチューニング](https://docs.swanlab.cn/en/examples/stable_diffusion.html)
- [LLM 事前学習（Pre-training）](https://docs.swanlab.cn/en/examples/pretrain_llm.html)
- [GLM-4 指示ファインチューニング](https://docs.swanlab.cn/en/examples/glm4-instruct.html)
- [Qwen 下流タスク学習](https://docs.swanlab.cn/en/examples/qwen_finetune.html)
- [固有表現抽出（NER）](https://docs.swanlab.cn/en/examples/ner.html)
- [Qwen3 医療ドメインモデル微調整](https://docs.swanlab.cn/en/examples/qwen3-medical.html)
- [Qwen2-VL マルチモーダルモデル微調整実践](https://docs.swanlab.cn/en/examples/qwen_vl_coco.html)
- [GRPO 大規模モデル強化学習](https://docs.swanlab.cn/en/examples/qwen_grpo.html)
- [Qwen3-SmVL-0.6B マルチモーダルモデル学習](https://docs.swanlab.cn/en/examples/qwen3_smolvlm_muxi.html)
- [LeRobot による身体性 AI（Embodied AI）入門](https://docs.swanlab.cn/en/examples/robot/lerobot-guide.html)
- [GLM-4.5-Air-LoRA と SwanLab 可視化ロギング](https://github.com/datawhalechina/self-llm/blob/master/models/GLM-4.5-Air/03-GLM-4.5-Air-Lora%20%E5%8F%8A%20Swanlab%20%E5%8F%AF%E8%A7%86%E5%8C%96%E5%BE%AE%E8%B0%83.md)
- [RAG の構築方法：SwanLab ドキュメントアシスタントのオープンソース実装](https://docs.swanlab.cn/course/prompt_engineering_course/11-swanlab_rag/1.swanlab-rag.html)

🌟 掲載したいチュートリアルがありましたら、ぜひ PR をお送りください！

<br>

## 🎮 ハードウェア監視

SwanLab は AI 学習中に使用される**ハードウェアスペック**と**リアルタイムのリソース使用率**を自動記録します。対応状況は以下の通りです：

| ハードウェア      | スペック記録 | リソース監視 | プローブスクリプト                                                                                                                             |
| ----------------- | ------------ | ------------ | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| NVIDIA GPU        | ✅           | ✅           | [nvidia.py](https://github.com/SwanHubX/SwanLab/blob/main/swanlab/sdk/internal/probe_python/hardware_vendor/accelerator/nvidia.py)             |
| AMD GPU           | ✅           | ✅           | [amd.py](https://github.com/SwanHubX/SwanLab/blob/main/swanlab/sdk/internal/probe_python/hardware_vendor/accelerator/amd.py)                   |
| Huawei Ascend NPU | ✅           | ✅           | [ascend.py](https://github.com/SwanHubX/SwanLab/blob/main/swanlab/sdk/internal/probe_python/hardware_vendor/accelerator/huawei.py)             |
| Cambricon MLU     | ✅           | ✅           | [cambricon.py](https://github.com/SwanHubX/SwanLab/blob/main/swanlab/sdk/internal/probe_python/hardware_vendor/accelerator/cambricon.py)       |
| Kunlunxin XPU     | ✅           | ✅           | [kunlunxin.py](https://github.com/SwanHubX/SwanLab/blob/main/swanlab/sdk/internal/probe_python/hardware_vendor/accelerator/kunlunxin.py)       |
| Moore Threads GPU | ✅           | ✅           | [moorethreads.py](https://github.com/SwanHubX/SwanLab/blob/main/swanlab/sdk/internal/probe_python/hardware_vendor/accelerator/moorethreads.py) |
| MetaX GPU         | ✅           | ✅           | [metax.py](https://github.com/SwanHubX/SwanLab/blob/main/swanlab/sdk/internal/probe_python/hardware_vendor/accelerator/metax.py)               |
| Iluvatar GPU      | ✅           | ✅           | [iluvatar.py](https://github.com/SwanHubX/SwanLab/blob/main/swanlab/sdk/internal/probe_python/hardware_vendor/accelerator/iluvatar.py)         |
| Hygon DCU         | ✅           | ✅           | [hygon.py](https://github.com/SwanHubX/SwanLab/blob/main/swanlab/sdk/internal/probe_python/hardware_vendor/accelerator/hygon.py)               |
| CPU               | ✅           | ✅           | [cpu.py](https://github.com/SwanHubX/SwanLab/blob/main/swanlab/sdk/internal/probe_python/hardware_vendor/cpu.py)                               |
| システムメモリ    | ✅           | ✅           | [memory.py](https://github.com/SwanHubX/SwanLab/blob/main/swanlab/sdk/internal/probe_python/hardware_vendor/memory.py)                         |

その他のハードウェアのサポートをご希望の場合は、お気軽に Issue または PR をご提出ください！

<br>

## 🚗 フレームワーク統合

使い慣れたお気に入りのフレームワークを SwanLab と組み合わせて利用できます。  
以下は統合済みのフレームワーク一覧です。新たに対応してほしいフレームワークがあれば、[Issue](https://github.com/swanhubx/swanlab/issues) でぜひお知らせください。

**基盤フレームワーク**

- [PyTorch](https://docs.swanlab.cn/en/guide_cloud/integration/integration-pytorch.html)
- [MindSpore](https://docs.swanlab.cn/en/guide_cloud/integration/integration-ascend.html)
- [Keras](https://docs.swanlab.cn/en/guide_cloud/integration/integration-keras.html)

**LLM 学習フレームワーク**

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

**LLM 強化学習フレームワーク**

- [veRL](https://docs.swanlab.cn/en/guide_cloud/integration/integration-verl.html)
- [HuggingFace TRL](https://docs.swanlab.cn/en/guide_cloud/integration/integration-huggingface-trl.html)
- [NVIDIA NeMo RL](https://docs.swanlab.cn/en/guide_cloud/integration/integration-nvidia-nemo-rl.html)
- [EasyR1](https://docs.swanlab.cn/en/guide_cloud/integration/integration-easyr1.html)
- [AReaL](https://docs.swanlab.cn/en/guide_cloud/integration/integration-areal.html)
- [ROLL](https://docs.swanlab.cn/en/guide_cloud/integration/integration-roll.html)

**ロボティクス・身体性 AI フレームワーク**

- [RLinf](https://docs.swanlab.cn/en/guide_cloud/integration/integration-rlinf.html)

**画像・動画生成（Text-to-Image / Video）フレームワーク**

- [DiffSynth Studio](https://docs.swanlab.cn/en/guide_cloud/integration/integration-diffsynth-studio.html)

**ディープラーニングフレームワーク**

- [PyTorch Lightning](https://docs.swanlab.cn/en/guide_cloud/integration/integration-pytorch-lightning.html)
- [MMEngine](https://docs.swanlab.cn/en/guide_cloud/integration/integration-mmengine.html)
- [fastai](https://docs.swanlab.cn/en/guide_cloud/integration/integration-fastai.html)

**コンピュータビジョン**

- [Ultralytics](https://docs.swanlab.cn/en/guide_cloud/integration/integration-ultralytics.html)
- [MMDetection](https://docs.swanlab.cn/en/guide_cloud/integration/integration-mmdetection.html)
- [MMSegmentation](https://docs.swanlab.cn/en/guide_cloud/integration/integration-mmsegmentation.html)
- [PaddleDetection](https://docs.swanlab.cn/en/guide_cloud/integration/integration-paddledetection.html)
- [PaddleYOLO](https://docs.swanlab.cn/en/guide_cloud/integration/integration-paddleyolo.html)

**古典的機械学習・勾配ブースティング**

- [LightGBM](https://docs.swanlab.cn/en/guide_cloud/integration/integration-lightgbm.html)
- [XGBoost](https://docs.swanlab.cn/en/guide_cloud/integration/integration-xgboost.html)
- [CatBoost](https://docs.swanlab.cn/en/guide_cloud/integration/integration-catboost.html)

**モデル評価フレームワーク**

- [EvalScope](https://docs.swanlab.cn/en/guide_cloud/integration/integration-evalscope.html)

**強化学習フレームワーク**

- [Stable-Baselines3](https://docs.swanlab.cn/en/guide_cloud/integration/integration-sb3.html)

**その他のツール・フレームワーク：**

- [TensorBoard](https://docs.swanlab.cn/en/guide_cloud/integration/integration-tensorboard.html)
- [Weights & Biases](https://docs.swanlab.cn/en/guide_cloud/integration/integration-wandb.html)
- [MLflow](https://docs.swanlab.cn/en/guide_cloud/integration/integration-mlflow.html)
- [HuggingFace Accelerate](https://docs.swanlab.cn/en/guide_cloud/integration/integration-huggingface-accelerate.html)
- [Ray](https://docs.swanlab.cn/en/guide_cloud/integration/integration-ray.html)
- [Hydra](https://docs.swanlab.cn/en/guide_cloud/integration/integration-hydra.html)
- [OmegaConf](https://docs.swanlab.cn/en/guide_cloud/integration/integration-omegaconf.html)
- [SpecForge](https://docs.swanlab.cn/en/guide_cloud/integration/integration-specforge.html)

[その他の統合を見る](https://docs.swanlab.cn/en/guide_cloud/integration/integration-any-library.html)

<br>

## 🔌 プラグインと API

プラグインを活用して SwanLab の機能を拡張し、実験管理ワークフローをさらに効率化しましょう！

- [カスタムプラグインの作成](https://docs.swanlab.cn/en/plugin/custom-plugin.html)
- [メール通知](https://docs.swanlab.cn/en/plugin/notification-email.html)
- [Lark 通知](https://docs.swanlab.cn/en/plugin/notification-lark.html)
- [DingTalk 通知](https://docs.swanlab.cn/en/plugin/notification-dingtalk.html)
- [WeCom（企業微信）通知](https://docs.swanlab.cn/en/plugin/notification-wxwork.html)
- [Discord 通知](https://docs.swanlab.cn/en/plugin/notification-discord.html)
- [Slack 通知](https://docs.swanlab.cn/en/plugin/notification-slack.html)
- [Bark 通知](https://docs.swanlab.cn/en/plugin/notification-bark.html)
- [CSV ロガー](https://docs.swanlab.cn/en/plugin/writer-csv.html)

公開 API：

- [OpenAPI リファレンス](https://docs.swanlab.cn/en/api/py-openapi.html)

<br>

## 🆚 既存ツールとの比較

### TensorBoard vs SwanLab

- **☁️ クラウド同期とリモートアクセス**：
  SwanLab は学習実験をクラウド上で自動同期・保存できるため、リモートからの進捗確認、過去プロジェクトの管理、実験リンクの共有、リアルタイム通知の送信、マルチデバイスでの閲覧が簡単に行えます。一方、TensorBoard はローカル環境向けのオフライントラッキングツールです。

- **👥 チームコラボレーション**：
  複数人やチーム横断での機械学習プロジェクトにおいて、SwanLab は共同ワークスペースでのプロジェクト管理、実験リンクの共有、結果のディスカッションをスムーズに行えます。TensorBoard は主に個人利用を前提としており、複数人での共有や共同作業には不向きです。

- **💻 永続的かつ一元化されたダッシュボード**：
  ローカル PC、研究室のクラスタ、パブリッククラウドの GPU インスタンスなど、どこでモデルを学習しても、すべての結果が単一の集中型ダッシュボードに集約されます。TensorBoard のように各マシンから `tfevents` ファイルを手動でコピー・管理する手間がかかりません。

- **💪 スケーラブルなテーブル分析**：
  SwanLab のテーブルビューでは、数千件規模の実験結果を検索・ソート・フィルタリングして横断的に比較でき、タスクに最適なモデルを素早く見つけ出せます。TensorBoard は大規模な実験プロジェクトの管理には適していません。

### Weights & Biases vs SwanLab

- **Weights & Biases** はインターネット接続を前提としたプロプライエタリ（クローズドソース）な MLOps プラットフォームです。

- **SwanLab** はクラウド利用だけでなく、**オープンソースかつ無料で利用できるセルフホスティング版**も提供しており、オンプレミスやオフライン環境でも導入可能です。

<br>

## 👥 コミュニティ

### 関連リポジトリ

- [self-hosted](https://github.com/swanhubx/self-hosted)：セルフホスト（オンプレミス）デプロイ用スクリプトリポジトリ
- [SwanLab-Docs](https://github.com/swanhubx/swanlab-docs)：公式ドキュメントリポジトリ
- [SwanLab-Dashboard](https://github.com/swanhubx/swanlab-dashboard)：`swanlab watch` で起動する軽量オフラインダッシュボードの Web フロントエンドリポジトリ

### コミュニティとサポート

- [GitHub Issues](https://github.com/SwanHubX/SwanLab/issues)：バグ報告や機能リクエスト
- [メールサポート](mailto:contact@swanlab.cn)：SwanLab の利用に関するお問い合わせ・フィードバック
- <a href="https://docs.swanlab.cn/en/guide_cloud/community/online-support.html">WeChat コミュニティ</a>：SwanLab の使い方に関する質問や最新 AI 技術のディスカッション

### SwanLab README バッジ

プロジェクトで SwanLab をご利用いただいている場合は、ぜひ README に SwanLab バッジを追加してください：

[![][tracking-swanlab-shield]][tracking-swanlab-shield-link] [![][visualize-swanlab-shield]][visualize-swanlab-shield-link]

```markdown
[![](https://raw.githubusercontent.com/SwanHubX/assets/main/badge2.svg)](your experiment url)
[![](https://raw.githubusercontent.com/SwanHubX/assets/main/badge1.svg)](your experiment url)
```

その他のロゴ・デザインアセット：[assets](https://github.com/SwanHubX/assets)

### 論文での SwanLab の引用

研究において SwanLab が役立ちましたら、以下の BibTeX 形式での引用をご検討ください：

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

### SwanLab へのコントリビュート

SwanLab の開発に参加してみませんか？まずは [コントリビューションガイド](CONTRIBUTING.md) をご覧ください。

また、SNS や技術イベント、カンファレンス等で SwanLab を紹介していただけることも大変励みになります。ご支援ありがとうございます！

<br>

**Contributors**

<a href="https://github.com/swanhubx/swanlab/graphs/contributors">
  <img src="https://contrib.rocks/image?repo=swanhubx/swanlab" />
</a>

<br>

<img src="./readme_files/swanlab-and-user.png" width="50%" />

## 📃 ライセンス

本リポジトリは [Apache 2.0 License](https://github.com/SwanHubX/SwanLab/blob/main/LICENSE) のもとで公開されています。

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
[swanlab-cloud-shield]: https://img.shields.io/badge/Product-SwanLabクラウド版-636a3f?labelColor=black&style=flat-square
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

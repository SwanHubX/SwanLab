<div align="center">

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="readme_files/swanlab-logo-single-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="readme_files/swanlab-logo-single.svg">
  <img alt="SwanLab" src="readme_files/swanlab-logo-single.svg" width="70" height="70">
</picture>


<h1>SwanLab</h1>

オープンソースでモダンなデザインのディープラーニングトレーニング追跡・可視化ツール  
クラウド/オフライン使用に対応し、30以上の主要フレームワークと互換性があり、実験コードと簡単に統合可能

<a href="https://swanlab.cn">🔥SwanLab オンライン版</a> · <a href="https://docs.swanlab.cn">📃 ドキュメント</a> · <a href="https://github.com/swanhubx/swanlab/issues">問題を報告</a> · <a href="https://geektechstudio.feishu.cn/share/base/form/shrcnyBlK8OMD0eweoFcc2SvWKc">フィードバックを提案</a> · <a href="https://docs.swanlab.cn/zh/guide_cloud/general/changelog.html">更新履歴</a>

[![][release-shield]][release-link]
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

👋 [WeChatグループ](https://docs.swanlab.cn/zh/guide_cloud/community/online-support.html)に参加する

</div>


## 目次

- [🌟 最近の更新](#-最近の更新)
- [👋🏻 SwanLabとは](#-swanlabとは)
- [📃 オンラインデモ](#-オンラインデモ)
- [🏁 クイックスタート](#-クイックスタート)
- [💻 セルフホスティング](#-セルフホスティング)
- [🚗 フレームワーク統合](#-フレームワーク統合)
- [🆚 既存ツールとの比較](#-既存ツールとの比較)
- [👥 コミュニティ](#-コミュニティ)
- [📃 ライセンス](#-ライセンス)

<br/>

## 🌟 最近の更新

- 2025.03.04: MLFlowの統合を追加し、MLFlow実験をSwanLab実験に変換する機能をサポートしました。[使用ガイド](https://docs.swanlab.cn/guide_cloud/integration/integration-mlflow.html)

- 2025.03.01：新機能として、実験の移動が追加されました。

- 2025.02.24：我們與[EasyR1](https://github.com/hiyouga/EasyR1)完成了聯合集成，[使用指引](https://github.com/hiyouga/EasyR1?tab=readme-ov-file#merge-checkpoint-in-hugging-face-format)

- 2025.02.18：我們與 [Swift](https://github.com/modelscope/ms-swift) 完成了聯合集成，現在你可以在Swift的CLI/WebUI中使用SwanLab來**跟踪和可視化大模型微調實驗**，[使用指引](https://docs.swanlab.cn/guide_cloud/integration/integration-swift.html)。

- 2025.02.16：新機能として、チャートの移動グループ化とグループ作成が追加されました。

- 2025.02.09: 我們與 [veRL](https://github.com/volcengine/verl) 完成了聯合集成，現在你可以在veRL中使用SwanLab來**跟踪和可視化大模型強化學習實驗**，[使用指引](https://docs.swanlab.cn/guide_cloud/integration/integration-verl.html)。

- 2025.02.05：`swanlab.log`はネストされた辞書をサポートし、Jaxフレームワークの特性に適応 [#812](https://github.com/SwanHubX/SwanLab/pull/812)；`name`と`notes`パラメータをサポート

- 2025.01.22：`sync_tensorboardX`と`sync_tensorboard_torch`機能を追加し、この2つのTensorBoardフレームワークとの実験追跡の同期をサポート

- 2025.01.17：`sync_wandb`機能を追加し、[ドキュメント](https://docs.swanlab.cn/en/guide_cloud/integration/integration-wandb.html)、Weights & Biases実験追跡との同期をサポート；ログレンダリング性能を大幅に最適化


<details><summary>完全な更新履歴</summary>

- 2025.01.11：クラウド版はプロジェクトテーブルのパフォーマンスを大幅に最適化し、ドラッグ＆ドロップ、並べ替え、フィルタリングなどのインタラクションをサポートしました。

- 2025.01.01：折れ線グラフの**永続的スムージング**、折れ線グラフのドラッグによるサイズ変更を追加し、チャート閲覧体験を最適化

- 2024.12.22：[LLaMA Factory](https://github.com/hiyouga/LLaMA-Factory)との統合を完了し、LLaMA FactoryでSwanLabを使用して**大規模モデルのファインチューニング実験を追跡・可視化**できるようになりました。[使用ガイド](https://github.com/hiyouga/LLaMA-Factory?tab=readme-ov-file#use-swanlab-logger)

- 2024.12.15：**ハードウェア監視（0.4.0）**機能をリリースし、CPU、NPU（Ascend）、GPU（Nvidia）のシステム情報の記録と監視をサポート

- 2024.12.06：[LightGBM](https://docs.swanlab.cn/guide_cloud/integration/integration-lightgbm.html)、[XGBoost](https://docs.swanlab.cn/guide_cloud/integration/integration-xgboost.html)の統合を追加；ログ記録の1行あたりの長さ制限を引き上げ

- 2024.11.26：環境タブのハードウェアセクションで**華為昇騰NPU**と**鯤鵬CPU**の識別をサポート；クラウドプロバイダーセクションで**青雲基石智算**の識別をサポート

</details>


<br>

## 👋🏻 SwanLabとは

SwanLabは、オープンソースで軽量なAIモデルトレーニング追跡・可視化ツールで、実験の追跡、記録、比較、コラボレーションのためのプラットフォームを提供します。

SwanLabはAI研究者向けに設計され、使いやすいPython APIと美しいUIを提供し、**トレーニングの可視化、自動ログ記録、ハイパーパラメータ記録、実験比較、複数人でのコラボレーション**などの機能を提供します。SwanLabでは、研究者が直感的な可視化チャートを通じてトレーニングの問題を発見し、複数の実験を比較して研究のインスピレーションを得ることができます。また、**オンラインページ**での共有や組織ベースの**複数人でのトレーニング**を通じて、チーム間のコミュニケーションの壁を打破し、組織のトレーニング効率を向上させます。

以下はその主な機能リストです：

**1. 📊 実験指標とハイパーパラメータの追跡**: 機械学習パイプラインに簡単に組み込めるコードで、トレーニングのキー指標を追跡

- **クラウド**使用をサポート（Weights & Biasesのような）、どこからでもトレーニングの進捗を確認可能。[携帯で実験を見る方法](https://docs.swanlab.cn/guide_cloud/general/app.html)
- **ハイパーパラメータ記録**とテーブル表示をサポート
- **サポートするメタデータタイプ**: スカラー指標、画像、音声、テキスト、...
- **サポートするチャートタイプ**: 折れ線グラフ、メディアグラフ（画像、音声、テキスト）、...
- **バックグラウンドでの自動記録**: ログ記録、ハードウェア環境、Gitリポジトリ、Python環境、Pythonライブラリリスト、プロジェクト実行ディレクトリ

**2. ⚡️ 幅広いフレームワーク統合**: PyTorch、🤗HuggingFace Transformers、PyTorch Lightning、🦙LLaMA Factory、MMDetection、Ultralytics、PaddleDetetion、LightGBM、XGBoost、Keras、Tensorboard、Weights&Biases、OpenAI、Swift、XTuner、Stable Baseline3、Hydraなど**30以上**のフレームワーク

![](readme_files/integrations.png)

**3. 💻 ハードウェア監視**: CPU、NPU（昇騰Ascend）、GPU（Nvidia）、メモリのシステムレベルのハードウェア指標をリアルタイムで記録・監視

**4. 📦 実験管理**: トレーニングシーン向けに設計された集中型ダッシュボードで、全体ビューを一目で確認し、複数のプロジェクトと実験を迅速に管理

**4. 🆚 結果の比較**: オンラインテーブルと比較チャートを使用して、異なる実験のハイパーパラメータと結果を比較し、イテレーションのインスピレーションを発掘

**5. 👥 オンラインコラボレーション**: チームと協力してトレーニングを行い、実験をリアルタイムで同期させ、チームのトレーニング記録をオンラインで確認し、結果に基づいて意見や提案を共有可能

**6. ✉️ 結果の共有**: 各実験の永続的なURLをコピーして送信し、パートナーに簡単に送信したり、オンラインノートに埋め込んだり可能

**7. 💻 セルフホスティング対応**: オフライン環境での使用をサポートし、セルフホスティングのコミュニティ版でもダッシュボードの表示と実験の管理が可能

> \[!IMPORTANT]
>
> **プロジェクトをスター**すると、GitHubからすべてのリリース通知を遅延なく受け取れます～ ⭐️

![star-us](readme_files/star-us.png)

<br>

## 📃 オンラインデモ

SwanLabのオンラインデモをご覧ください：

| [ResNet50 猫犬分類][demo-cats-dogs] | [Yolov8-COCO128 物体検出][demo-yolo] |
| :--------: | :--------: |
| [![][demo-cats-dogs-image]][demo-cats-dogs] | [![][demo-yolo-image]][demo-yolo] |
| 猫犬データセットでの簡単なResNet50モデルのトレーニングを追跡。 | Yolov8を使用してCOCO128データセットで物体検出タスクを追跡。 |

| [Qwen2 指示ファインチューニング][demo-qwen2-sft] | [LSTM Google株価予測][demo-google-stock] |
| :--------: | :--------: |
| [![][demo-qwen2-sft-image]][demo-qwen2-sft] | [![][demo-google-stock-image]][demo-google-stock] |
| Qwen2大規模言語モデルの指示ファインチューニングを追跡。 | 簡単なLSTMモデルを使用してGoogle株価データセットでトレーニングし、将来の株価を予測。 |

| [ResNeXt101 音声分類][demo-audio-classification] | [Qwen2-VL COCOデータセットファインチューニング][demo-qwen2-vl] |
| :--------: | :--------: |
| [![][demo-audio-classification-image]][demo-audio-classification] | [![][demo-qwen2-vl-image]][demo-qwen2-vl] |
| ResNetからResNeXtへの音声分類タスクの進化的実験プロセス | Qwen2-VL多モーダル大規模モデルを使用してCOCO2014データセットでLoraファインチューニング。 |

[その他の例](https://docs.swanlab.cn/zh/examples/mnist.html)

<br>

## 🏁 クイックスタート

### 1.インストール

```bash
pip install swanlab
```

<details><summary>ソースコードからのインストール</summary>

最新の機能を体験したい場合は、ソースコードからインストールすることができます。

**ステップ1**: プロジェクトをクローンする

```bash
git clone https://github.com/SwanHubX/SwanLab.git
cd SwanLab
```

**ステップ2**: `swanlab/package.json` の `version` フィールドを変更します。例えば、`0.10.0` にします。

**ステップ3**: インストールする

```bash
pip install -e .
```

</details>

### 2.ログインしてAPIキーを取得

1. 無料で[アカウント登録](https://swanlab.cn)

2. アカウントにログインし、ユーザー設定 > [APIキー](https://swanlab.cn/settings)でAPIキーをコピー

3. ターミナルを開き、以下を入力：

```bash
swanlab login
```

プロンプトが表示されたら、APIキーを入力してEnterを押し、ログインを完了。

### 3.SwanLabをコードに統合

```python
import swanlab

# 新しいSwanLab実験を初期化
swanlab.init(
    project="my-first-ml",
    config={'learning-rate': 0.003},
)

# 指標を記録
for i in range(10):
    swanlab.log({"loss": i, "acc": i})
```

完了！[SwanLab](https://swanlab.cn)で最初のSwanLab実験を確認できます。

<br>

## 💻 セルフホスティング

セルフホスティングのコミュニティ版は、オフラインでSwanLabダッシュボードを表示できます。

### オフライン実験追跡

swanlab.initで`logir`と`mode`パラメータを設定することで、オフラインで実験を追跡可能：

```python
...

swanlab.init(
    logdir='./logs',
    mode='local',
)

...
```

- `mode`パラメータを`local`に設定し、実験のクラウド同期を無効化

- `logdir`パラメータはオプションで、SwanLabログファイルの保存場所を指定（デフォルトは`swanlog`フォルダ）

  - ログファイルは実験の追跡中に作成・更新され、オフラインダッシュボードの起動もこれらのログファイルに基づく

その他の部分はクラウド使用と完全に同じです。

### オフラインダッシュボードの起動

ターミナルを開き、以下のコマンドを実行してSwanLabダッシュボードを起動：

```bash
swanlab watch ./logs
```

実行が完了すると、SwanLabはローカルのURLリンクを提供します（デフォルトは[http://127.0.0.1:5092](http://127.0.0.1:5092)）

このリンクにアクセスすると、ブラウザでオフラインダッシュボードを使用して実験を確認できます。

<br>

## 🚗 フレームワーク統合

お気に入りのフレームワークをSwanLabと統合しましょう！  
以下は既に統合されているフレームワークのリストです。統合してほしいフレームワークがあれば、[Issue](https://github.com/swanhubx/swanlab/issues)を提出してフィードバックをお願いします。

**基本フレームワーク**
- [PyTorch](https://docs.swanlab.cn/guide_cloud/integration/integration-pytorch.html)
- [MindSpore](https://docs.swanlab.cn/guide_cloud/integration/integration-ascend.html)
- [Keras](https://docs.swanlab.cn/guide_cloud/integration/integration-keras.html)

**専用/ファインチューニングフレームワーク**
- [PyTorch Lightning](https://docs.swanlab.cn/guide_cloud/integration/integration-pytorch-lightning.html)
- [HuggingFace Transformers](https://docs.swanlab.cn/guide_cloud/integration/integration-huggingface-transformers.html)
- [OpenMind](https://modelers.cn/docs/zh/openmind-library/1.0.0/basic_tutorial/finetune/finetune_pt.html#%E8%AE%AD%E7%BB%83%E7%9B%91%E6%8E%A7)
- [LLaMA Factory](https://docs.swanlab.cn/guide_cloud/integration/integration-llama-factory.html)
- [Modelscope Swift](https://docs.swanlab.cn/guide_cloud/integration/integration-swift.html)
- [Sentence Transformers](https://docs.swanlab.cn/guide_cloud/integration/integration-sentence-transformers.html)
- [Torchtune](https://docs.swanlab.cn/guide_cloud/integration/integration-pytorch-torchtune.html)
- [XTuner](https://docs.swanlab.cn/guide_cloud/integration/integration-xtuner.html)
- [MMEngine](https://docs.swanlab.cn/guide_cloud/integration/integration-mmengine.html)
- [FastAI](https://docs.swanlab.cn/guide_cloud/integration/integration-fastai.html)
- [LightGBM](https://docs.swanlab.cn/guide_cloud/integration/integration-lightgbm.html)
- [XGBoost](https://docs.swanlab.cn/guide_cloud/integration/integration-xgboost.html)


**コンピュータビジョン**
- [Ultralytics](https://docs.swanlab.cn/guide_cloud/integration/integration-ultralytics.html)
- [MMDetection](https://docs.swanlab.cn/guide_cloud/integration/integration-mmdetection.html)
- [MMSegmentation](https://docs.swanlab.cn/guide_cloud/integration/integration-mmsegmentation.html)
- [PaddleDetection](https://docs.swanlab.cn/guide_cloud/integration/integration-paddledetection.html)
- [PaddleYOLO](https://docs.swanlab.cn/guide_cloud/integration/integration-paddleyolo.html)

**強化学習**
- [Stable Baseline3](https://docs.swanlab.cn/guide_cloud/integration/integration-sb3.html)
- [veRL](https://docs.swanlab.cn/guide_cloud/integration/integration-verl.html)
- [HuggingFace trl](https://docs.swanlab.cn/guide_cloud/integration/integration-huggingface-trl.html)
- [EasyR1](https://docs.swanlab.cn/guide_cloud/integration/integration-easyr1.html)

**その他のフレームワーク：**
- [Tensorboard](https://docs.swanlab.cn/guide_cloud/integration/integration-tensorboard.html)
- [Weights&Biases](https://docs.swanlab.cn/guide_cloud/integration/integration-wandb.html)
- [MLFlow](https://docs.swanlab.cn/guide_cloud/integration/integration-mlflow.html)
- [HuggingFace Accelerate](https://docs.swanlab.cn/guide_cloud/integration/integration-huggingface-accelerate.html)
- [Unsloth](https://docs.swanlab.cn/guide_cloud/integration/integration-unsloth.html)
- [Hydra](https://docs.swanlab.cn/guide_cloud/integration/integration-hydra.html)
- [Omegaconf](https://docs.swanlab.cn/guide_cloud/integration/integration-omegaconf.html)
- [OpenAI](https://docs.swanlab.cn/guide_cloud/integration/integration-openai.html)
- [ZhipuAI](https://docs.swanlab.cn/guide_cloud/integration/integration-zhipuai.html)

[その他の統合](https://docs.swanlab.cn/zh/guide_cloud/integration/integration-pytorch-lightning.html)

<br>

## 🆚 既存ツールとの比較

### Tensorboard vs SwanLab

- **☁️ オンライン使用をサポート**：
  SwanLabを使用すると、トレーニング実験をクラウド上で簡単に同期・保存でき、リモートでトレーニングの進捗を確認したり、過去のプロジェクトを管理したり、実験リンクを共有したり、リアルタイムのメッセージ通知を送信したり、複数のデバイスで実験を確認したりできます。一方、Tensorboardはオフラインの実験追跡ツールです。

- **👥 複数人でのコラボレーション**：
  複数人やチーム間での機械学習コラボレーションを行う場合、SwanLabを使用すると、複数人のトレーニングプロジェクトを簡単に管理し、実験リンクを共有し、異なるスペースで議論や意見交換ができます。一方、Tensorboardは主に個人向けに設計されており、複数人でのコラボレーションや実験の共有が難しいです。

- **💻 永続的で集中型のダッシュボード**：
  どこでモデルをトレーニングしても、ローカルマシン、ラボのクラスター、またはパブリッククラウドのGPUインスタンスであっても、結果は同じ集中型ダッシュボードに記録されます。一方、TensorBoardを使用する場合、異なるマシンからTFEventファイルをコピーして管理するのに時間がかかります。

- **💪 より強力なテーブル**：
  SwanLabのテーブルを使用すると、異なる実験からの結果を表示、検索、フィルタリングでき、数千のモデルバージョンを簡単に確認し、さまざまなタスクに最適なパフォーマンスのモデルを見つけることができます。
  TensorBoardは大規模なプロジェクトには適していません。

### Weights and Biases vs SwanLab

- Weights and Biasesは、インターネット接続が必須のクローズドソースのMLOpsプラットフォームです。

- SwanLabは、インターネット接続をサポートするだけでなく、オープンソースで無料のセルフホスティング版も提供しています。

<br>

## 👥 コミュニティ

### コミュニティとサポート

- [GitHub Issues](https://github.com/SwanHubX/SwanLab/issues)：SwanLab使用中に発生したエラーや問題
- [メールサポート](zeyi.lin@swanhub.co)：SwanLabの使用に関する問題のフィードバック
- <a href="https://docs.swanlab.cn/guide_cloud/community/online-support.html">WeChatグループ</a>：SwanLabの使用に関する問題の議論や最新のAI技術の共有

### SwanLab READMEバッジ

SwanLabを仕事で使用している場合、SwanLabバッジをREADMEに追加してください：

[![][tracking-swanlab-shield]][tracking-swanlab-shield-link]、[![][visualize-swanlab-shield]][visualize-swanlab-shield-link]

```
[![](https://raw.githubusercontent.com/SwanHubX/assets/main/badge2.svg)](your experiment url)
[![](https://raw.githubusercontent.com/SwanHubX/assets/main/badge1.svg)](your experiment url)
```

さらにデザイン素材：[assets](https://github.com/SwanHubX/assets)

### 論文でSwanLabを引用

SwanLabがあなたの研究の旅に役立った場合は、以下の形式で引用を検討してください：

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

### SwanLabへの貢献

SwanLabに貢献したいですか？まず、[貢献ガイド](CONTRIBUTING.md)をお読みください。

また、ソーシャルメディア、イベント、カンファレンスでの共有を通じてSwanLabをサポートすることも大歓迎です。心から感謝します！

<br>

**貢献者**

<a href="https://github.com/swanhubx/swanlab/graphs/contributors">
  <img src="https://contrib.rocks/image?repo=swanhubx/swanlab" />
</a>

<br>

## 📃 ライセンス

このリポジトリは [Apache 2.0 License](https://github.com/SwanHubX/SwanLab/blob/main/LICENSE) オープンソースライセンスに従います

## Star History

[![Star History Chart](https://api.star-history.com/svg?repos=swanhubx/swanlab&type=Date)](https://star-history.com/#swanhubx/swanlab&Date)

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
# 为 SwanLab 作出贡献

有兴趣为 SwanLab 做出贡献吗？我们欢迎社区的贡献！本指南讨论`swanlab`的开发工作流和内部结构。

## 📦 目录

- [标准开发流程](#标准开发流程)
- [本地调试](#本地调试)
    - [IDE 与插件](#IDE与插件)
    - [配置 Python 环境](#配置python环境)
    - [调试脚本](#调试脚本)
- [本地测试](#本地测试)
    - [python 脚本调试](#python-脚本调试)
    - [单元测试](#单元测试)

## 标准开发流程

1. 浏览 GitHub 上的[Issues](https://github.com/SwanHubX/SwanLab/issues)，查看你愿意添加的功能或修复的错误，以及它们是否已被
   Pull Request。

    - 如果没有，请创建一个[新 Issue](https://github.com/SwanHubX/SwanLab/issues/new/choose)——这将帮助项目跟踪功能请求和错误报告，并确保不重复工作。

2. 如果你是第一次为开源项目贡献代码，请转到 [本项目首页](https://github.com/SwanHubX/SwanLab) 并单击右上角的"Fork"
   按钮。这将创建你用于开发的仓库的个人副本。

    - 将 Fork 的项目克隆到你的计算机，并添加指向`swanlab`项目的远程链接：

   ```bash
   git clone https://github.com/<your-username>/swanlab.git
   cd swanlab
   git remote add upstream https://github.com/swanhubx/swanlab.git
   ```

3. 开发你的贡献

    - 确保您的 Fork 与主存储库同步：

   ```bash
   git checkout main
   git pull upstream main
   ```

    - 创建一个`git`分支，您将在其中发展您的贡献。为分支使用合理的名称，例如：

   ```bash
   git checkout -b <username>/<short-dash-seperated-feature-description>
   ```

    - 当你取得进展时，在本地提交你的改动，例如：

   ```bash
   git add changed-file.py tests/test-changed-file.py
   git commit -m "feat(integrations): Add integration with the `awesomepyml` library"
   ```

4. 发起贡献：

    - [Github Pull Request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests)
    - 当您的贡献准备就绪后，将您的分支推送到 GitHub：

   ```bash
   git push origin <username>/<short-dash-seperated-feature-description>
   ```

    - 分支上传后， `GitHub`将打印一个 URL，用于将您的贡献作为拉取请求提交。在浏览器中打开该 URL，为您的拉取请求编写信息丰富的标题和详细描述，然后提交。

    - 请将相关 Issue（现有 Issue 或您创建的 Issue）链接到您的 PR。请参阅 PR 页面的右栏。或者，在 PR
      描述中提及“修复问题链接” - GitHub 将自动进行链接。

    - 我们将审查您的贡献并提供反馈。要合并审阅者建议的更改，请将编辑提交到您的分支，然后再次推送到分支（无需重新创建拉取请求，它将自动跟踪对分支的修改），例如：

   ```bash
   git add tests/test-changed-file.py
   git commit -m "test(sdk): Add a test case to address reviewer feedback"
   git push origin <username>/<short-dash-seperated-feature-description>
   ```

    - 一旦您的拉取请求被审阅者批准，它将被合并到存储库的主分支中。

## 本地调试

### IDE 与插件

1. **使用 VSCode 作为你的开发 IDE**

   SwanLab 仓库已经配好了[VSCode](https://code.visualstudio.com/)的环境、插件与调试脚本（位于`.vscode`
   文件夹中），使用 VSCode 开发 SwanLab 会有最好的体验。

2. **安装 VSCode 插件（可选）**

   用 VSCode 打开项目，进入 [扩展] ，在搜索框输入“@recommended”，会出现一系列推荐插件，推荐全部安装这些插件。

   ![vscode-recommend](/readme_files/contribution_images/vscode_recommend.png)

### 配置 Python 环境

SwanLab 项目环境需要`python>=3.8`的支持。

必须性的 python 依赖集中记录在项目根目录下的 `requirements.txt`。

同样在项目根目录启动终端，运行以下命令安装依赖：

```Bash
# swanlab所依赖的包
pip install -r requirements.txt
pip install -r requirements-media.txt
```

编译、开发、单元测试等工作需要使用以下命令额外安装依赖：

```Bash
# 编译、单元测试等功能需要使用的包
pip install -r requirements-dev.txt
```

### 调试脚本

1. **VSCode 调试脚本**

在 VSCode-运行和调试 中，项目配置好了一系列调试脚本：

![img](/readme_files/contribution_images/debug.png)

- **开启一个实验**：运行`test/create_experiment.py`脚本

- **运行当前文件**：使用配置好的 Python 环境运行你选中的文件

- **测试当前文件**：使用 debug 模式测试你选中的文件

- **进行所有单元测试**：运行`test/unit`中的脚本对 swanlab 基础功能进行完整单元测试

- **(跳过云)进行所有单元测试**：运行`test/unit`中的脚本对 swanlab 基础功能进行完整单元测试，但是跳过云测试

- **构建项目**：打包项目为 whl 文件（pip 安装包格式）

Ps: 如果你不想使用 VSCode 进行开发，可以前往`.vscode/launch.json`，查看每个调试项对应的命令，了解其配置。

## 本地测试

进行测试的前提是你已经安装完毕所有的所需依赖。

### python 脚本调试

在完成你的改动后，可以将你用于测试的 python 脚本放到根目录或`test`文件夹下，然后通过[VSCode 脚本](#调试脚本)中的"
运行当前文件"来运行你的 Python 测试脚本, 这样你的脚本运行将使用到已改动后的 swanlab。

### 单元测试

可以通过[VSCode 脚本](#调试脚本)或者在项目根目录下运行以下命令进行单元测试：

```Bash
export PYTHONPATH=. && pytest test/unit
```

由于 swanlab 涉及与云端的交互，而云端部分是闭源的，所以如果你是第一次贡献代码，最简单的方式是只进行本地测试。
针对这种情况，请在本地根目录下创建`.env`文件，并填写如下环境变量配置：

```dotenv
SWANLAB_RUNTIME=test-no-cloud
```

这样就可以跳过云端测试，只进行本地的部分功能测试。 如果想进行完整的测试，请在`.env`中补充如下信息：

```dotenv
SWANLAB_RUNTIME=test
SWANLAB_API_KEY=<你的API KEY>
SWANLAB_API_HOST=https://swanlab.cn/api
SWANLAB_WEB_HOST=https://swanlab.cn
```

*注意：在进行云端版测试时会在您的云端账号下生成一些无用的测试实验数据，需要手动删除*

配置完后即可运行完整测试
# 为 SwanLab 作出贡献

有兴趣为 SwanLab 做出贡献吗？我们欢迎社区的贡献！本指南讨论`swanlab`的开发工作流和内部结构。

## 📦 目录

- [标准开发流程](#标准开发流程)
- [本地调试](#本地调试)
    - [IDE与插件](#IDE与插件)
    - [配置Python环境](#配置python环境)
    - [调试脚本](#调试脚本)
- [本地测试](#本地测试)

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

### IDE与插件

1. **使用VSCode作为你的开发IDE**

   SwanLab仓库已经配好了[VSCode](https://code.visualstudio.com/)的环境、插件与调试脚本（位于`.vscode`
   文件夹中），使用VSCode开发SwanLab会有最好的体验。

2. **安装VSCode插件（可选）**

   用VSCode打开项目，进入 [扩展] ，在搜索框输入“@recommended”，会出现一系列推荐插件，推荐全部安装这些插件。

   ![vscode-recommend](/readme_files/contribution_images/vscode_recommend.png)

### 配置Python环境

SwanLab项目环境需要`python>=3.8`的支持。

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

1. **VSCode调试脚本**

在 VSCode-运行和调试 中，项目配置好了一系列调试脚本：

![img](/readme_files/contribution_images/debug.png)

- **开启一个实验**：运行`test/create_experiment.py`脚本

- **运行当前文件**：使用配置好的Python环境运行你选中的文件

- **测试当前文件**：使用debug模式测试你选中的文件

- **进行所有单元测试**：运行`test/unit`中的脚本对swanlab基础功能进行完整单元测试

- **构建项目**：打包项目为whl文件（pip安装包格式）

Ps: 如果你不想使用VSCode进行开发，可以前往`.vscode/launch.json`，查看每个调试项对应的命令，了解其配置。

## 本地测试



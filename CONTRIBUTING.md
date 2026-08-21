# 为 SwanLab 作出贡献

有兴趣为 SwanLab 做出贡献吗？我们欢迎社区的贡献！本指南讨论`swanlab`的开发工作流和内部结构。

## 📦 目录

- [标准开发流程](#标准开发流程)
- [分支模型与版本维护](#分支模型与版本维护)
  - [分支模型](#分支模型)
  - [修复策略](#修复策略)
  - [Backport](#backport)
  - [发版流程](#发版流程)
- [本地开发](#本地开发)
  - [环境配置](#环境配置)
  - [测试脚本](#测试脚本)
- [本地测试](#本地测试)
  - [单元测试](#单元测试)
  - [代码风格与类型检查](#代码风格与类型检查)

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

   - 创建一个`git`分支，您将在其中发展您的贡献。分支使用`<type>/<short-dash-separated-description>`的命名方式，例如：

   ```bash
   git checkout -b fix/env-quotes
   git checkout -b feat/define-metric
   ```

   - 当你取得进展时，在本地提交你的改动，提交信息遵循 Conventional Commits 风格（`feat`、`fix`、`chore`、`refactor`
     等作为前缀），例如：

   ```bash
   git add changed-file.py tests/test-changed-file.py
   git commit -m "feat(integrations): add integration with the `awesomepyml` library"
   ```

4. 发起贡献：

   - [Github Pull Request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests)
   - 当您的贡献准备就绪后，将您的分支推送到 GitHub：

   ```bash
   git push origin fix/env-quotes
   ```

   - 分支上传后， `GitHub`将打印一个 URL，用于将您的贡献作为拉取请求提交。在浏览器中打开该 URL，为您的拉取请求编写信息丰富的标题和详细描述，然后提交。

   - 请将相关 Issue（现有 Issue 或您创建的 Issue）链接到您的 PR。请参阅 PR 页面的右栏。或者，在 PR
     描述中提及"修复问题链接" - GitHub 将自动进行链接。

   - PR 的 base 分支如何选择（`main` 还是某个 `release/vX.Y`），参见[分支模型与版本维护](#分支模型与版本维护)。

   - 我们将审查您的贡献并提供反馈。要合并审阅者建议的更改，请将编辑提交到您的分支，然后再次推送到分支（无需重新创建拉取请求，它将自动跟踪对分支的修改），例如：

   ```bash
   git add tests/test-changed-file.py
   git commit -m "test(sdk): add a test case to address reviewer feedback"
   git push origin fix/env-quotes
   ```

   - 一旦您的拉取请求被审阅者批准，它将被合并到存储库的主分支中。

## 分支模型与版本维护

### 分支模型

| 分支                       | 用途                                    | 接收的改动             |
| -------------------------- | --------------------------------------- | ---------------------- |
| `main`                     | 活跃开发线，承载下一个大版本            | 新功能、bug 修复、重构 |
| `release/vX.Y`             | 历史版本维护分支（如 `release/v0.7`）   | 仅 bug 修复            |
| `backport/<branch>/pr-<N>` | backport 产生的临时分支（工具自动创建） | 合并后即可删除         |

> 注：表中 `release/v0.7` 仅作为分支命名示例，当前实际维护的 release 分支以仓库为准。

**版本开发思想**：`main` 只有一条，始终承载最新开发成果，通常处于两种状态之一：

- **当前活跃版本**：`main` 即线上最新的版本线，以 bug 修复与优化为主，补丁版本的 tag（`vX.Y.Z`）直接打在 `main` 上发布；
- **下一个开发中的大版本**：新功能与重构持续合入 `main`，为下一次大版本发布积累内容。

新大版本发布时，从当时的 `main` 切出对应的 `release/vX.Y` 分支承接历史版本维护，此后 `main` 进入下一版本的活跃开发。`release/vX.Y` 分支只做缺陷修复，不接受新功能，以保证历史版本的稳定。

### 修复策略

提交 PR 前，先判断 bug 影响的范围：

- **通用 bug**（`main` 与历史版本都受影响）：PR 的 base 选择 `main`，合并后再通过 [Backport](#backport)
  将修复合入需要维护的 `release/vX.Y` 分支。
- **分支特有 bug**（仅某个 `release/vX.Y` 受影响，通常由该版本特有的代码引入）：PR 的 base 直接选择该
  `release/vX.Y` 分支，无需 backport。

### Backport

将 `main` 上已合并的修复合入 `release/vX.Y` 有两种方式。

**方式一：label 自动化（推荐）**

在 PR 合并前为其添加 `backport release/vX.Y` label（如 `backport release/v0.8`，可同时添加多个）：

- PR 合并后，[Backport workflow](.github/workflows/backport.yml) 会自动 cherry-pick 并向目标分支创建标题带
  `[backport]` 前缀的 PR，审核合并即可。
- cherry-pick 发生冲突时，workflow 会在原 PR 下留言说明，此时改用方式二手动处理。

**方式二：手动执行（兜底）**

适用于忘了打 label、或自动 backport 遇到冲突的场景。使用仓库提供的 Makefile 目标：

```bash
make backport PR=<pr-number> BRANCH=<release-branch>
# 例如：make backport PR=1744 BRANCH=release/v0.8
```

脚本会校验 PR 已合并、从 `origin/<BRANCH>` 切出 `backport/<BRANCH>/pr-<N>` 分支并 cherry-pick，成功后自动推送并创建
PR；同时拦截重复 backport（目标分支已包含该 PR 的改动时报错退出）。

若 cherry-pick 出现冲突，脚本会中止并保留冲突现场，按提示手动完成后续步骤：

1. 解决冲突后执行 `git add <files> && git cherry-pick --continue`
2. `git push -u origin backport/<BRANCH>/pr-<N>`
3. 使用脚本打印的 `gh pr create` 命令创建 PR

也可以随时通过 `git cherry-pick --abort && git branch -D backport/<BRANCH>/pr-<N>` 放弃本次 backport。

> **维护者提示**：新切出 `release/vX.Y` 分支时，需要同步在仓库创建对应的 `backport release/vX.Y` label，
> label 自动化才会对该分支生效。

### 发版流程

面向维护者，以发布 `vX.Y.Z` 为例：

1. 确认目标 `release/vX.Y` 分支上要发布的修复均已合入；
2. 在该分支上打 tag 并推送：

   ```bash
   git checkout release/vX.Y
   git tag vX.Y.Z
   git push origin vX.Y.Z
   ```

3. tag 推送后，[发布 workflow](.github/workflows/publish-to-pypi.yml) 自动完成构建（构建时会以 tag 版本覆写
   `swanlab/package.json` 中的版本号）→ 发布到 PyPI → 创建 GitHub Release → 发送通知。

> 注：`release/vX.Y` 分支上 `swanlab/package.json` 中残留的 `x.y.0-dev` 版本号是切分支时的遗留，
> 实际发布版本以 tag 为准。

## 本地开发

### 环境配置

SwanLab 使用 [uv](https://docs.astral.sh/uv/) 管理依赖和测试，使用 [hatch](https://hatch.pypa.io/) 进行构建，
请确保本机安装了 `uv`、`make` 与 Python 3.9 及以上版本。

首次配置开发环境（安装全部依赖并挂载 pre-commit 钩子）：

```bash
make init
```

后续同步依赖：

```bash
make sync
```

pre-commit 钩子会在提交时运行 `ruff-check`（不自动修复）与 `ruff-format`，Python 侧不合规的提交会被直接拦截。

### 测试脚本

本地调试时，建议将临时验证脚本放在 `playground/` 目录下（该目录已被 gitignore，不会被误提交），脚本运行时使用的即是
本地改动后的 swanlab：

```bash
uv run python playground/your_script.py
```

## 本地测试

进行测试的前提是你已经安装完毕所有的所需依赖。

### 单元测试

运行全部单元测试：

```bash
make unit
# 或：uv run pytest
```

运行单个测试文件或用例：

```bash
uv run pytest tests/unit/path/to/test_file.py
uv run pytest tests/unit/path/to/test_file.py::TestClass::test_name
```

运行性能基准测试：

```bash
make bench
```

### 代码风格与类型检查

格式化代码（import 排序 + ruff format）：

```bash
make format
```

CI 会在每个 PR 上运行 `ruff check` 与 `basedpyright` 类型检查（见
[test-python-pr.yml](.github/workflows/test-python-pr.yml)），提交前建议本地先行验证：

```bash
uv run ruff check .
uv run basedpyright
```

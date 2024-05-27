# ⚡️更新日志

## v0.3.4 - 2024.5.27

**🚀新增功能**
- `swanlab.init`增加参数`mode`，支持新模式`disabled`
- 支持批量删除实验

**⚡️改进**
- 优化ultralytics集成代码

**👥集成**
- 与Stable Baseline3集成，[指引](/zh/guide_cloud/integration/integration-sb3.md)

## v0.3.3 - 2024.5.22

**👥集成**
- 与Weights & Biases集成，支持将wandb项目转换为`SwanLab`项目，[指引](/zh/guide_cloud/integration/integration-wandb.md)
- 与Ultralytics集成，[指引](/zh/guide_cloud/integration/integration-ultralytics.md)
- 与fastai集成，[指引](/zh/guide_cloud/integration/integration-fastai.md)

## v0.3.2 - 2024.5.17

**👥集成**
- 与Tensorboard集成，支持将Tensorboard日志文件转换为SwanLab实验，指引

**🚀新增功能**
- 支持下载折线图为PNG图像
- SwanLab实验可以被嵌入到在线文档中了（飞书/Notion等支持嵌入网页的在线文档）
- 表格视图支持导出CSV
- 表格视图支持仅看指标

**⚡️改进**
- 优化了折线图与表格视图的数值显示

**⚙️修复问题**
- 修复了在Windows系统下，swanlab.config载入hydra配置文件时，config表格的显示Bug
- 解决SwanLab在jupyter Notebook中的登录问题


## v0.3.1 - 2024.5.3

**⚡️改进**
- `swanlog`日志文件夹默认增加一个`.gitignore`

**⚙️修复问题**
- 修复`swanlab.init`的config不兼容Omegaconfig等类型的问题


## v0.3.0 云端版上线！ - 2024.5.1

**🚀新增功能**
- [云端版](https://dev101.swanlab.cn)发布
- `swanlab.init`支持用配置文件初始化
- “环境”增加对Apple M系列芯片的记录

**👥集成**
- 与🤗HuggingFace Transformers集成
- 与PyTorch Lightning集成
- 与Jupyter Notebook集成
- 与MMEngine集成
- 与openai集成
- 与zhipuai集成
- 与Hydra集成

**⚡️改进**
- 优化折线图在大数据量时的渲染表现
- 优化在Jupyter使用的表现
- 修复历史版本的大量问题


## 开源库 v0.2.4 - 2024.3.17

**⚡️改进**
- 改进图像图表的显示效果
- `swanlab.Image`支持传入pytorch tensor类型

**🔧修复问题**
- 修复多媒体图表captiop参数问题
- 修复在`swanlab.init`初始化的中途停止进程引发的问题
- 修复jupyter notebook初始化SwanLab会报错的问题


## 开源库 v0.2.3 - 2024.3.12

**🚀新增功能**
- 增加折线图平滑功能，支持3种平滑算法
- 增加图表Pin
- 增加图表隐藏

**⚡️改进**
- 支持记住图表分组的折叠状态
- 多媒体图表默认查看最后一个step

## 开源库 v0.2.2 - 2024.3.4

**⚡️改进**
- 多媒体图表支持通过方向键切换查看

**🔧修复问题**
- 修复部分Bug

## 开源库 v0.2.1 - 2024.3.1
**🚀新增功能**
- 文本图表：支持log文本类型
- 多实验图表支持图像图表和音频图表

**⚡️改进**
- 折线图组件改进性能

## 开源库 v0.2.0 - 2024.2.8

**新增功能**
- 多实验对比图表：支持项目下多个实验的日志数据在一张图表中对比
- 图像图表：支持log图像类型（支持文件、numpy array、PIL.Image、matplotlib）
- 音频图表：支持log音频类型（支持文件、numpt array）
- 支持查看1个Chart的信息时，自动查看其他Chart在相同位置的信息

**改进**
- swanlab.init增加suffix参数支持自定义实验后缀
- swanlab.log将loggings参数改为logger参数，将支持以字典的形式控制自动打印的内容
- 将默认实验名称格式改为：'%b%d-%h-%m-%s'(example:'Feb03_14-45-37')
- 将Environment中的logdir项变更为具体实验的日志文件路径
- 增加硬件数据监看类
- 改进大量UI细节

**修复问题**
- 修复部分由swanlab.log->step参数引起的折线图显示错误
- 修复部分情况下logs加载失败的问题

开源库 v0.1.6    2024.1.25
**新增功能**
- swanlab.init与swanlab.log增加APIloggings，开启loggings时将自动打印指标到终端
- 新的Config/Summary表格组件，支持参数搜索

**优化**
- 优化网页字体
- 优化网页Header
- 优化swanlab.log出现NaN的情况

**修复**
- 修复在python3.8版本时运行swanlab watch会报错的问题
- 修复当swanlab.log出现非兼容类型数据时会导致GridView和Summary组件崩溃的问题


## 开源库 v0.1.5 - 2024.1.23

**新增功能**
- 🚨使用SQLite数据库和Peewee库替代了之前的基础配置信息读写方案（#114）
  - 这是个极大有利于项目未来的改动，但缺陷是不兼容旧版本（swanlab<=v0.1.4）的日志数据文件，所以如需可视化旧版本产生的日志文件，请使用[转换脚本](https://github.com/SwanHubX/SwanLab/blob/main/script/transfer_logfile_0.1.4.py)
- 实验列表支持快速导出CSV文件
- 实验列表支持“仅看Summary”
- 实验列表支持搜索
- 环境项增加“快捷复制”交互
- 自动环境记录增加 logdir、Run path
- 新的APIswanlab.config

**优化**
- 优化部分UI
- 折线图增加Y轴线
- 优化报错信息：swanlab watch读取的日志文件路径不存在

**修复**
- 修复在hydra库传入参数时会报错的问题
- 修复swanlab.log传入字典的key在出现空格时报错的问题
- 修复运行训练脚本的路径下没有初始化git时会出现报错的问题


## 开源库 v0.1.4 - 2024.1.13

**新增功能**

- 全新的UI & UX
- 页面响应式：优化手机、平板等各种分辨率设备查看实验看板的体验
- 自动环境记录增加
  - Command 
  - Git Branch
  - Git Commit
  - 内存
  - Requirements-支持自动记录当前训练环境的pip环境列表
- Logs、Requirements支持搜索、复制与下载

**新API**
  - swanlab.init增加APIlogdir：支持设定日志文件的保存位置，
  - swanlab watch增加API--logdir：支持指定要读取的日志文件位置，
- swanlab.init->config支持形如Argparse的方式调用

**优化**
- 优化Charts横坐标的显示逻辑
- 优化Charts的自动刷新逻辑，提高性能
- 优化Charts、Summary对极小数、极大数的记录方式
- 增加单元测试

**问题修复**
- 修复记录Requirements时的问题
- 修复Charts记录取整的浮点数时的显示问题


## 开源库 v0.1.3 - 2024.1.7

**新增功能**

- 项目和实验支持在网页端删除
- 项目和实验支持在网页端修改名称和描述
- 增加实验时间后缀名
- 默认项目名为运行训练脚本的文件夹名

**优化**
- 完善实验表格Table组件的交互与性能
- 优化pypi展示信息
- 增加消息弹窗组件、占位符
- 完善消息弹窗组件样式
- 增加git分支名和最新commit hash的函数
- 将OverView 改名为 Project Dashboard

**问题修复**
- 修复记录终端日志Logs时出现的格式错乱问题
- 修复重名问题
- 修复实验看板端口占用报错问题
- 完善错误页面跳转
- 修复运行实验看板时的错误warning问题


## 开源库 v0.1.2 - 2024.1.1

- 提高Charts图表的数据更新频率，添加了请求状态映射
- 修复在训练完成时，项目和实验的OverView状态未更新的问题
- 修复当swanlab.log的data参数的key中带斜杠时会发生错误的问题
- 修复当swanlab.log的data参数的value为字符串等非标准类型时的显示问题
- 优化Chart的值域范围，现在值域范围由记录数据的最大值和最小值决定
- 优化实验表格组件
- 增加记录command、requirements的工具函数
- 添加 git hooks 用于执行格式化代码或提交测试
- 修复实验表格在上传带斜杠的key时出现的重复列问题


## 开源库 v0.1.1 - 2023.12.26

**功能更新**
- 增加【实验指标对比表】，对比各项实验配置与结果，找到最关键的训练优化点！
- 新增更多实验参数记录（Python解释器目录、系统硬件、Git仓库路径、Swanlab版本）

**优化**
- 更换全新的表格组件，得到更好的交互体验
- 时间显示格式调整

**API更新**
- swanlab.log增加参数step



## 开源库 v0.1.0 - 2023.12.22

**功能更新**
- 实验看板增加【Log】功能，自动记录终端打印信息（包含报错信息）

**优化**
- 默认语言改为【English】
- 左侧边栏文本“社区版”改为“版本号”
- 优化图表数据显示格式


## 开源库 v0.0.2    2024.12.16

**功能更新**
- 增加实验状态：【停止】
- 增加实验看板显示【配置】和【摘要】信息
- 增加【自动收集终端打印日志】功能，日志文件位置：swanlog/实验文件夹/console/yy-mm-dd.log，在下个版本我们计划将日志显示在实验看板中
- 增加实验看板的【错误页面】
- 增加swanlab watch新API：host，详细介绍：实验看板 

**重大更新**
- 日志文件结构调整（不再兼容先前的日志结构），从v0.0.2版本开始，默认的日志文件夹为swanlog

**优化**
- 优化图表样式：标题加粗、支持放大图表、显示横坐标
- 优化运行训练程序时的终端打印信息
- 优化图表颜色选用机制
- 优化图表默认分组名称与样式
- 优化默认项目名称
- 优化实验排序：默认倒序，最新进行的实验排在前面

"""
ISSUE: https://github.com/SwanHubX/SwanLab/issues/437

支持转换Tensorboard日志文件（即tfevent文件）

转换方式1 - python脚本内转换：
```python
from swanlab.converter improt TFBConverter

tfb_converter = TFBConverter(logdir="...")  # 这里也可以填一些project等参数
tfb_converter.run()
```

转换方式2 - 命令行：

```bash
swanlab convert tensorboard --logdir="..."
```
"""

from ._utils import find_tfevents, get_tf_events_tags_type, get_tf_events_tags_data
import os


class TFBConverter:
    def __init__(self, logdir: str, project: str = None):
        self.logdir = logdir
        self.project = project

    def run(self, depth=3):
        # 找到所有TFEvent文件, 生成一个路径字典
        path_dict = find_tfevents(self.logdir, depth=depth)
        if path_dict:
            print("Found file path dictionary:", path_dict)
        else:
            print("No TFEvent file found in the path, please check the path.")
            return

        for dir, paths in path_dict.items():
            """
            打印字典，了解结构，例如：
            Dir: a; Path:['path/to/logdir/a/events.out.tfevents.1715134571.qingdeMBP']
            Dir: b; Path:['path/to/logdir/b/events.out.tfevents.1715134202.qingdeMBP']
            """
            # print("=====>>>> Dir: {};\tPath:{}".format(dir, paths))

            for path in paths:
                filename = os.path.basename(path)
                # 获取所有的tag及其类型
                """
                {'training_loss': 'scalar', 'fake_image': 'image', 'fake audio': 'audio', 'fake text/text_summary': 'text'}
                """
                tags = get_tf_events_tags_type(path)
                # print("=====>>>> Tags: {} \n".format(tags))

                # 根据tag提取数据
                """
                step, value, wall_time
                {
                    'training_loss': [
                        (0, 0.0, '2024-05-08 10:16:11'),
                        (1, 0.019999999552965164, '2024-05-08 10:16:11'),
                        (2, 0.03999999910593033, '2024-05-08 10:16:11')
                        ]
                }
                """
                if tags:
                    import swanlab
                    from datetime import datetime

                    run = swanlab.init(
                        project=f"Tensorboard-Conversion-{datetime.now()}" if self.project is None else self.project,
                        experiment_name=f"{dir}/{filename}",
                        config={"tfevent_path": path},
                    )

                    data_by_tag = get_tf_events_tags_data(path, tags)
                    # print("=====>>>> Data by tag: {} \n".format(data_by_tag))

                    if data_by_tag:
                        # 打印并转换数据到SwanLab
                        for tag, data in data_by_tag.items():
                            # print(f"Data for {tag} ({tags[tag]}):")
                            for step, value, time in data:
                                # print(f"Step: {step}, Time: {time}, Value: {value}")
                                # 如果是标量
                                if tags[tag] == "scalar":
                                    swanlab.log({tag: value}, step=step)
                                # 如果是图片
                                elif tags[tag] == "image":
                                    swanlab.log({tag: swanlab.Image(value)}, step=step)
                                # 如果是音频
                                elif tags[tag] == "audio":
                                    swanlab.log({tag: swanlab.Audio(value[0], sample_rate=value[1])}, step=step)
                                # 如果是文本
                                elif tags[tag] == "text":
                                    swanlab.log({tag: swanlab.Text(value)}, step=step)

                    # TODO:
                    # 1. config增加一个runtime
                    # 2. config增加一个转换时间
                    # 3. ./_utils.py增加一个库判断报错信息

                    run.finish()

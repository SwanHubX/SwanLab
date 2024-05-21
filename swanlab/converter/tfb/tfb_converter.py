import os
import swanlab
from datetime import datetime
from ._utils import find_tfevents, get_tf_events_tags_type, get_tf_events_tags_data
from swanlab.log import swanlog as swl


class TFBConverter:
    def __init__(
        self,
        convert_dir: str,
        project: str = None,
        workspace: str = None,
        cloud: bool = True,
        logdir: str = None,
        **kwargs,
    ):
        self.convert_dir = convert_dir
        self.project = project
        self.workspace = workspace
        self.cloud = cloud
        self.logdir = logdir

    def run(self, depth=3):
        swl.info("Start converting TFEvent files to SwanLab format...")
        timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

        # 找到所有TFEvent文件, 生成一个路径字典
        path_dict = find_tfevents(self.convert_dir, depth=depth)
        if path_dict:
            swl.info("Found TFEvent file path dictionary.")
        else:
            swl.error(f"No TFEvent file found in {self.convert_dir}, please check the path.")
            return

        for dir, paths in path_dict.items():
            for path in paths:
                filename = os.path.basename(path)

                """
                获取所有的tag与其对应的类型, example:
                type_by_tags = {'training/loss': 'scalar', 'fake image': 'image', 'fake audio': 'audio', 'fake text/text_summary': 'text'}
                """
                type_by_tags = get_tf_events_tags_type(path)

                # 如果有tag（即该日志文件有记录指标，而非空文件）
                if type_by_tags:
                    # 初始化一个SwanLab实验
                    run = swanlab.init(
                        project=(f"Tensorboard-Conversion-{timestamp}" if self.project is None else self.project),
                        experiment_name=f"{dir}/{filename}",
                        workspace=self.workspace,
                        config={"tfevent_path": path},
                        cloud=self.cloud,
                        logdir=self.logdir,
                    )

                    """
                    根据tag提取数据, 格式为{tag: [(step, value, wall_time), ...]}, example:
                    data_by_tags = {
                        'training_loss': [
                            (0, 0.0, 1715839693),
                            (1, 0.019999999552965164, 1715839711),
                            (2, 0.03999999910593033, 1715839717)
                            ],
                        ...
                    }
                    """
                    data_by_tags = get_tf_events_tags_data(path, type_by_tags)

                    times = []
                    # 遍历数据
                    if data_by_tags:
                        # 打印并转换数据到SwanLab
                        for tag, data in data_by_tags.items():
                            for step, value, time in data:
                                times.append(time)
                                # 如果是标量
                                if type_by_tags[tag] == "scalar":
                                    swanlab.log({tag: value}, step=step)
                                # 如果是图片
                                elif type_by_tags[tag] == "image":
                                    swanlab.log({tag: swanlab.Image(value)}, step=step)
                                # 如果是音频
                                elif type_by_tags[tag] == "audio":
                                    swanlab.log({tag: swanlab.Audio(value[0], sample_rate=value[1])}, step=step)
                                # 如果是文本
                                elif type_by_tags[tag] == "text":
                                    swanlab.log({tag: swanlab.Text(value)}, step=step)
                                # TODO: 随着SwanLab的发展，支持转换更多类型

                    # 计算完整的运行时间
                    runtime = max(times) - min(times)
                    swanlab.config.update({"RunTime(s)": runtime})

                    # 结束当前实验
                    run.finish()

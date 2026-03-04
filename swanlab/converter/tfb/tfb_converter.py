import os
import swanlab
from datetime import datetime
from ._utils import find_tfevents, get_tf_events_tags_type, get_tf_events_tags_data
from rich.progress import Progress, BarColumn, TextColumn, TimeElapsedColumn, TimeRemainingColumn
from swanlab.log import swanlog as swl
from swanlab.data.porter import DataPorter
from swanlab.env import create_time
from swanlab.toolkit import ColumnInfo, ChartType
import time


SUPPORTED_TYPES = ["scalar", "image", "audio", "text"]

class TFBConverter:
    def __init__(
        self,
        convert_dir: str,
        project: str = None,
        workspace: str = None,
        mode: str = "cloud",
        logdir: str = None,
        types: str = None,
        **kwargs,
    ):
        self.convert_dir = convert_dir
        self.project = project
        self.workspace = workspace
        self.mode = mode
        self.logdir = logdir
        self.types = types
        if self.types is None:
            self.types = SUPPORTED_TYPES
        else:
            self.types = self.types.split(",")
            self.types = [t.strip().lower() for t in self.types]
            self.types = list(set(self.types))
            if not all(type in SUPPORTED_TYPES for type in self.types):
                raise ValueError(f"Unsupported types: {self.types}")

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
                        mode=self.mode,
                        logdir=self.logdir,
                    )
                    _porter = DataPorter._instance
                    _pre = [0]
                    if _porter is not None and _porter._pool is not None:
                        def _pre_cb(n): _pre[0] += n
                        _porter._pool.collector.upload_callback = _pre_cb
                    _conv_column_kids = {}
                    _conv_epoch_counters = {}

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
                    total_points = 0
                    if data_by_tags:
                        handlers = {
                            "scalar": lambda v: v,
                            "image": lambda v: swanlab.Image(v),
                            "audio": lambda v: swanlab.Audio(v[0], sample_rate=v[1]),
                            "text": lambda v: swanlab.Text(v),
                        }

                        # 计算总数据点数量用于进度条
                        total_points = sum(len(data) for tag, data in data_by_tags.items()
                                          if type_by_tags.get(tag) in self.types)

                        progress = None
                        task_id = None
                        if Progress is not None:
                            progress = Progress(
                                TextColumn("[bold blue]{task.description}"),
                                BarColumn(bar_width=40),
                                TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
                                TimeElapsedColumn(),
                                TimeRemainingColumn(),
                            )
                            progress.start()
                            task_id = progress.add_task(f"Parsing  {os.path.basename(path)}", total=total_points)

                        index = 0
                        for tag, data in data_by_tags.items():
                            tag_type = type_by_tags[tag]
                            if tag_type not in self.types:
                                continue
                            handler = handlers[tag_type]
                            index += 1
                            # Register column once per scalar tag
                            if tag_type == "scalar" and _porter is not None and tag not in _conv_column_kids:
                                kid = len(_conv_column_kids)
                                _conv_column_kids[tag] = str(kid)
                                split_key = tag.split("/")
                                sname = split_key[0] if len(split_key) > 1 and split_key[0] else None
                                _porter.trace_column(ColumnInfo(
                                    key=tag, kid=str(kid), name=tag, cls='CUSTOM',
                                    chart_type=ChartType.LINE, chart_reference='STEP',
                                    section_name=sname, section_type="PUBLIC",
                                ))
                            for step, value, t in data:
                                times.append(t)
                                if tag_type == "scalar" and _porter is not None:
                                    _conv_epoch_counters[tag] = _conv_epoch_counters.get(tag, 0) + 1
                                    _porter.trace_scalars_step(step, {tag: value}, dict(_conv_epoch_counters), create_time())
                                else:
                                    swanlab.log({tag: handler(value)}, step=step)
                                if progress is not None and task_id is not None:
                                    progress.update(task_id, advance=1)
                            if index % 5 == 0:
                                time.sleep(1)

                        if progress is not None:
                            progress.stop()

                    # 计算完整的运行时间
                    runtime = max(times) - min(times)
                    swanlab.config.update({"RunTime(s)": runtime})

                    # 结束当前实验（显示上传进度条）
                    _pool = _porter._pool if _porter is not None else None
                    if _pool is not None:
                        up = Progress(
                            TextColumn("[bold green]{task.description}"),
                            BarColumn(bar_width=40),
                            TextColumn("{task.completed}/{task.total} items"),
                            TimeElapsedColumn(),
                            TimeRemainingColumn(),
                        )
                        up.start()
                        t = up.add_task("Uploading to SwanLab", total=total_points, completed=_pre[0])
                        def _ucb(n): up.update(t, advance=n)
                        _pool.collector.upload_callback = _ucb
                        _oi = swl.info
                        swl.info = lambda *a, **k: None
                        try:
                            run.finish()
                        finally:
                            swl.info = _oi
                        up.update(t, completed=total_points)
                        up.stop()
                    else:
                        run.finish()

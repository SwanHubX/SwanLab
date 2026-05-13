import os
from datetime import datetime
from typing import Optional

from swanlab.converter.base import BaseConverter
from swanlab.converter.tfb.utils import find_tfevents, get_tf_events_tags_data, get_tf_events_tags_type

SUPPORTED_TYPES = ["scalar", "image", "audio", "text"]


class TFBConverter(BaseConverter):
    """Convert local TensorBoard TFEvent files to SwanLab runs."""

    def __init__(
        self,
        project: Optional[str] = None,
        workspace: Optional[str] = None,
        mode: str = "online",
        log_dir: Optional[str] = None,
        logdir: Optional[str] = None,
        types: Optional[str] = None,
    ):
        super().__init__(project=project, workspace=workspace, mode=mode, log_dir=log_dir, logdir=logdir)
        if types is None:
            self.types = SUPPORTED_TYPES
        else:
            parsed = [t.strip().lower() for t in types.split(",")]
            self.types = list(set(parsed))
            if not all(t in SUPPORTED_TYPES for t in self.types):
                raise ValueError(f"Unsupported types: {self.types}")

    def run(self, convert_dir: str = ".", depth: int = 3, **kwargs) -> None:
        import swanlab

        print("Start converting TFEvent files to SwanLab format...")
        timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

        path_dict = find_tfevents(convert_dir, depth=depth)
        if not path_dict:
            raise FileNotFoundError(f"No TFEvent file found in {convert_dir}")

        handlers = {
            "scalar": lambda v: v,
            "image": lambda v: swanlab.Image(v),
            "audio": lambda v: swanlab.Audio(v[0], sample_rate=v[1]),
            "text": lambda v: swanlab.Text(v),
        }

        for dir_name, paths in path_dict.items():
            for path in paths:
                filename = os.path.basename(path)
                type_by_tags = get_tf_events_tags_type(path)

                if not type_by_tags:
                    print(f"Skipping empty TFEvent file: {path}")
                    continue

                data_by_tags = get_tf_events_tags_data(path, type_by_tags)
                swanlab.login(relogin=True)
                run = swanlab.init(
                    project=(f"Tensorboard-Conversion-{timestamp}" if self.project is None else self.project),
                    name=f"{dir_name}/{filename}",
                    workspace=self.workspace,
                    config={"tfevent_path": path},
                    mode=self.mode,  # type: ignore[arg-type]
                    log_dir=self.log_dir,
                    reinit=True,
                )

                times = []
                for tag, data in data_by_tags.items():
                    tag_type = type_by_tags[tag]
                    if tag_type not in self.types:
                        continue
                    handler = handlers[tag_type]
                    for step, value, wall_time in data:
                        times.append(wall_time)
                        run.log({tag: handler(value)}, step=step)

                if times:
                    runtime = max(times) - min(times)
                    run.config.update({"RunTime(s)": runtime})

                run.finish()
                print(f"  Finished converting TFEvent file: {path}")

        print("All TFEvent files converted successfully.")

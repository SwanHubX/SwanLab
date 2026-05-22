"""
@author: cunyue
@description: swanlab local 模式
原本swanlab dashboard的设计理念已经不太符合swanlab 0.8.0 版本以后的架构设计
并且 dashboard 功能也不再维护，未来 dashboard 相关的功能会由 swanlab-core 重新设计并替代
本部分仅作为过渡方案（by vibe coding），使用 swanlab callbacker 机制强行兼容原本的设计
"""

import json
import math
import os
import re
import shutil
import socket
import sys
import tempfile
import time
from datetime import datetime, timezone
from importlib import import_module
from importlib.metadata import version
from pathlib import Path
from typing import Any, Optional, TextIO, cast

import click

from swanlab.sdk import settings
from swanlab.sdk.internal import pkg
from swanlab.sdk.internal.settings import Settings
from swanlab.sdk.protocol.callbacker import Callback


def _get_free_port(address: str = "0.0.0.0", default_port: int = 5092) -> int:
    """获取一个可用端口。默认返回 5092，如果被占用，返回一个随机可用端口。"""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        try:
            s.bind((address, default_port))
        except OSError:
            pass
        else:
            return default_port
    print(f"Port {default_port} is occupied, trying to find a free port...")
    sock = socket.socket()
    sock.bind((address, 0))
    _, port = sock.getsockname()
    sock.close()
    return port


@click.command()
@click.argument(
    "path",
    type=click.Path(
        exists=True,
        dir_okay=True,
        file_okay=False,
        resolve_path=True,
        readable=True,
    ),
    nargs=1,
    required=False,
    default=None,
)
@click.option(
    "--host",
    "-h",
    default=settings.integration.dashboard.host,
    type=str,
    nargs=1,
    help="The host of swanlab web, default by 127.0.0.1",
)
@click.option(
    "--port",
    "-p",
    default=settings.integration.dashboard.port,
    nargs=1,
    type=click.IntRange(1, 65535),
    help="The port of swanlab web, default by 5092",
)
@click.option(
    "--log-dir",
    "-l",
    default=None,
    nargs=1,
    type=click.Path(
        exists=True,
        dir_okay=True,
        file_okay=False,
        resolve_path=True,
        readable=True,
    ),
    help="Specify the folder to store Swanlog. Deprecated: use `swanlab watch <LOG PATH>` instead.",
)
@click.option(
    "--logdir",
    default=None,
    nargs=1,
    type=click.Path(
        exists=True,
        dir_okay=True,
        file_okay=False,
        resolve_path=True,
        readable=True,
    ),
    help="Deprecated: use --log-dir instead.",
    hidden=True,
)
def watch(path: str, host: str, port: int, log_dir: str, logdir: str):
    """Run this command to turn on the SwanLab dashboard service."""
    # logdir 兼容旧参数（向后兼容）
    if logdir is not None:
        click.echo("Warning: The option `--logdir` is deprecated, use `--log-dir` instead.")
        log_dir = logdir

    # log-dir 覆盖 path（向后兼容）
    if log_dir is not None:
        click.echo(
            "Warning: The option `--log-dir` is deprecated, use `swanlab watch [PATH]` to specify the path instead."
        )
        path = log_dir

    if path is not None:
        path = os.path.abspath(path)

    port = _get_free_port(default_port=port)

    try:
        # noinspection PyPackageRequirements
        from swanboard import SwanBoardRun
        from swanboard.utils import get_swanlog_dir
    except ModuleNotFoundError:
        raise
        click.echo("Please install the swanboard package: `pip install swanlab[dashboard]`")
        return sys.exit(1)
    # ----- 校验path，path如果被输入，已经由上层校验已存在，可读，是一个文件夹 -----
    if logdir is not None:
        pkg.console.warning(
            "The option `--logdir` will be deprecated in the future, "
            "you can just use `swanlab watch [PATH]` to specify the path."
        )
    # logdir 覆盖 path，接下来统一处理path而不再管logdir
    path = logdir if logdir is not None else path
    if path is not None:
        path = os.path.abspath(path)
        os.environ["SWANLAB_LOG_DIR"] = path
    # 为None时从环境变量中获取
    try:
        path = get_swanlog_dir()
        # 产品经理要求无论日志文件夹是否存在都不报错 🤡 ，给用户以开启web服务的“爽感”
        # if not os.path.exists(path):
        #     raise FileNotFoundError
    except ValueError as e:
        click.BadParameter(str(e))
        return sys.exit(3)
    except NotADirectoryError:
        click.BadParameter("SWANLAB_LOG_DIR must be a directory")
        return sys.exit(4)
    except FileNotFoundError:
        click.BadParameter(f"The log folder `{path}` was not found")
        return sys.exit(5)
    # ----- 校验host和port -----
    try:
        SwanBoardRun.is_valid_port(port)
        SwanBoardRun.is_valid_ip(host)
    except ValueError as e:
        click.BadParameter(str(e))
        return sys.exit(6)
    # ---- 启动服务 ----
    SwanBoardRun.run(
        path=path,
        host=host,
        port=port,
    )


class _ChartTypeItem:
    def __init__(self, chart_type: str, column_type: str):
        self.chart_type = chart_type
        self.column_type = column_type


class _ChartType:
    def __init__(self, chart_type: str, column_type: str):
        self.value = _ChartTypeItem(chart_type, column_type)


class _ColumnInfo:
    def __init__(self, key: str, kid: str, chart_type: _ChartType):
        self.key = key
        self.kid = kid
        self.name = key
        self.cls = "CUSTOM"
        self.chart_type = chart_type
        self.chart_reference = "STEP"
        self.section_name = "default"
        self.section_type = "CUSTOM"
        self.section_sort = 0
        self.error = None
        self.config = None


class _MetricContext:
    _SLICE_SIZE = 1000

    def __init__(self, key: str, kid: str, scalar: bool):
        self.key = key
        self.kid = kid
        self.scalar = scalar
        self.steps: set[int] = set()
        self._step_epochs: dict[int, int] = {}
        self._step_summary_values: dict[int, Optional[float]] = {}
        self._summary: dict[str, Any] = {}

    @property
    def summary(self) -> dict[str, Any]:
        return self._summary

    def add(
        self,
        step: int,
        data: Any,
        create_time: str,
        more: Optional[list[Any]] = None,
    ) -> tuple[dict[str, Any], bool, Path]:
        overwrite = step in self.steps
        if not overwrite:
            self.steps.add(step)
            self._step_epochs[step] = len(self.steps)

        metric = {
            "index": int(step),
            "data": data,
            "create_time": create_time,
        }
        if more is not None:
            metric["more"] = more

        if self.scalar:
            self._update_scalar_summary(step, data, overwrite)
        else:
            self._summary["num"] = len(self.steps)

        epoch = self._step_epochs[step]
        slice_name = f"{math.ceil(epoch / self._SLICE_SIZE) * self._SLICE_SIZE}.log"
        return metric, overwrite, Path(self.kid) / slice_name

    def _update_scalar_summary(self, step: int, value: Any, overwrite: bool) -> None:
        try:
            scalar = float(value)
        except (TypeError, ValueError):
            scalar = None

        if scalar is not None and (math.isnan(scalar) or math.isinf(scalar)):
            scalar = None

        needs_rebuild = False
        if overwrite:
            max_step = self._summary.get("max_step")
            min_step = self._summary.get("min_step")
            current_max = self._summary.get("max")
            current_min = self._summary.get("min")
            if scalar is None:
                needs_rebuild = step == max_step or step == min_step
            else:
                needs_rebuild = (
                    step == max_step
                    and current_max is not None
                    and scalar < current_max
                    or step == min_step
                    and current_min is not None
                    and scalar > current_min
                )

        self._step_summary_values[step] = scalar
        if needs_rebuild:
            self._rebuild_summary()
            return

        if scalar is not None:
            if self._summary.get("max") is None or scalar > self._summary["max"]:
                self._summary["max"] = scalar
                self._summary["max_step"] = step
            if self._summary.get("min") is None or scalar < self._summary["min"]:
                self._summary["min"] = scalar
                self._summary["min_step"] = step
        self._summary["num"] = len(self.steps)

    def _rebuild_summary(self) -> None:
        summary: dict[str, Any] = {"num": len(self.steps)}
        for step, _epoch in sorted(self._step_epochs.items(), key=lambda item: item[1]):
            value = self._step_summary_values.get(step)
            if value is None:
                continue
            if summary.get("max") is None or value > summary["max"]:
                summary["max"] = value
                summary["max_step"] = step
            if summary.get("min") is None or value < summary["min"]:
                summary["min"] = value
                summary["min_step"] = step
        self._summary = summary


class LocalCallbacker(Callback):
    _SUPPORTED_CHART_TYPES = {"line", "image", "audio", "text"}
    _COLUMN_TYPE_MAP = {
        1: _ChartType("line", "FLOAT"),
        2: _ChartType("image", "IMAGE"),
        3: _ChartType("audio", "AUDIO"),
        4: _ChartType("text", "TEXT"),
    }

    def __init__(self) -> None:
        self.exp: Any = None
        self._run_dir: Optional[Path] = None
        self._swanlog_dir: Optional[Path] = None
        self._console_dir: Optional[Path] = None
        self._media_dir: Optional[Path] = None
        self._metrics: dict[str, _MetricContext] = {}
        self._next_kid = 1
        self._log_file: Optional[TextIO] = None
        self._created_columns: set[str] = set()
        self._dashboard_ready = False

    @property
    def name(self) -> str:
        return "SwanLabLocalCallbacker"

    def on_run_initialized(self, run_dir: Path, path: str, settings: Settings, *args, **kwargs) -> None:
        _ = path

        log_dir = Path(run_dir).parent
        if not log_dir.exists():
            log_dir.mkdir(parents=True, exist_ok=True)
        os.environ["SWANLAB_LOG_DIR"] = str(log_dir)
        assert settings.project.name is not None, "Project name must be specified in settings for local mode"
        self._run_dir = Path(run_dir)
        self._swanlog_dir = self._run_dir / "logs"
        self._console_dir = self._run_dir / "console"
        self._media_dir = self._run_dir / "media"
        self._swanlog_dir.mkdir(parents=True, exist_ok=True)
        self._console_dir.mkdir(parents=True, exist_ok=True)
        self._media_dir.mkdir(parents=True, exist_ok=True)

        try:
            from swanboard.db import ExistedError, connect
            from swanboard.db.models import Experiment, Project
        except ImportError as e:
            raise ImportError("Please install swanboard to use 'local' mode: pip install 'swanlab[dashboard]'") from e

        if os.getenv("SWANLAB_SKIP_SWANBOARD_VERSION_CHECK") != "1" and version("swanboard") != "0.1.10b2":
            raise ImportError(
                "Your swanboard version does not match, please use this command to install the matching version: "
                "pip install 'swanlab[dashboard]'. "
                "To skip this check temporarily, set SWANLAB_SKIP_SWANBOARD_VERSION_CHECK=1."
            )

        get_swanlog_dir = getattr(import_module("swanboard.utils"), "get_swanlog_dir")
        connect(autocreate=True, path=get_swanlog_dir())
        Project.init(settings.project.name)

        run_id = settings.run.id or run_dir.name
        exp_name = settings.experiment.name or run_id
        description = settings.experiment.description or ""
        colors = self._normalize_colors(settings.experiment.color)

        pattern = r"-\d+$"
        count = Experiment.filter(Experiment.name.startswith(exp_name)).count()
        while True:
            try:
                self.exp = Experiment.create(
                    name=exp_name,
                    run_id=run_dir.name,
                    description=description,
                    num=1,
                    colors=colors,
                )
                break
            except ExistedError:
                count += 1
                if bool(re.search(pattern, exp_name)):
                    arr = exp_name.split("-")
                    arr[-1] = str(count)
                    exp_name = "-".join(arr)
                else:
                    exp_name = f"{exp_name}-{count}"
                time.sleep(0.2)
        self._dashboard_ready = True

    def on_scalar_flush(self, scalar_records: list[Any], *args, **kwargs) -> None:
        _ = (args, kwargs)
        for record in scalar_records:
            self._handle_metric(record, scalar=True)

    def on_media_flush(self, media_records: list[Any], *args, **kwargs) -> None:
        _ = (args, kwargs)
        for record in media_records:
            self._handle_metric(record, scalar=False)

    def on_log_flush(self, log_records: list[Any], *args, **kwargs) -> None:
        _ = (args, kwargs)
        if not log_records:
            return
        self._check_exp_status()
        console_dir = self._require_console_dir()
        log_name = f"{datetime.now().strftime('%Y-%m-%d')}.log"
        log_path = console_dir / log_name
        if self._log_file is None or Path(self._log_file.name) != log_path:
            if self._log_file is not None:
                self._log_file.close()
            self._log_file = open(log_path, "a", encoding="utf-8")
        for record in log_records:
            self._log_file.write(f"{record.line}\n")
        self._log_file.flush()

    def on_run_finished(self, state: str, error: Optional[str] = None, *args, **kwargs) -> None:
        _ = (state, args, kwargs)
        if error is not None:
            if self._console_dir is not None:
                console_dir = self._require_console_dir()
                with open(console_dir / "error.log", "a", encoding="utf-8") as f:
                    print(datetime.now(), file=f)
                    print(error, file=f)
        if self._dashboard_ready and self.exp is not None:
            self.exp.update_status(-1 if error is not None else 1)
        if self._log_file is not None:
            self._log_file.close()
            self._log_file = None

    def _handle_metric(self, record: Any, scalar: bool) -> None:
        self._check_exp_status()
        ctx = self._get_metric_context(record, scalar=scalar)
        if record.key not in self._created_columns:
            self._create_column(record, ctx)
            self._created_columns.add(record.key)

        data, more = (record.value.number, None) if scalar else self._media_data(record, ctx)
        metric, overwrite, relative_path = ctx.add(
            step=int(record.step),
            data=data,
            create_time=self._record_time(record),
            more=more,
        )
        metric_dir = self._require_swanlog_dir() / ctx.kid
        metric_dir.mkdir(parents=True, exist_ok=True)

        metric_path = self._require_swanlog_dir() / relative_path
        summary_path = metric_dir / "_summary.json"
        with open(summary_path, "w", encoding="utf-8") as f:
            f.write(json.dumps(ctx.summary, ensure_ascii=False))

        if overwrite:
            self._rewrite_metric_file(metric_path, metric)
        else:
            with open(metric_path, "a", encoding="utf-8") as f:
                f.write(json.dumps(metric, ensure_ascii=False) + "\n")

    def _get_metric_context(self, record: Any, scalar: bool) -> _MetricContext:
        if record.key in self._metrics:
            return self._metrics[record.key]
        kid = str(self._next_kid)
        self._next_kid += 1
        ctx = _MetricContext(record.key, kid, scalar=scalar)
        self._metrics[record.key] = ctx
        return ctx

    def _create_column(self, record: Any, ctx: _MetricContext) -> None:
        if not self._dashboard_ready or self.exp is None:
            return
        chart_type = self._COLUMN_TYPE_MAP.get(int(record.type))
        if chart_type is None:
            return
        if chart_type.value.chart_type not in self._SUPPORTED_CHART_TYPES:
            return

        try:
            from swanboard.db import ChartTypeError, ExistedError, add_multi_chart
            from swanboard.db.models import Chart, Display, Namespace, Source, Tag
        except ImportError as e:
            raise ImportError("Please install swanboard to use 'local' mode: pip install 'swanlab[dashboard]'") from e

        column_info = _ColumnInfo(record.key, ctx.kid, chart_type)
        chart = Chart.create(
            column_info.key,
            experiment_id=self.exp,
            type=column_info.chart_type.value.chart_type,
            reference=column_info.chart_reference.lower(),
        )
        try:
            namespace = Namespace.create(
                name=column_info.section_name,
                experiment_id=self.exp.id,
                sort=column_info.section_sort,
            )
        except ExistedError:
            namespace = Namespace.get(name=column_info.section_name, experiment_id=self.exp.id)
        chart_id = cast(Any, chart).id
        namespace_id = cast(Any, namespace).id
        Display.create(chart_id=chart_id, namespace_id=namespace_id)
        tag = Tag.create(
            experiment_id=cast(Any, self.exp).id,
            name=column_info.key,
            type=column_info.chart_type.value.chart_type,
            folder=column_info.kid,
        )
        tag_id = cast(Any, tag).id
        Source.create(tag_id=tag_id, chart_id=chart_id, error=cast(Any, None))
        try:
            add_multi_chart(tag_id=tag_id, chart_id=chart_id)
        except ChartTypeError:
            pass

    def _check_exp_status(self) -> None:
        if not self._dashboard_ready or self.exp is None:
            return
        try:
            from swanboard.db import NotExistedError
            from swanboard.db.models import Experiment
        except ImportError as e:
            raise ImportError("Please install swanboard to use 'local' mode: pip install 'swanlab[dashboard]'") from e

        try:
            exp = Experiment.get(id=self.exp.id)
        except NotExistedError:
            raise KeyboardInterrupt("The experiment has been deleted by the user")
        if exp.status != 0:
            raise KeyboardInterrupt("The experiment has been stopped by the user")

    def _require_swanlog_dir(self) -> Path:
        if self._swanlog_dir is None:
            raise RuntimeError("LocalCallbacker has not been initialized.")
        return self._swanlog_dir

    def _require_console_dir(self) -> Path:
        if self._console_dir is None:
            raise RuntimeError("LocalCallbacker has not been initialized.")
        self._console_dir.mkdir(parents=True, exist_ok=True)
        return self._console_dir

    def _require_media_dir(self) -> Path:
        if self._media_dir is None:
            raise RuntimeError("LocalCallbacker has not been initialized.")
        self._media_dir.mkdir(parents=True, exist_ok=True)
        return self._media_dir

    @staticmethod
    def _normalize_colors(color: Any) -> tuple[str, str]:
        if isinstance(color, (tuple, list)) and len(color) >= 2:
            return str(color[0]), str(color[1])
        if isinstance(color, str) and color:
            return color, color
        return "#000000", "#000000"

    @staticmethod
    def _record_time(record: Any) -> str:
        try:
            dt = record.timestamp.ToDatetime(tzinfo=timezone.utc)
            return dt.isoformat()
        except Exception:
            return datetime.now(timezone.utc).isoformat()

    def _media_data(self, record: Any, ctx: _MetricContext) -> tuple[list[str], list[Optional[dict[str, Any]]]]:
        filenames: list[str] = []
        more: list[Optional[dict[str, Any]]] = []
        media_dir = self._require_media_dir()
        medium = self._medium_name(record)
        legacy_dir = media_dir / ctx.kid
        legacy_dir.mkdir(parents=True, exist_ok=True)
        for item in record.value.items:
            source = media_dir / medium / item.filename
            if medium == "text":
                filenames.append(self._read_text_media(source, item.filename))
            else:
                filenames.append(item.filename)
                target = legacy_dir / item.filename
                self._ensure_legacy_media_link(source, target)
            more.append({"caption": item.caption} if item.caption else None)
        return filenames, more

    @staticmethod
    def _read_text_media(source: Path, fallback: str) -> str:
        try:
            return source.read_text(encoding="utf-8")
        except OSError:
            return fallback

    @staticmethod
    def _ensure_legacy_media_link(source: Path, target: Path) -> None:
        if target.exists() or target.is_symlink() or not source.exists():
            return
        target.parent.mkdir(parents=True, exist_ok=True)
        try:
            relative_source = os.path.relpath(source, start=target.parent)
            target.symlink_to(relative_source)
        except OSError:
            shutil.copy2(source, target)

    @staticmethod
    def _medium_name(record: Any) -> str:
        return {
            2: "image",
            3: "audio",
            4: "text",
            5: "video",
            6: "echarts",
            7: "object3d",
            8: "molecule",
        }.get(int(record.type), "unknown")

    @staticmethod
    def _rewrite_metric_file(metric_path: Path, metric: dict[str, Any]) -> None:
        serialized = json.dumps(metric, ensure_ascii=False) + "\n"
        metric_path.parent.mkdir(parents=True, exist_ok=True)
        try:
            f_in = open(metric_path, "r", encoding="utf-8")
        except FileNotFoundError:
            with open(metric_path, "w", encoding="utf-8") as f:
                f.write(serialized)
            return

        with (
            f_in,
            tempfile.NamedTemporaryFile(
                mode="w",
                encoding="utf-8",
                dir=metric_path.parent,
                delete=False,
            ) as tmp,
        ):
            tmp_path = Path(tmp.name)
            replaced = False
            for line in f_in:
                try:
                    existing = json.loads(line)
                except json.JSONDecodeError:
                    tmp.write(line)
                    continue
                if existing.get("index") != metric["index"]:
                    tmp.write(line)
                elif not replaced:
                    tmp.write(serialized)
                    replaced = True
            if not replaced:
                tmp.write(serialized)
        os.replace(tmp_path, metric_path)

# -*- coding: utf-8 -*-
"""
This code references: https://github.com/yuz1wan/localWandb2cloudSwanlab, thanks to @yuz1wan.

------example.py------
from swanlab.converter import WandbLocalConverter

wb_local_converter = WandbLocalConverter()
wb_local_converter.run(root_wandb_dir="WANDB_LOCAL_DIR", wandb_run_dir="WANDB_LOCAL_RUN_DIR")

------command------
swanlab convert -t wandb-local --wb-dir ./wandb --wb-run-dir run-1234567890
"""

# GC frequency: collect garbage every N records
GC_INTERVAL = 10000

import argparse
import gc
import time
import glob
import json

try:
    import orjson as _orjson
    _json_loads = _orjson.loads
except ImportError:
    _json_loads = json.loads

import os
import re
import traceback
from typing import Optional, List

import yaml
import swanlab
from swanlab.log import swanlog as swl
from swanlab.data.porter import DataPorter
from swanlab.env import create_time
from swanlab.toolkit import ColumnInfo, ChartType

from rich.progress import Progress, BarColumn, TextColumn, TimeElapsedColumn, TimeRemainingColumn
from rich.progress import ProgressColumn
from rich.text import Text


class SizeColumn(ProgressColumn):
    """自定义列：显示文件大小（如 84.3M/7.61G）"""

    def __init__(self):
        super().__init__()

    def render(self, task: "Task") -> "Text":  # type: ignore
        completed = task.completed
        total = task.total
        if total is None:
            return Text("-/-", style="progress.filesize")
        return Text(f"{_format_size(int(completed))}/{_format_size(int(total))}", style="progress.filesize")


def _format_size(size_bytes: int) -> str:
    """格式化字节为人类可读格式"""
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size_bytes < 1024:
            return f"{size_bytes:.1f}{unit}"
        size_bytes /= 1024
    return f"{size_bytes:.1f}PB"

try:
    from wandb.sdk.internal.datastore import DataStore
    from wandb.proto import wandb_internal_pb2
except ImportError:
    # We add a print here because if this fails, swanlog might not be available.
    print("Error: Wandb is required to parse local log files. Please install it with 'pip install wandb'.")
    exit(1)


class WandbLocalConverter:
    """
    A robust converter to migrate local wandb run data to SwanLab.

    This class reliably parses all data from a .wandb file by using the
    `ds.scan_data()` method and converting Protobuf messages to dictionaries
    with `MessageToDict`.
    """
    def __init__(
        self,
        project: Optional[str] = None,
        workspace: Optional[str] = None,
        mode: Optional[str] = "cloud",
        tags: Optional[List[str]] = None,
        logdir: Optional[str] = None,
        **kwargs,
    ):
        """
        Initializes the converter with SwanLab project details.

        Args:
            project (Optional[str]): The SwanLab project name. If None, it will be inferred from the wandb run.
            workspace (Optional[str]): The SwanLab workspace.
            logdir (Optional[str]): The local directory to save SwanLab logs.
        """
        self.project = project
        self.workspace = workspace
        self.mode = mode
        self.tags = tags
        self.logdir = logdir

    def _find_run_dirs(self, root_wandb_dir: str, wandb_run_dir: Optional[str] = None) -> list[str]:
        """Finds all wandb run directories within a given root directory."""
        if wandb_run_dir:
            patterns = [
                os.path.join(root_wandb_dir, wandb_run_dir)
            ]
            found_path = os.path.join(root_wandb_dir, wandb_run_dir)
        else:
            patterns = [
                os.path.join(root_wandb_dir, "run-*/"),
                os.path.join(root_wandb_dir, "offline-run-*/")
            ]
            found_path = root_wandb_dir

        run_dirs = []
        for pattern in patterns:
            run_dirs.extend(glob.glob(pattern))
        if not run_dirs:
            swl.warning(f"No wandb run directories found in '{found_path}'.")
        return run_dirs

    def _proto_items_to_dict(self, items) -> dict:
        """Directly read proto repeated item fields (key/nested_key/value_json) into a dict.

        This avoids the expensive MessageToDict conversion by accessing proto fields directly.
        """
        mapping = {}
        for item in items:
            key = item.key or '/'.join(item.nested_key)
            if not key:
                continue
            try:
                mapping[key] = _json_loads(item.value_json)
            except (ValueError, Exception):
                swl.warning(f"Could not decode json for key '{key}': {item.value_json}")
        return mapping

    def _parse_run(self, run_dir: str):
        """Parses a single wandb run directory and converts it to a SwanLab run."""
        wandb_files = glob.glob(os.path.join(run_dir, "*.wandb"))
        if not wandb_files:
            swl.error(f"No .wandb file found in {run_dir}. Skipping.")
            return

        wandb_file_path = wandb_files[0]
        files_root_dir = os.path.join(run_dir, "files")
        run_basename = os.path.basename(run_dir.rstrip('/\\'))

        # 获取文件大小用于进度条
        wandb_file_size = os.path.getsize(wandb_file_path)

        ds = DataStore()
        ds.open_for_scan(wandb_file_path)

        swanlab_run = None
        run_metadata = {'id': None, 'name': None, 'notes': None, 'project': None}
        run_config = {}

        # Direct-write state (bypasses swanlab_run.log() overhead for scalar metrics)
        porter = None           # set after swanlab.init()
        _conv_column_kids = {}  # key -> kid str
        _conv_epoch_counters = {}  # key -> cumulative epoch count
        _upload_pre_count = [0]  # items uploaded before upload progress bar appears

        def _log_scalars_direct(scalars, step):
            """Write float scalars directly to the porter, bypassing log() overhead."""
            if not scalars or porter is None:
                return
            for key in scalars:
                if key not in _conv_column_kids:
                    kid = len(_conv_column_kids)
                    _conv_column_kids[key] = str(kid)
                    split_key = key.split("/")
                    sname = split_key[0] if len(split_key) > 1 and split_key[0] else None
                    porter.trace_column(ColumnInfo(
                        key=key, kid=str(kid), name=key, cls='CUSTOM',
                        chart_type=ChartType.LINE, chart_reference='STEP',
                        section_name=sname, section_type="PUBLIC",
                    ))
            for key in scalars:
                _conv_epoch_counters[key] = _conv_epoch_counters.get(key, 0) + 1
            porter.trace_scalars_step(step, scalars, dict(_conv_epoch_counters), create_time())

        def _finish_with_progress():
            """Run swanlab_run.finish() while showing a Rich upload progress bar."""
            _pool = porter._pool if porter is not None else None
            if _pool is None:
                swanlab_run.finish()
                return
            total = len(_conv_column_kids) + sum(_conv_epoch_counters.values())
            up = Progress(
                TextColumn("[bold green]{task.description}"),
                BarColumn(bar_width=40),
                TextColumn("{task.completed}/{task.total} items"),
                TimeElapsedColumn(),
                TimeRemainingColumn(),
            )
            up.start()
            t = up.add_task("Uploading to SwanLab", total=total, completed=_upload_pre_count[0])
            def _upload_cb(n):
                up.update(t, advance=n)
            _pool.collector.upload_callback = _upload_cb
            # Suppress swanlog.info during finish() to hide the "View project/run" URL reprints
            _orig_info = swl.info
            swl.info = lambda *a, **k: None
            try:
                swanlab_run.finish()
            finally:
                swl.info = _orig_info
            up.update(t, completed=total)
            up.stop()

        def initialize_swanlab_run_if_needed():
            """Lazy initializer for the SwanLab run."""
            nonlocal swanlab_run, progress, task_id, porter
            if swanlab_run is not None:
                return

            # 暂停进度条，避免与 swanlab.init 中的 rich Status 冲突
            current_progress = 0
            if progress is not None and task_id is not None:
                current_progress = int(progress.tasks[task_id].completed)
                progress.stop()

            if run_metadata['id'] is None:
                match = re.search(r"-([a-zA-Z0-9]+)$", os.path.basename(run_dir.rstrip('/\\')))
                if match:
                    run_metadata['id'] = match.group(1)
                else:
                    raise RuntimeError(f"Fatal: Cannot determine run_id for {run_dir}.")
            if run_metadata['name'] is None:
                run_metadata['name'] = run_metadata['id']

            swl.info(f"Converting run: {run_metadata['name']} (ID: {run_metadata['id']})")
            swanlab_run = swanlab.init(
                project=self.project or run_metadata['project'] or "wandb-imports",
                workspace=self.workspace,
                experiment_name=run_metadata['name'],
                description=run_metadata['notes'],
                logdir=self.logdir,
                mode=self.mode,
                tags=self.tags,
            )
            porter = DataPorter._instance
            # Set pre-count callback: track items uploaded during parse before the upload bar appears
            if porter is not None and porter._pool is not None:
                def _pre_upload_cb(n):
                    _upload_pre_count[0] += n
                porter._pool.collector.upload_callback = _pre_upload_cb

            # 恢复进度条
            if progress is not None and task_id is not None:
                try:
                    progress.start()
                    progress.update(task_id, completed=current_progress)
                except Exception:
                    # 如果 start 失败，重新创建进度条并恢复进度
                    progress = Progress(
                        TextColumn("[bold blue]{task.description}"),
                        BarColumn(bar_width=40),
                        SizeColumn(),
                        TimeElapsedColumn(),
                        TimeRemainingColumn(),
                    )
                    progress.start()
                    task_id = progress.add_task(
                        f"Parsing  {run_basename}",
                        total=wandb_file_size,
                        completed=current_progress
                    )

            try:
                config_path = os.path.join(files_root_dir, "config.yaml")
                with open(config_path, 'r', encoding='utf-8') as f:
                    wandb_config_yaml = yaml.safe_load(f)
                    cleaned_config = {k: v['value'] for k, v in wandb_config_yaml.items() if isinstance(v, dict) and 'value' in v}
                    run_config.update(cleaned_config)
            except FileNotFoundError:
                pass
            
            final_config = {"wandb_run_id": run_metadata['id'], "wandb_local_path": run_dir}
            final_config.update(run_config)
            swanlab_run.config.update(final_config)

        # Core Logic: Scan the .wandb file record by record
        # 初始化进度条
        progress = None
        task_id = None
        record_count = 0
        last_summary = {}
        last_summary_step = 0

        # Performance timing stats
        timing_stats = {
            "scan_data": 0.0,
            "parse_proto": 0.0,
            "process_run": 0.0,
            "process_config": 0.0,
            "process_history": 0.0,
            "process_summary": 0.0,
            "swanlab_log": 0.0,
            "gc_collect": 0.0,
        }

        if Progress is not None:
            progress = Progress(
                TextColumn("[bold blue]{task.description}"),
                BarColumn(bar_width=40),
                SizeColumn(),
                TimeElapsedColumn(),
                TimeRemainingColumn(),
            )
            progress.start()
            task_id = progress.add_task(f"Parsing  {run_basename}", total=wandb_file_size)

        while True:
            t0 = time.perf_counter()
            record_bin = ds.scan_data()
            timing_stats["scan_data"] += time.perf_counter() - t0
            if record_bin is None:
                break

            # 更新进度条
            if progress is not None and task_id is not None:
                progress.update(task_id, advance=len(record_bin))

            t0 = time.perf_counter()
            record_pb = wandb_internal_pb2.Record()
            record_pb.ParseFromString(record_bin)
            timing_stats["parse_proto"] += time.perf_counter() - t0

            record_type = record_pb.WhichOneof("record_type")

            # Process record by directly accessing protobuf fields (no MessageToDict)
            if record_type == "run":
                t0 = time.perf_counter()
                run = record_pb.run
                run_metadata.update({
                    'id': run.run_id or None,
                    'name': run.display_name or None,
                    'notes': run.notes or None,
                    'project': run.project or None,
                })
                run_config.update(self._proto_items_to_dict(run.config.update))
                timing_stats["process_run"] += time.perf_counter() - t0
            elif record_type == "config":
                t0 = time.perf_counter()
                run_config.update(self._proto_items_to_dict(record_pb.config.update))
                timing_stats["process_config"] += time.perf_counter() - t0
            elif record_type == "history":
                t0 = time.perf_counter()
                initialize_swanlab_run_if_needed()
                scalar_dict = {}
                media_dict = {}
                step = 0
                # Single-pass: parse proto items directly, fast-path float() for scalars
                for item in record_pb.history.item:
                    key = item.key or '/'.join(item.nested_key)
                    if not key:
                        continue
                    vj = item.value_json
                    if key == '_step':
                        try:
                            step = int(float(vj))
                        except (ValueError, TypeError):
                            pass
                        continue
                    if key.startswith('_'):
                        continue
                    # Fast path: direct float conversion (avoids full JSON parse for scalars)
                    try:
                        scalar_dict[key] = float(vj)
                        continue
                    except (ValueError, TypeError):
                        pass
                    # Slow path: full JSON parse for complex types (media, etc.)
                    try:
                        value = _json_loads(vj)
                    except (ValueError, Exception):
                        continue
                    if isinstance(value, int):
                        scalar_dict[key] = float(value)
                    elif isinstance(value, dict) and "_type" in value:
                        media_type = value["_type"]
                        path = os.path.join(files_root_dir, value.get("path", ""))
                        if os.path.exists(path) and media_type == "image-file":
                            media_dict[key] = swanlab.Image(path)

                timing_stats["process_history"] += time.perf_counter() - t0

                if scalar_dict or media_dict:
                    t1 = time.perf_counter()
                    if scalar_dict:
                        _log_scalars_direct(scalar_dict, step)
                    if media_dict:
                        swanlab_run.log(media_dict, step=step)
                    timing_stats["swanlab_log"] += time.perf_counter() - t1
                    last_summary_step = step
                del scalar_dict, media_dict
            elif record_type == "summary":
                t0 = time.perf_counter()
                # Accumulate into last_summary; only log once after the loop ends.
                # wandb writes a summary record after every step, so calling log() here
                # would result in N redundant log calls (N = number of steps).
                for item in record_pb.summary.update:
                    key = item.key or '/'.join(item.nested_key)
                    if not key or key.startswith('_'):
                        continue
                    try:
                        last_summary[key] = float(item.value_json)
                    except (ValueError, TypeError):
                        pass
                timing_stats["process_summary"] += time.perf_counter() - t0

            # 清理公共变量，释放内存
            del record_pb
            del record_bin

            # GC every GC_INTERVAL records to reduce overhead
            record_count += 1
            if record_count % GC_INTERVAL == 0:
                t0 = time.perf_counter()
                gc.collect()
                timing_stats["gc_collect"] += time.perf_counter() - t0

        # 关闭进度条
        if progress is not None:
            progress.stop()

        # Log 最终的 summary 数据（只 log 一次，避免对每条 summary record 都 log）
        if last_summary:
            initialize_swanlab_run_if_needed()
            if swanlab_run:
                t0 = time.perf_counter()
                _log_scalars_direct(last_summary, last_summary_step)
                timing_stats["swanlab_log"] += time.perf_counter() - t0

        # 打印性能时间统计
        total_time = sum(timing_stats.values())
        swl.info("=" * 30)
        swl.info("Performance Timing Stats:")
        for name, secs in timing_stats.items():
            pct = (secs / total_time * 100) if total_time > 0 else 0
            swl.info(f"  {name:20s}: {secs:8.3f}s ({pct:5.1f}%)")
        swl.info(f"  {'TOTAL':20s}: {total_time:8.3f}s")
        swl.info("=" * 30)

        if swanlab_run:
            swl.info(f"Finished converting run: {run_metadata['name']}")
            _finish_with_progress()
        else:
            try:
                initialize_swanlab_run_if_needed()
                if swanlab_run:
                    swl.warning(f"Run in {run_dir} has no metrics, but its config was saved.")
                    _finish_with_progress()
            except RuntimeError as e:
                swl.error(str(e))

        # 清理 DataStore 对象，释放文件句柄
        del ds
        gc.collect()

    def run(self, root_wandb_dir: str, wandb_run_dir: Optional[str] = None):
        """The main entry point to start the conversion process."""
        swl.info(f"Starting import from wandb directory: {os.path.abspath(root_wandb_dir)}")
        run_dirs = self._find_run_dirs(root_wandb_dir, wandb_run_dir)
        
        if not run_dirs:
            swl.info("No runs found to import. Exiting.")
            return

        swl.info(f"Found {len(run_dirs)} runs to import.")
        
        for i, run_dir in enumerate(run_dirs):
            print("\n")
            swl.info(f"--- Processing run {i+1}/{len(run_dirs)} ---")
            try:
                self._parse_run(run_dir)
            except Exception as e:
                swl.error(f"FATAL: Failed to convert run in {run_dir}: {e}")
                traceback.print_exc()
        
        swl.info("="*20)
        swl.info("All runs processed.")
        swl.info("="*20)


def main():
    """Main function to parse arguments and run the converter."""
    parser = argparse.ArgumentParser(
        description="Import local Wandb runs into SwanLab.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument(
        "--wandb_dir",
        "-d",
        type=str,
        default="./wandb",
        required=True,
        help="Path to the local 'wandb' directory containing run folders (e.g., './wandb')."
    )
    
    parser.add_argument(
        "--wandb_run_dir",
        "-r",
        type=str,
        default=None,
        help="Path to the local 'wandb' directory containing run folders (e.g., './wandb')."
    )
    
    parser.add_argument(
        "--project",
        "-p",
        type=str,
        default=None,
        help="Specify the SwanLab project name. If not set, it's inferred from the wandb run."
    )
    
    parser.add_argument(
        "--workspace",
        "-w",
        type=str,
        default=None,
        help="Specify the SwanLab workspace."
    )
    
    parser.add_argument(
        "--logdir",
        "-l",
        type=str,
        default=None,
        help="Specify the root directory for SwanLab logging. Defaults to './swanlog'."
    )
    
    args = parser.parse_args()
    
    swl.info("SwanLab Importer for Wandb")
    swl.info("---------------------------")
    swl.info(f"Source Wandb Directory: {os.path.abspath(args.wandb_dir)}")
    swl.info(f"Source Wandb Run Directory: {args.wandb_run_dir or '(all runs)'}")
    swl.info(f"Target SwanLab Project: {args.project or '(inferred from run)'}")
    swl.info(f"Target SwanLab Workspace: {args.workspace or '(default)'}")
    swl.info(f"SwanLab Log Directory: {args.logdir or './swanlog'}")
    swl.info("---------------------------")

    converter = WandbLocalConverter(
        project=args.project,
        workspace=args.workspace,
        logdir=args.logdir
    )
    converter.run(root_wandb_dir=args.wandb_dir, wandb_run_dir=args.wandb_run_dir)


if __name__ == '__main__':
    main()
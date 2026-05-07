import csv
import os
import time
from typing import Optional

from swanlab.sdk.internal.pkg import console
from swanlab.sdk.protocol import Callback
from swanlab.sdk.typings.run.callback import RunInfo


class CSVWriter(Callback):
    """CSV experiment summary writer.

    Each run produces one row in CSV: metadata + config + final metric values.
    """

    def __init__(self, dir: Optional[str] = None, filename: str = "swanlab_run.csv"):
        if dir is None:
            dir = os.getcwd()
        self._dir = dir
        os.makedirs(dir, exist_ok=True)
        self._filepath = os.path.join(dir, filename)
        self._run_info: Optional[RunInfo] = None

    @property
    def name(self) -> str:
        return "CSVWriter"

    def on_run_initialized(self, run_dir, path, *args, **kwargs):
        run_info = kwargs.get("run_info")
        if isinstance(run_info, RunInfo):
            self._run_info = run_info

    def on_run_finished(self, state, error=None, **kwargs):
        if self._run_info is None:
            return

        config = kwargs.get("config", {})
        metrics = kwargs.get("metrics", {})

        row, headers = self._build_row(state, error, config, metrics)

        if not os.path.isfile(self._filepath):
            self._write_new_file(headers, row)
        else:
            self._append_to_existing(headers, row)
        console.info(f"✅ CSVWriter written to {self._filepath}")

    def _build_row(self, state, error, config, metrics):
        """构建一行数据和对应的列头。"""
        info = self._run_info
        assert info is not None
        timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

        headers = [
            "project",
            "exp_name",
            "description",
            "datetime",
            "run_id",
            "workspace",
            "logdir",
            "url",
            "state",
            "error",
        ]
        row = [
            info.project,
            info.experiment_name,
            info.description or "",
            timestamp,
            info.run_id,
            info.workspace,
            str(info.run_dir),
            info.url or "",
            state,
            error or "",
        ]

        for key in sorted(config.keys()):
            headers.append(f"config/{key}")
            row.append(config[key])

        for key in sorted(metrics.keys()):
            headers.append(key)
            row.append(metrics[key].get("latest"))

        return row, headers

    def _write_new_file(self, headers, row):
        with open(self._filepath, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(headers)
            writer.writerow(row)

    def _append_to_existing(self, new_headers, new_row):
        """读取已有 CSV，合并列，追加行。"""
        with open(self._filepath, "r", newline="") as f:
            reader = csv.reader(f)
            existing_headers = next(reader)
            existing_rows = list(reader)

        header_set = {h: i for i, h in enumerate(existing_headers)}
        for h in new_headers:
            if h not in header_set:
                header_set[h] = len(existing_headers)
                existing_headers.append(h)

        for old_row in existing_rows:
            while len(old_row) < len(existing_headers):
                old_row.append("")

        aligned_row = [""] * len(existing_headers)
        for h, val in zip(new_headers, new_row):
            aligned_row[header_set[h]] = val

        tmp_path = self._filepath + ".tmp"
        with open(tmp_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(existing_headers)
            for r in existing_rows:
                writer.writerow(r)
            writer.writerow(aligned_row)
        os.replace(tmp_path, self._filepath)

    def __str__(self):
        return "CSVWriter"

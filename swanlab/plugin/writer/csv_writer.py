"""
CSV writer callback for SwanLab.

One row per run — records metadata, config, and the final scalar values for all
metrics, matching the SwanLab table view.
"""

import csv
import os
import time
from pathlib import Path
from typing import TYPE_CHECKING, Dict, List, Optional

from swanlab.proto.swanlab.metric.data.v1.data_pb2 import ScalarRecord
from swanlab.sdk.internal.pkg import safe
from swanlab.sdk.internal.pkg.helper import fmt_run_path
from swanlab.sdk.protocol import Callback

if TYPE_CHECKING:
    from swanlab.sdk.internal.settings import Settings

_METADATA_COLUMNS = [
    "run_id",
    "experiment_name",
    "description",
    "project",
    "workspace",
    "timestamp",
    "logdir",
    "url",
]


def _build_url(settings, path: Optional[str]) -> Optional[str]:
    if path and settings.mode in ("online", "offline"):
        return f"{settings.web_host}{fmt_run_path(path)}"
    return None


def _load_config(run_dir: Path) -> Dict[str, str]:
    config_path = run_dir / "files" / "config.yaml"
    if not config_path.is_file():
        return {}
    import yaml

    with open(config_path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f)
    if not isinstance(data, dict):
        return {}
    result: Dict[str, str] = {}
    for key, val in data.items():
        if isinstance(val, dict) and "value" in val:
            result[str(key)] = str(val["value"])
        else:
            result[str(key)] = str(val)
    return result


def _extract_value(record: ScalarRecord) -> float:
    return record.value.number


class CSVWriter(Callback):
    """CSV writer callback — one row per run with metadata, config, and final scalar values."""

    def __init__(self, dir: str = ".", filename: str = "swanlab_runs.csv"):
        self._save_dir = dir
        self._filename = filename
        self._save_path = os.path.join(dir, filename)

        self._run_dir: Path = Path(".")
        self._metadata: Dict[str, Optional[str]] = {}
        self._config: Dict[str, str] = {}
        self._latest_scalars: Dict[str, float] = {}

    @property
    def name(self) -> str:
        return "CSVWriter"

    def on_run_initialized(self, run_dir: Path, path: str, settings: "Settings", *args, **kwargs) -> None:
        self._run_dir = run_dir
        self._metadata = {}
        self._config = {}
        self._latest_scalars = {}

        self._metadata = {
            "run_id": settings.run.id or "",
            "experiment_name": settings.experiment.name or "",
            "description": settings.experiment.description or "",
            "project": settings.project.name or "",
            "workspace": settings.project.workspace or "",
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
            "logdir": str(run_dir),
            "url": _build_url(settings, path) or "",
        }

    def on_scalar_flush(self, scalar_records: List[ScalarRecord], *args, **kwargs) -> None:
        for record in scalar_records:
            self._latest_scalars[record.key] = _extract_value(record)

    def on_run_finished(self, state: str, error: Optional[str] = None, *args, **kwargs) -> None:
        with safe.block(message="CSVWriter failed to load config"):
            self._config = _load_config(self._run_dir)
        with safe.block(message="CSVWriter failed to write row"):
            self._write_row()

    # ------------------------------------------------------------------
    # CSV writing
    # ------------------------------------------------------------------

    def _write_row(self) -> None:
        os.makedirs(self._save_dir, exist_ok=True)

        config_columns = [f"config/{k}" for k in self._config]
        scalar_columns = list(self._latest_scalars.keys())

        file_exists = os.path.isfile(self._save_path) and os.path.getsize(self._save_path) > 0

        if file_exists:
            existing_headers = self._read_headers()
            new_cols = [c for c in config_columns + scalar_columns if c not in existing_headers]
            if new_cols:
                self._update_header(existing_headers + new_cols)
            all_columns = existing_headers + new_cols
        else:
            all_columns = list(_METADATA_COLUMNS) + config_columns + scalar_columns

        row = self._build_row(all_columns)

        mode = "a" if file_exists else "w"
        with open(self._save_path, mode, newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            if not file_exists:
                writer.writerow(all_columns)
            writer.writerow(row)

    def _build_row(self, columns: List[str]) -> List[str]:
        row: List[str] = []
        for col in columns:
            if col in self._metadata:
                val = self._metadata[col]
                row.append(val if val is not None else "")
            elif col.startswith("config/"):
                row.append(self._config.get(col[7:], ""))
            elif col in self._latest_scalars:
                row.append(str(self._latest_scalars[col]))
            else:
                row.append("")
        return row

    def _read_headers(self) -> List[str]:
        with open(self._save_path, "r", newline="", encoding="utf-8") as f:
            return next(csv.reader(f), [])

    def _update_header(self, new_headers: List[str]) -> None:
        with open(self._save_path, "r", newline="", encoding="utf-8") as f:
            reader = csv.reader(f)
            try:
                next(reader)  # skip old header
            except StopIteration:
                return
            rows = list(reader)
        with open(self._save_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(new_headers)
            writer.writerows(rows)

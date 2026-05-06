from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass(frozen=True)
class RunInfo:
    """Snapshot of experiment metadata passed to Callback hooks."""

    project: str
    workspace: str
    experiment_name: str
    description: Optional[str]
    run_id: str
    run_dir: Path
    path: str
    mode: str
    url: Optional[str]

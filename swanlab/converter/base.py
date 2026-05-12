from abc import ABC, abstractmethod
from typing import List, Optional


class BaseConverter(ABC):
    def __init__(
        self,
        project: Optional[str] = None,
        workspace: Optional[str] = None,
        mode: str = "online",
        log_dir: Optional[str] = None,
        logdir: Optional[str] = None,
        tags: Optional[List[str]] = None,
        resume: bool = False,
    ):
        if logdir is not None:
            import warnings

            warnings.warn(
                "The `logdir` parameter is deprecated, use `log_dir` instead.",
                DeprecationWarning,
                stacklevel=2,
            )
            log_dir = logdir

        self.project = project
        self.workspace = workspace
        self.mode = mode
        self.log_dir = log_dir
        self.tags = tags
        self.resume = resume

    @abstractmethod
    def run(self, **kwargs) -> None: ...

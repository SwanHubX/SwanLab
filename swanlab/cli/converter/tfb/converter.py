from typing import Optional

from swanlab.cli.converter.base import BaseConverter


class TFBConverter(BaseConverter):
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
        self.types = types

    def run(self, convert_dir: str = ".", depth: int = 3, **kwargs) -> None:
        raise NotImplementedError("TensorBoard converter is not yet implemented.")

from swanlab.converter.base import BaseConverter


class MLFlowConverter(BaseConverter):
    def run(self, **kwargs) -> None:
        raise NotImplementedError("MLflow converter is not yet implemented.")

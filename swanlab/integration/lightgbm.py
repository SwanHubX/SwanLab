from __future__ import annotations

from typing import TYPE_CHECKING

import swanlab
from swanlab import Callback

if TYPE_CHECKING:
    from typing import Any, Dict, List, NamedTuple, Tuple, Union

    _EvalResultTuple = Union[Tuple[str, str, float, bool], Tuple[str, str, float, bool, float]]

    class CallbackEnv(NamedTuple):
        model: Any
        params: Dict
        iteration: int
        begin_iteration: int
        end_iteration: int
        evaluation_result_list: List[_EvalResultTuple]


class SwanLabCallback(Callback):
    def __init__(self, log_params: bool = True) -> None:
        self._log_params = log_params
        self._initialized = False

    @property
    def name(self) -> str:
        return "swanlab-integration-lightgbm"

    def _init(self, env: CallbackEnv) -> None:
        if self._initialized:
            return
        run = swanlab.get_run()
        if run is None:
            return
        self._initialized = True
        run.config["FRAMEWORK"] = "lightgbm"
        if self._log_params:
            run.config.update(env.params)

    def __call__(self, env: CallbackEnv) -> None:
        self._init(env)

        for item in env.evaluation_result_list:
            if len(item) == 4:
                if isinstance(item[1], str):
                    data_name, eval_name, result = item[:3]
                    swanlab.log({f"{data_name}_{eval_name}": result}, step=env.iteration)
                else:
                    eval_name, result, _, stdv = item
                    swanlab.log(
                        {f"{eval_name}-mean": result, f"{eval_name}-stdv": stdv},
                        step=env.iteration,
                    )
            elif len(item) == 5:
                data_name, eval_name, result, _, stdv = item
                swanlab.log(
                    {
                        f"{data_name}_{eval_name}-mean": result,
                        f"{data_name}_{eval_name}-stdv": stdv,
                    },
                    step=env.iteration,
                )

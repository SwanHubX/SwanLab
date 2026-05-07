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
        self._initialized = True
        run = swanlab.get_run()
        if run is None:
            return
        run.config["FRAMEWORK"] = "lightgbm"
        if self._log_params:
            run.config.update(env.params)

    def __call__(self, env: CallbackEnv) -> None:
        self._init(env)

        for item in env.evaluation_result_list:
            if len(item) == 4:
                data_name, eval_name, result = item[:3]
                swanlab.log({data_name + "_" + eval_name: result}, step=env.iteration)
            else:
                data_name, eval_name = item[1].split()
                res_mean = item[2]
                res_stdv = item[4]
                swanlab.log(
                    {
                        data_name + "_" + eval_name + "-mean": res_mean,
                        data_name + "_" + eval_name + "-stdv": res_stdv,
                    },
                    step=env.iteration,
                )

from typing import TYPE_CHECKING
import lightgbm  # type: ignore
from lightgbm import Booster
import swanlab


if TYPE_CHECKING:
    from typing import Any, Dict, List, NamedTuple, Tuple, Union

    # Note: upstream lightgbm has this defined incorrectly
    _EvalResultTuple = Union[
        Tuple[str, str, float, bool], Tuple[str, str, float, bool, float]
    ]

    class CallbackEnv(NamedTuple):
        model: Any
        params: Dict
        iteration: int
        begin_interation: int
        end_iteration: int
        evaluation_result_list: List[_EvalResultTuple]



class SwanLabCallback:
    def __init__(self, log_params: bool = True) -> None:
        self.order = 20
        self.before_iteration = False
        self.log_params = log_params

    def _init(self, env: "CallbackEnv") -> None:
        swanlab.config["FRAMEWORK"] = "lightgbm"
        if self.log_params:
            swanlab.config.update(env.params)

    def __call__(self, env: "CallbackEnv") -> None:
        if env.iteration == env.begin_iteration:  # type: ignore
            self._init(env)

        for item in env.evaluation_result_list:
            if len(item) == 4:
                data_name, eval_name, result = item[:3]
                swanlab.log(
                    {data_name + "_" + eval_name: result},
                )
            else:
                data_name, eval_name = item[1].split()
                res_mean = item[2]
                res_stdv = item[4]
                swanlab.log(
                    {
                        data_name + "_" + eval_name + "-mean": res_mean,
                        data_name + "_" + eval_name + "-stdv": res_stdv,
                    },
                )

        swanlab.log({"iteration": env.iteration})
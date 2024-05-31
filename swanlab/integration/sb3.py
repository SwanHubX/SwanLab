"""
For adaptation to the stable_baseline3 framework. Detailed usage are as follows:
------trian.py in stable_baseline3------
import swanlab
import gymnasium as gym
from stable_baselines3 import PPO
from stable_baselines3.common.monitor import Monitor
from stable_baselines3.common.vec_env import DummyVecEnv
from swanlab.integration.sb3 import SwanLabCallback

config = {
    "policy_type": "MlpPolicy",
    "total_timesteps": 25000,
    "env_name": "CartPole-v1",
}

def make_env():
    env = gym.make(config["env_name"], render_mode="rgb_array")
    env = Monitor(env)
    return env

env = DummyVecEnv([make_env])
model = PPO(
    config["policy_type"],
    env,
    verbose=1,
)

model.learn(
    total_timesteps=config["total_timesteps"],
    callback=SwanLabCallback(
        project="PPO",
        experiment_name="MlpPolicy",
        verbose=2,
    ),
)

swanlab.finish()
---------------------------------
"""

import swanlab
from typing import Optional, Dict, Any, Union, Tuple
from stable_baselines3.common.callbacks import BaseCallback
from stable_baselines3.common.logger import KVWriter, Logger


class SwanLabOutputFormat(KVWriter):
    def __init__(self, swanlab_callback):
        self.swanlab_callback = swanlab_callback

    def write(
        self,
        key_values: Dict[str, Any],
        key_excluded: Dict[str, Union[str, Tuple[str, ...]]],
        step: int = 0,
    ) -> None:
        for (key, value), (_, excluded) in zip(sorted(key_values.items()), sorted(key_excluded.items())):
            # 如果是标量指标
            if isinstance(value, (int, float)):
                # 记录指标
                self.swanlab_callback.experiment.log({key: value}, step=step)


class SwanLabCallback(BaseCallback):
    def __init__(
        self,
        project: Optional[str] = None,
        workspace: Optional[str] = None,
        experiment_name: Optional[str] = None,
        description: Optional[str] = None,
        logdir: Optional[str] = None,
        mode: Optional[bool] = None,
        verbose: int = 0,
        **kwargs: Any,
    ):
        super().__init__(verbose)
        self._run = None

        self._swanlab_init: Dict[str, Any] = {
            "project": project,
            "workspace": workspace,
            "experiment_name": experiment_name,
            "description": description,
            "logdir": logdir,
            "mode": mode,
        }

        self._swanlab_init.update(**kwargs)

        self._project = self._swanlab_init.get("project")
        self._workspace = self._swanlab_init.get("workspace")
        self._experiment_name = self._swanlab_init.get("experiment_name")
        self._description = self._swanlab_init.get("decsription")
        self._logdir = self._swanlab_init.get("logdir")
        self._mode = self._swanlab_init.get("mode")

    def _init_callback(self) -> None:
        args = {"algo": type(self.model).__name__}
        for key in self.model.__dict__:
            if type(self.model.__dict__[key]) in [float, int, str]:
                args[key] = self.model.__dict__[key]
            else:
                args[key] = str(self.model.__dict__[key])

        self.setup(config=args)

        loggers = Logger(
            folder=None,
            output_formats=[SwanLabOutputFormat(self)],
        )

        self.model.set_logger(loggers)

    @property
    def experiment(self):
        if swanlab.get_run() is None:
            self.setup()
        return self._run

    def setup(self, config=None):
        if swanlab.get_run() is None:
            self._run = swanlab.init(**self._swanlab_init)
        else:
            self._run = swanlab.get_run()

        if config:
            self._run.config.update(config)

    def _on_step(self) -> bool:
        return True

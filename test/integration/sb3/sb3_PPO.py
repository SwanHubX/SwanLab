from tutils import open_dev_mode
import swanlab

swanlab.login(open_dev_mode())

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

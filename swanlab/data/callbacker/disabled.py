"""
@author: cunyue
@file: disabled.py
@time: 2025/6/24 14:21
@description: disabled 模式回调
"""

from swanlab.data.callbacker.callback import SwanLabRunCallback

from .. import namer as N


class DisabledCallback(SwanLabRunCallback):
    def on_init(self, proj_name: str, workspace: str, public: bool = None, logdir: str = None, *args, **kwargs):
        self.run_store.run_name = "run-disabled"
        self.run_store.run_colors = N.generate_colors(0)
        self.run_store.run_id = N.generate_run_id()
        self.run_store.new = True

    def __str__(self):
        return "SwanLabDisabledCallback"

    def on_stop(self, error: str = None, *args, **kwargs):
        self.porter.close()

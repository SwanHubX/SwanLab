"""
Writer plugin for SwanLab.
Used for writing experiment metadata to CSV and online notes.
"""
from swankit.callback import SwanKitCallback
import swanlab
from typing import Optional, Tuple
import csv
import os
import time

class CSVWriter(SwanKitCallback):
    """CSV writer for SwanLab."""

    def __init__(self, dir: str = None, filename: str = "swanlab_metadata.csv"):
        if dir is None:
            dir = os.getcwd()
        self.save_path = os.path.join(dir, filename)
        self.csv_file = open(self.save_path, "w")
        self.csv_writer = csv.writer(self.csv_file)

    def on_init(self, proj_name: str, workspace: str, logdir: str, **kwargs):
        self.project = proj_name
        self.workspace = workspace
    
    def before_init_experiment(self, run_id: str, exp_name: str, description: str, num: int, colors: Tuple[str, str]):
        self.run_id = run_id
        self.exp_name = exp_name
        self.description = description
    
    def on_run(self):
        self.csv_writer.writerow(["datetime", "run_id", "project", "workspace", "exp_name", "description","logdir", "url"])
        timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        run = swanlab.get_run()
        self.logdir = run.public.swanlog_dir
        if run.public.cloud.project_url is not None:
            # self.project_url = run.public.cloud.project_url
            self.experiment_url = run.public.cloud.experiment_url
        
        self.csv_writer.writerow([timestamp, self.run_id, self.project, self.workspace, self.exp_name, self.description, self.logdir, self.experiment_url])

    def on_stop(self, error: Optional[str] = None):
        """关闭CSV文件"""
        self.csv_file.close()
    
    def __str__(self):
        return "CSVWriter"

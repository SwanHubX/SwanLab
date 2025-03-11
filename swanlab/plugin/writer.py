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

    def __init__(self, dir: str = None, filename: str = "swanlab_run_metadata.csv"):
        if dir is None:
            dir = os.getcwd()
        self.save_path = os.path.join(dir, filename)
        
        # Check if file exists to determine mode
        file_exists = os.path.isfile(self.save_path)
        self.csv_file = open(self.save_path, "a" if file_exists else "w")
        self.csv_writer = csv.writer(self.csv_file)
        self.file_exists = file_exists

    def on_init(self, proj_name: str, workspace: str, logdir: str, **kwargs):
        self.project = proj_name
        self.workspace = workspace
    
    def before_init_experiment(self, run_id: str, exp_name: str, description: str, num: int, colors: Tuple[str, str]):
        self.run_id = run_id
        self.exp_name = exp_name
        self.description = description
    
    def on_run(self):
        # Get run information
        run = swanlab.get_run()
        config = run.config
        self.logdir = run.public.swanlog_dir
        timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        
        # Set experiment URL if available
        self.experiment_url = None
        if run.public.cloud.project_url is not None:
            self.experiment_url = run.public.cloud.experiment_url
        
        # Prepare headers and base row data
        headers = ["project",  "exp_name", "description", "datetime", "run_id",  "workspace", "logdir", "url"]
        row_data = [self.project, self.exp_name, self.description, timestamp, self.run_id, self.workspace, self.logdir, self.experiment_url]
        
        # Add config keys as columns
        # TODO: Consider the case where the config contains commas.
        if config:
            for key in config:
                headers.append(key)
                row_data.append(config[key])
                
        # Write to CSV - only write headers if file is new
        if not self.file_exists:
            self.csv_writer.writerow(headers)
        self.csv_writer.writerow(row_data)


    def on_stop(self, error: Optional[str] = None):
        """关闭CSV文件"""
        self.csv_file.close()
    
    def __str__(self):
        return "CSVWriter"

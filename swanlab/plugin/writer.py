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
try:
    import pandas as pd
except ImportError:
    raise ImportError("pandas is required for CSVWriter. Please install it using 'pip install pandas'.")

class CSVWriter(SwanKitCallback):
    """CSV writer for SwanLab."""
    def __init__(self, dir: str = None, filename: str = "swanlab_run.csv"):
        
        if dir is None:
            dir = os.getcwd()
        
        # Create directory if it doesn't exist
        os.makedirs(dir, exist_ok=True)
        
        self.save_path = os.path.join(dir, filename)
        self.file_exists = os.path.isfile(self.save_path)
        self.df = self._initialize_dataframe()
        self.original_headers = list(self.df.columns) if self.file_exists else []

    def _initialize_dataframe(self):
        """Initialize the DataFrame based on file existence."""
        if self.file_exists:
            try:
                df = pd.read_csv(self.save_path)
                if df.empty:
                    return self._create_empty_dataframe()
                return df
            except pd.errors.EmptyDataError:
                return self._create_empty_dataframe()
        else:
            return self._create_empty_dataframe()

    def _create_empty_dataframe(self):
        """Create an empty DataFrame with default columns."""
        self.file_exists = False  # Treat as a new file
        return pd.DataFrame(columns=["project", "exp_name", "description", "datetime", "run_id", "workspace", "logdir", "url"])

    def on_init(self, proj_name: str, workspace: str, logdir: str = None, *args, **kwargs):
        """Initialize project and workspace information."""
        self.project = proj_name
        self.workspace = workspace
    
    def before_init_experiment(
        self,
        run_id: str,
        exp_name: str,
        description: str,
        num: int,
        colors: Tuple[str, str],
        *args,
        **kwargs,
    ):
        """Set experiment details before initialization."""
        self.run_id = run_id
        self.exp_name = exp_name
        self.description = description
    
    def on_run(self, *args, **kwargs):
        """Handle actions to perform on run."""
        run = swanlab.get_run()
        config = run.config
        self.logdir = run.public.swanlog_dir
        timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        
        # Set experiment URL if available
        self.experiment_url = run.public.cloud.experiment_url if run.public.cloud.project_url else None
        
        headers = ["project", "exp_name", "description", "datetime", "run_id", "workspace", "logdir", "url"]
        row_data = [self.project, self.exp_name, self.description, timestamp, self.run_id, self.workspace, self.logdir, self.experiment_url]
        
        if not self.file_exists:
            self._handle_new_file(config, headers, row_data)
        else:
            self._handle_existing_file(config, headers, row_data)

    def _handle_new_file(self, config, headers, row_data):
        """Handle writing to a new file."""
        if config:
            for key in config:
                headers.append("config/" + key)
                row_data.append(config[key])
        
        self.df = pd.DataFrame([], columns=headers)
        updated_df = pd.DataFrame([row_data], columns=headers)
        updated_df.to_csv(self.save_path, index=False)
        
        self.headers = headers
        self.last_row_data = row_data

    def _handle_existing_file(self, config, headers, row_data):
        """Handle writing to an existing file."""
        headers = self.original_headers.copy()
        headers_metadata = headers[:8]
        headers_config = [header for header in headers[8:] if header.startswith("config/")]
        headers_config_dict = {key: {"value": " ", "index": headers_config.index(key)} for key in headers_config}
        
        for key in config:
            if "config/" + key not in headers_config_dict:
                headers_config_dict["config/" + key] = {"value": config[key], "index": len(headers_config)}
                headers_config.append("config/" + key)
            else:
                headers_config_dict["config/" + key]["value"] = config[key]
                
        headers_config_list = [headers_config_dict[header]["value"] for header in headers_config]
        headers = headers_metadata + headers_config
        row_data = row_data + headers_config_list
        
        new_row_df = pd.DataFrame([row_data], columns=headers)
        updated_df = pd.concat([self.df, new_row_df], ignore_index=True)
        updated_df.to_csv(self.save_path, index=False)
        
        self.headers = headers
        self.last_row_data = row_data

    def on_log(self, data: dict, step: int | None, *args, **kwargs):
        """Log data during the run."""
        if not data:
            return
        
        new_headers = self.headers.copy()
        for key in data:
            if key not in new_headers:
                new_headers.append(key)
        self.headers = new_headers
        
        new_row = self.last_row_data.copy()
        
        for key in data:
            if key in self.headers:
                idx = self.headers.index(key)
                if idx < len(new_row):
                    new_row[idx] = data[key]
                else:
                    while len(new_row) < idx:
                        new_row.append(None)
                    new_row.append(data[key])
            else:
                self.headers.append(key)
                new_row.append(data[key])
        
        new_row_df = pd.DataFrame([new_row], columns=self.headers)
        updated_df = pd.concat([self.df, new_row_df], ignore_index=True)
        updated_df.to_csv(self.save_path, index=False)
    
    def __str__(self):
        return "CSVWriter"
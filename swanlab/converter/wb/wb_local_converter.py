# -*- coding: utf-8 -*-
"""
This code references: https://github.com/yuz1wan/localWandb2cloudSwanlab, thanks to @yuz1wan.

------example.py------
from swanlab.converter import WandbLocalConverter

wb_local_converter = WandbLocalConverter()
wb_local_converter.run(root_wandb_dir="WANDB_LOCAL_DIR", wandb_run_dir="WANDB_LOCAL_RUN_DIR")

------command------
swanlab convert -t wandb-local --wb-dir ./wandb --wb-run-dir run-1234567890
"""
import argparse
import glob
import json
import os
import re
import traceback
from typing import Optional, List

import yaml
import swanlab
from swanlab.log import swanlog as swl

# Dependency for converting Protobuf objects to dictionaries
import google.protobuf.json_format as protobuf_json

try:
    from wandb.sdk.internal.datastore import DataStore
    from wandb.proto import wandb_internal_pb2
except ImportError:
    # We add a print here because if this fails, swanlog might not be available.
    print("Error: Wandb is required to parse local log files. Please install it with 'pip install wandb'.")
    exit(1)


class WandbLocalConverter:
    """
    A robust converter to migrate local wandb run data to SwanLab.

    This class reliably parses all data from a .wandb file by using the
    `ds.scan_data()` method and converting Protobuf messages to dictionaries
    with `MessageToDict`.
    """
    def __init__(
        self,
        project: Optional[str] = None,
        workspace: Optional[str] = None,
        mode: Optional[str] = "cloud",
        tags: Optional[List[str]] = None,
        logdir: Optional[str] = None,
        **kwargs,
    ):
        """
        Initializes the converter with SwanLab project details.

        Args:
            project (Optional[str]): The SwanLab project name. If None, it will be inferred from the wandb run.
            workspace (Optional[str]): The SwanLab workspace.
            logdir (Optional[str]): The local directory to save SwanLab logs.
        """
        self.project = project
        self.workspace = workspace
        self.mode = mode
        self.tags = tags
        self.logdir = logdir

    def _find_run_dirs(self, root_wandb_dir: str, wandb_run_dir: Optional[str] = None) -> list[str]:
        """Finds all wandb run directories within a given root directory."""
        if wandb_run_dir:
            patterns = [
                os.path.join(root_wandb_dir, wandb_run_dir)
            ]
            found_path = os.path.join(root_wandb_dir, wandb_run_dir)
        else:
            patterns = [
                os.path.join(root_wandb_dir, "run-*/"),
                os.path.join(root_wandb_dir, "offline-run-*/")
            ]
            found_path = root_wandb_dir

        run_dirs = []
        for pattern in patterns:
            run_dirs.extend(glob.glob(pattern))
        if not run_dirs:
            swl.warning(f"No wandb run directories found in '{found_path}'.")
        return run_dirs

    def _unpack_key_value_json_list(self, items: list) -> dict:
        """A helper function to unpack wandb's specific key-value_json list format."""
        if not (isinstance(items, list) and items and 'value_json' in items[0]):
            return items

        mapping = {}
        for item in items:
            key = item.get('key') or '/'.join(item.get('nested_key', []))
            if not key:
                continue
            try:
                value = json.loads(item['value_json'])
                mapping[key] = value
            except json.JSONDecodeError:
                swl.warning(f"Could not decode json for key '{key}': {item['value_json']}")
        return mapping

    def _parse_run(self, run_dir: str):
        """Parses a single wandb run directory and converts it to a SwanLab run."""
        wandb_files = glob.glob(os.path.join(run_dir, "*.wandb"))
        if not wandb_files:
            swl.error(f"No .wandb file found in {run_dir}. Skipping.")
            return

        wandb_file_path = wandb_files[0]
        files_root_dir = os.path.join(run_dir, "files")

        ds = DataStore()
        ds.open_for_scan(wandb_file_path)

        swanlab_run = None
        run_metadata = {'id': None, 'name': None, 'notes': None, 'project': None}
        run_config = {}

        def initialize_swanlab_run_if_needed():
            """Lazy initializer for the SwanLab run."""
            nonlocal swanlab_run
            if swanlab_run is not None:
                return

            if run_metadata['id'] is None:
                match = re.search(r"-([a-zA-Z0-9]+)$", os.path.basename(run_dir.rstrip('/\\')))
                if match:
                    run_metadata['id'] = match.group(1)
                else:
                    raise RuntimeError(f"Fatal: Cannot determine run_id for {run_dir}.")
            if run_metadata['name'] is None:
                run_metadata['name'] = run_metadata['id']

            swl.info(f"Converting run: {run_metadata['name']} (ID: {run_metadata['id']})")
            swanlab_run = swanlab.init(
                project=self.project or run_metadata['project'] or "wandb-imports",
                workspace=self.workspace,
                experiment_name=run_metadata['name'],
                description=run_metadata['notes'],
                logdir=self.logdir,
                mode=self.mode,
                tags=self.tags,
            )

            try:
                config_path = os.path.join(files_root_dir, "config.yaml")
                with open(config_path, 'r', encoding='utf-8') as f:
                    wandb_config_yaml = yaml.safe_load(f)
                    cleaned_config = {k: v['value'] for k, v in wandb_config_yaml.items() if isinstance(v, dict) and 'value' in v}
                    run_config.update(cleaned_config)
            except FileNotFoundError:
                pass
            
            final_config = {"wandb_run_id": run_metadata['id'], "wandb_local_path": run_dir}
            final_config.update(run_config)
            swanlab_run.config.update(final_config)

        # Core Logic: Scan the .wandb file record by record
        while True:
            record_bin = ds.scan_data()
            if record_bin is None:
                break

            record_pb = wandb_internal_pb2.Record()
            record_pb.ParseFromString(record_bin)

            record_dict = protobuf_json.MessageToDict(record_pb, preserving_proto_field_name=True)
            record_type = record_pb.WhichOneof("record_type")
            data_dict = record_dict.get(record_type, {})

            if not data_dict:
                continue

            # Process the converted dictionary data based on record type
            if record_type == "run":
                run_metadata.update({
                    'id': data_dict.get('run_id'),
                    'name': data_dict.get('display_name'),
                    'notes': data_dict.get('notes'),
                    'project': data_dict.get('project')
                })
                config_items = data_dict.get('config', {}).get('update', [])
                run_config.update(self._unpack_key_value_json_list(config_items))
            elif record_type == "config":
                config_items = data_dict.get('update', [])
                run_config.update(self._unpack_key_value_json_list(config_items))
            elif record_type == "history":
                initialize_swanlab_run_if_needed()
                log_dict = {}
                history_data = self._unpack_key_value_json_list(data_dict.get('item', []))
                
                for key, value in history_data.items():
                    # 把指标中的runtime、_step、_timestamp都放到_wandb分组里
                    if key == "_runtime" or key == "_step" or key == "_timestamp":
                        key = "_wandb/" + key
                    if isinstance(value, (int, float)):
                        log_dict[key] = value
                    elif isinstance(value, dict) and "_type" in value:
                        media_type = value["_type"]
                        path = os.path.join(files_root_dir, value.get("path", ""))
                        if not os.path.exists(path):
                            continue
                        
                        if media_type == "image-file":
                            log_dict[key] = swanlab.Image(path)
                
                step = int(data_dict.get('_step', 0))
                if log_dict:
                    swanlab_run.log(log_dict, step=step)
            elif record_type == "summary":
                initialize_swanlab_run_if_needed()
                log_dict = {}
                summary_data = self._unpack_key_value_json_list(data_dict.get('update', []))
                
                for key, value in summary_data.items():
                    if key == "_runtime" or key == "_step" or key == "_timestamp":
                        key = "_wandb/" + key
                    if isinstance(value, (int, float)):
                        log_dict[key] = value
                if log_dict:
                    swanlab_run.log(log_dict)

        if swanlab_run:
            swl.info(f"Finished converting run: {run_metadata['name']}")
            swanlab_run.finish()
        else:
            try:
                initialize_swanlab_run_if_needed()
                if swanlab_run:
                    swl.warning(f"Run in {run_dir} has no metrics, but its config was saved.")
                    swanlab_run.finish()
            except RuntimeError as e:
                swl.error(str(e))

    def run(self, root_wandb_dir: str, wandb_run_dir: Optional[str] = None):
        """The main entry point to start the conversion process."""
        swl.info(f"Starting import from wandb directory: {os.path.abspath(root_wandb_dir)}")
        run_dirs = self._find_run_dirs(root_wandb_dir, wandb_run_dir)
        
        if not run_dirs:
            swl.info("No runs found to import. Exiting.")
            return

        swl.info(f"Found {len(run_dirs)} runs to import.")
        
        for i, run_dir in enumerate(run_dirs):
            swl.info(f"--- Processing run {i+1}/{len(run_dirs)} ---")
            try:
                self._parse_run(run_dir)
            except Exception as e:
                swl.error(f"FATAL: Failed to convert run in {run_dir}: {e}")
                traceback.print_exc()
        
        swl.info("="*20)
        swl.info("All runs processed.")
        swl.info("="*20)


def main():
    """Main function to parse arguments and run the converter."""
    parser = argparse.ArgumentParser(
        description="Import local Wandb runs into SwanLab.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument(
        "--wandb_dir",
        "-d",
        type=str,
        default="./wandb",
        required=True,
        help="Path to the local 'wandb' directory containing run folders (e.g., './wandb')."
    )
    
    parser.add_argument(
        "--wandb_run_dir",
        "-r",
        type=str,
        default=None,
        help="Path to the local 'wandb' directory containing run folders (e.g., './wandb')."
    )
    
    parser.add_argument(
        "--project",
        "-p",
        type=str,
        default=None,
        help="Specify the SwanLab project name. If not set, it's inferred from the wandb run."
    )
    
    parser.add_argument(
        "--workspace",
        "-w",
        type=str,
        default=None,
        help="Specify the SwanLab workspace."
    )
    
    parser.add_argument(
        "--logdir",
        "-l",
        type=str,
        default=None,
        help="Specify the root directory for SwanLab logging. Defaults to './swanlog'."
    )
    
    args = parser.parse_args()
    
    swl.info("SwanLab Importer for Wandb")
    swl.info("---------------------------")
    swl.info(f"Source Wandb Directory: {os.path.abspath(args.wandb_dir)}")
    swl.info(f"Source Wandb Run Directory: {args.wandb_run_dir or '(all runs)'}")
    swl.info(f"Target SwanLab Project: {args.project or '(inferred from run)'}")
    swl.info(f"Target SwanLab Workspace: {args.workspace or '(default)'}")
    swl.info(f"SwanLab Log Directory: {args.logdir or './swanlog'}")
    swl.info("---------------------------")

    converter = WandbLocalConverter(
        project=args.project,
        workspace=args.workspace,
        logdir=args.logdir
    )
    converter.run(root_wandb_dir=args.wandb_dir, wandb_run_dir=args.wandb_run_dir)


if __name__ == '__main__':
    main()
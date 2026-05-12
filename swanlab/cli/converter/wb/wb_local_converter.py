import gc
import glob
import os
import re
from typing import List, Optional

import click

from swanlab import vendor
from swanlab.cli.converter.base import BaseConverter
from swanlab.cli.converter.utils import json_loads, proto_items_to_dict, validate_path
from swanlab.sdk.internal.pkg import safe

GC_INTERVAL: int = 10_000


class WandbLocalConverter(BaseConverter):
    """Convert local W&B .wandb files to SwanLab via public SDK API."""

    def run(  # type: ignore[override]
        self, root_wandb_dir: str = "./wandb", wandb_run_dir: Optional[str] = None, wb_run_id: Optional[str] = None
    ) -> None:
        click.echo(f"Starting import from wandb directory: {os.path.abspath(root_wandb_dir)}")
        run_dirs = self._find_run_dirs(root_wandb_dir, wandb_run_dir)

        if not run_dirs:
            click.echo("No runs found to import.")
            return

        click.echo(f"Found {len(run_dirs)} run(s) to import.")

        for i, rd in enumerate(run_dirs):
            click.echo(f"\n--- Processing run {i + 1}/{len(run_dirs)} ---")
            with safe.block(message=f"Failed to convert run in {rd}"):
                self._parse_run(rd, wb_run_id)

        click.echo("\nAll runs processed.")

    def _find_run_dirs(self, root_wandb_dir: str, wandb_run_dir: str | None = None) -> List[str]:
        if wandb_run_dir:
            patterns = [os.path.join(root_wandb_dir, wandb_run_dir)]
            found_path = os.path.join(root_wandb_dir, wandb_run_dir)
        else:
            patterns = [
                os.path.join(root_wandb_dir, "run-*/"),
                os.path.join(root_wandb_dir, "offline-run-*/"),
            ]
            found_path = root_wandb_dir

        run_dirs: list[str] = []
        for pattern in patterns:
            run_dirs.extend(glob.glob(pattern))
        if not run_dirs:
            click.echo(f"No wandb run directories found in '{found_path}'.")
        return run_dirs

    def _filter_text_columns(self, columns, data):
        non_text_indices = {
            i for row in data for i, cell in enumerate(row) if isinstance(cell, dict) and "_type" in cell
        }
        text_indices = [i for i in range(len(columns)) if i not in non_text_indices]
        return (
            [columns[i] for i in text_indices],
            [[row[i] for i in text_indices if i < len(row)] for row in data],
        )

    def _make_table_echarts(self, columns, data):
        try:
            from pyecharts.components import Table as PyEchartsTable
        except ImportError:
            click.echo("Warning: pyecharts is required for table conversion. Skipping table.")
            return None

        filtered_cols, filtered_data = self._filter_text_columns(columns, data)
        if not filtered_cols:
            return None

        table = PyEchartsTable()
        table.add(filtered_cols, filtered_data)
        import swanlab

        return swanlab.ECharts(table)

    def _parse_run(self, run_dir: str, wb_run_id: Optional[str] = None) -> None:
        import swanlab

        DataStore = vendor.wandb.sdk.internal.datastore.DataStore  # type: ignore
        wandb_internal_pb2 = vendor.wandb.proto.wandb_internal_pb2  # type: ignore

        wandb_files = glob.glob(os.path.join(run_dir, "*.wandb"))
        if not wandb_files:
            click.echo(f"Error: No .wandb file found in {run_dir}. Skipping.")
            return

        wandb_file_path = wandb_files[0]
        files_root_dir = os.path.join(run_dir, "files")

        ds = DataStore()
        ds.open_for_scan(wandb_file_path)

        swanlab_run = None
        run_metadata: dict = {"id": None, "name": None, "notes": None, "project": None}
        run_config: dict = {}

        def initialize_swanlab_run():
            nonlocal swanlab_run
            if swanlab_run is not None:
                return

            if run_metadata["id"] is None:
                match = re.search(r"-([a-zA-Z0-9]+)$", os.path.basename(run_dir.rstrip("/\\")))
                if match:
                    run_metadata["id"] = match.group(1)
                else:
                    raise RuntimeError(f"Cannot determine run_id for {run_dir}.")
            if run_metadata["name"] is None:
                run_metadata["name"] = run_metadata["id"]

            click.echo(f"Converting run: {run_metadata['name']} (ID: {run_metadata['id']})")
            swanlab_run = swanlab.init(
                project=self.project or run_metadata["project"] or "wandb-imports",
                workspace=self.workspace,
                name=run_metadata["name"],
                description=run_metadata["notes"],
                log_dir=self.log_dir,
                mode=self.mode,  # type: ignore[arg-type]
                tags=self.tags,
                reinit=True,
                **({"resume": True, "id": wb_run_id} if (self.resume and wb_run_id) else {}),
            )

            # Load config.yaml from wandb run files (non-critical)
            with safe.block(message="Failed to load wandb config.yaml", write_to_tty=False):
                import yaml

                config_path = os.path.join(files_root_dir, "config.yaml")
                with open(config_path, "r", encoding="utf-8") as f:
                    wandb_config_yaml = yaml.safe_load(f)
                    cleaned = {
                        k: v["value"] for k, v in wandb_config_yaml.items() if isinstance(v, dict) and "value" in v
                    }
                    run_config.update(cleaned)

            final_config = {"wandb_run_id": run_metadata["id"], "wandb_local_path": run_dir}
            final_config.update(run_config)
            swanlab_run.config.update(final_config)

        record_count = 0

        while True:
            record_bin = ds.scan_data()
            if record_bin is None:
                break

            record_pb = wandb_internal_pb2.Record()
            record_pb.ParseFromString(record_bin)

            record_type = record_pb.WhichOneof("record_type")

            if record_type == "run":
                run = record_pb.run
                run_metadata.update(
                    {
                        "id": run.run_id or None,
                        "name": run.display_name or None,
                        "notes": run.notes or None,
                        "project": run.project or None,
                    }
                )
                run_config.update(proto_items_to_dict(run.config.update))

            elif record_type == "config":
                run_config.update(proto_items_to_dict(record_pb.config.update))

            elif record_type == "history":
                initialize_swanlab_run()
                scalar_dict: dict = {}
                media_dict: dict = {}
                grouped_items: dict = {}
                step = 0

                for item in record_pb.history.item:
                    key = item.key or "/".join(item.nested_key)
                    if not key:
                        continue
                    value_json = item.value_json

                    if key == "_step":
                        with safe.block(message=None, write_to_tty=False):
                            step = int(float(value_json))
                        continue
                    if key.startswith("_"):
                        continue

                    # Fast path: scalar
                    try:
                        scalar_dict[key] = float(value_json)
                        continue
                    except (ValueError, TypeError):
                        pass

                    # Slow path: media objects
                    with safe.block(message=None, write_to_tty=False):
                        value = json_loads(value_json)
                        if isinstance(value, int):
                            scalar_dict[key] = float(value)
                            continue
                        if isinstance(value, dict) and "path" in value:
                            validated_path = validate_path(files_root_dir, value["path"])
                            if validated_path and os.path.exists(validated_path):
                                if value.get("_type") == "image-file":
                                    media_dict[key] = swanlab.Image(validated_path)
                                    continue
                                elif value.get("_type") == "audio-file":
                                    media_dict[key] = swanlab.Audio(validated_path)
                                    continue

                    # Grouped path: tables etc.
                    if "/" in key:
                        parts = key.split("/", 1)
                        base_key = parts[0]
                        if base_key not in grouped_items:
                            grouped_items[base_key] = {}
                        with safe.block(message=None, write_to_tty=False):
                            grouped_items[base_key][parts[1]] = json_loads(value_json)

                # Process grouped items (tables, grouped media)
                for base_key, props in grouped_items.items():
                    if props.get("_type") == "table-file" and "path" in props:
                        validated_path = validate_path(files_root_dir, props["path"])
                        if validated_path and os.path.exists(validated_path):
                            with safe.block(message=f"Failed to parse table from {validated_path}"):
                                with open(validated_path, "rb") as f:
                                    table_data = json_loads(f.read())
                                columns = table_data.get("columns", [])
                                data = table_data.get("data", [])
                                echarts_obj = self._make_table_echarts(columns, data)
                                if echarts_obj is not None:
                                    media_dict[base_key] = echarts_obj
                    elif props.get("_type") == "image-file" and "path" in props:
                        validated_path = validate_path(files_root_dir, props["path"])
                        if validated_path and os.path.exists(validated_path):
                            media_dict[base_key] = swanlab.Image(validated_path)
                    elif props.get("_type") == "audio-file" and "path" in props:
                        validated_path = validate_path(files_root_dir, props["path"])
                        if validated_path and os.path.exists(validated_path):
                            media_dict[base_key] = swanlab.Audio(validated_path)

                if scalar_dict or media_dict:
                    if scalar_dict and swanlab_run:
                        swanlab_run.log(scalar_dict, step=step)
                    if media_dict and swanlab_run:
                        swanlab_run.log(media_dict, step=step)
                del scalar_dict, media_dict

            record_count += 1
            if record_count % GC_INTERVAL == 0:
                gc.collect()

        if swanlab_run:
            click.echo(f"Finished converting run: {run_metadata['name']}")
            swanlab_run.finish()
        else:
            try:
                initialize_swanlab_run()
                if swanlab_run:
                    click.echo(f"Warning: Run in {run_dir} has no metrics, config was saved.")
                    swanlab_run.finish()
            except RuntimeError as e:
                click.echo(f"Error: {e}")

        del ds
        gc.collect()

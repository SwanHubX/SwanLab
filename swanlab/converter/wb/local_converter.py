import gc
import glob
import os
import re
from typing import List, Optional

from swanlab import vendor
from swanlab.converter.base import BaseConverter
from swanlab.converter.wb.utils import json_loads, proto_items_to_dict, validate_path
from swanlab.sdk.internal.pkg import safe

GC_INTERVAL: int = 10_000


class WandbLocalConverter(BaseConverter):
    """Convert local W&B .wandb files to SwanLab via public SDK API.

    Example::

        wb = WandbLocalConverter(project="my-project")
        wb.run(root_wandb_dir="./wandb", wandb_run_dir="run-abc123")
    """

    def __init__(
        self,
        project: Optional[str] = None,
        workspace: Optional[str] = None,
        mode: str = "online",
        log_dir: Optional[str] = None,
        logdir: Optional[str] = None,
        tags: Optional[List[str]] = None,
        resume: bool = False,
        root_wandb_dir: str = "./wandb",
        wandb_run_dir: Optional[str] = None,
    ):
        super().__init__(
            project=project,
            workspace=workspace,
            mode=mode,
            log_dir=log_dir,
            logdir=logdir,
            tags=tags,
            resume=resume,
        )
        # Local W&B defaults — can be overridden at run() time
        self._root_wandb_dir = root_wandb_dir
        self._wandb_run_dir = wandb_run_dir

    def run(  # type: ignore[override]
        self,
        root_wandb_dir: Optional[str] = None,
        wandb_run_dir: Optional[str] = None,
        wb_run_id: Optional[str] = None,
    ) -> None:
        # Resolve: explicit arg > constructor default > fallback
        root_wandb_dir = root_wandb_dir or self._root_wandb_dir
        wandb_run_dir = wandb_run_dir or self._wandb_run_dir

        print(f"Starting import from wandb directory: {os.path.abspath(root_wandb_dir)}")
        run_dirs = self._find_run_dirs(root_wandb_dir, wandb_run_dir)

        if not run_dirs:
            print("No runs found to import.")
            return

        print(f"Found {len(run_dirs)} run(s) to import.")

        for i, rd in enumerate(run_dirs):
            print(f"\n--- Processing run {i + 1}/{len(run_dirs)} ---")
            with safe.block(message=f"Failed to convert run in {rd}"):
                self._parse_run(rd, wb_run_id)

        print("\nAll runs processed.")

    def _find_run_dirs(self, root_wandb_dir: str, wandb_run_dir: Optional[str] = None) -> List[str]:
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
            print(f"No wandb run directories found in '{found_path}'.")
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
            print("Warning: pyecharts is required for table conversion. Skipping table.")
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
            raise FileNotFoundError(f"No .wandb file found in {run_dir}")

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

            print(f"Converting run: {run_metadata['name']} (ID: {run_metadata['id']})")
            swanlab.login(relogin=True)
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

            # Load wandb-metadata.json from wandb run files (non-critical)
            with safe.block(message="Failed to load wandb metadata", write_to_tty=False):
                metadata_path = os.path.join(files_root_dir, "wandb-metadata.json")
                with open(metadata_path, "rb") as f:
                    wandb_metadata = json_loads(f.read())
                run_config.update({"wandb_metadata": wandb_metadata})

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
                config_update = proto_items_to_dict(run.config.update)
                run_config.update(config_update)
                if swanlab_run is not None:
                    swanlab_run.config.update(config_update)

            elif record_type == "config":
                config_update = proto_items_to_dict(record_pb.config.update)
                run_config.update(config_update)
                if swanlab_run is not None:
                    swanlab_run.config.update(config_update)

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
            print(f"Finished converting run: {run_metadata['name']}")
            swanlab_run.finish()
        else:
            with safe.block(message=f"Failed to initialize run for {run_dir}"):
                initialize_swanlab_run()
                if swanlab_run:
                    print(f"Warning: Run in {run_dir} has no metrics, config was saved.")
                    swanlab_run.finish()

        del ds
        gc.collect()

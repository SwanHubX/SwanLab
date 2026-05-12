from typing import Optional

import click

from swanlab import vendor
from swanlab.cli.converter.base import BaseConverter


class WandbConverter(BaseConverter):
    """Convert W&B cloud runs to SwanLab via the W&B Public API."""

    def run(  # type: ignore[override]
        self,
        wb_project: str,
        wb_entity: str,
        wb_run_id: Optional[str] = None,
    ) -> None:
        wandb = vendor.wandb
        import swanlab

        client = wandb.Api()

        if wb_run_id:
            wb_runs = [client.run(f"{wb_entity}/{wb_project}/{wb_run_id}")]
        else:
            wb_runs = list(client.runs(f"{wb_entity}/{wb_project}"))

        if not wb_runs:
            click.echo(f"No W&B runs found in {wb_entity}/{wb_project}.")
            return

        click.echo(f"Found {len(wb_runs)} W&B run(s) to convert.")

        for i, wb_run in enumerate(wb_runs):
            click.echo(f"[{i + 1}/{len(wb_runs)}] Converting W&B run: {wb_run.name or wb_run.id}")
            swanlab_run = swanlab.init(
                project=self.project or wb_project,
                workspace=self.workspace,
                mode=self.mode,  # type: ignore[arg-type]
                log_dir=self.log_dir,
                name=wb_run.name,
                description=wb_run.notes,
                tags=wb_run.tags,
                group=wb_run.group,
                job_type=wb_run.job_type,
                config=wb_run.config,
                reinit=True,
                **({"resume": True, "id": wb_run.id} if (self.resume and wb_run_id) else {}),
            )

            for row in wb_run.scan_history():
                log_data = {}
                step = None
                for key, value in row.items():
                    if value is None:
                        continue
                    if key == "_step":
                        if isinstance(value, (int, float)):
                            step = int(value)
                        continue
                    if key.startswith("_"):
                        continue
                    if isinstance(value, dict):
                        continue
                    log_data[key] = value
                if log_data:
                    swanlab_run.log(log_data, step=step)

            swanlab_run.finish()
            click.echo(f"  ✓ Finished converting run: {wb_run.name or wb_run.id}")

        click.echo(f"All {len(wb_runs)} W&B run(s) converted successfully.")

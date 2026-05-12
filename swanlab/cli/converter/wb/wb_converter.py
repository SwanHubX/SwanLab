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

            swan_run = swanlab.init(
                project=self.project or wb_project,
                workspace=self.workspace,
                mode=self.mode,  # type: ignore[arg-type]
                log_dir=self.log_dir,
                name=wb_run.name,
                description=wb_run.notes,
                tags=wb_run.tags,
                group=wb_run.group,
                job_type=wb_run.job_type,
                reinit=True,
                **({"resume": True, "id": wb_run.id} if (self.resume and wb_run_id) else {}),
            )

            if wb_run.config:
                swan_run.config.update(wb_run.config)

            for row in wb_run.scan_history():
                log_data = {}
                for key, value in row.items():
                    if value is None or isinstance(value, dict):
                        continue
                    log_data[key] = value
                if log_data:
                    step = row.get("_step")
                    swan_run.log(log_data, step=step if isinstance(step, int) else None)

            swan_run.finish()
            click.echo(f"  ✓ Finished converting run: {wb_run.name or wb_run.id}")

        click.echo(f"All {len(wb_runs)} W&B run(s) converted successfully.")

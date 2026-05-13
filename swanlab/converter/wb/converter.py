from typing import Optional

from swanlab import vendor
from swanlab.converter.base import BaseConverter


class WandbConverter(BaseConverter):
    """Convert W&B cloud runs to SwanLab via the W&B Public API.

    Example::

        wb = WandbConverter(project="my-swanlab-proj")
        wb.run(wb_project="wandb-proj", wb_entity="my-team")
    """

    def __init__(
        self,
        project: Optional[str] = None,
        workspace: Optional[str] = None,
        mode: str = "online",
        log_dir: Optional[str] = None,
        logdir: Optional[str] = None,
        tags: Optional[list] = None,
        resume: bool = False,
        wb_project: Optional[str] = None,
        wb_entity: Optional[str] = None,
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
        # W&B defaults — can be overridden at run() time
        self._wb_project = wb_project
        self._wb_entity = wb_entity

    def run(  # type: ignore[override]
        self,
        wb_project: Optional[str] = None,
        wb_entity: Optional[str] = None,
        wb_run_id: Optional[str] = None,
    ) -> None:
        # Resolve: explicit arg > constructor default > error
        wb_project = wb_project or self._wb_project
        wb_entity = wb_entity or self._wb_entity
        if not wb_project:
            raise ValueError("wb_project is required (pass to WandbConverter() or run())")
        if not wb_entity:
            raise ValueError("wb_entity is required (pass to WandbConverter() or run())")
        wandb = vendor.wandb
        import swanlab

        client = wandb.Api()

        if wb_run_id:
            wb_runs = [client.run(f"{wb_entity}/{wb_project}/{wb_run_id}")]
        else:
            wb_runs = list(client.runs(f"{wb_entity}/{wb_project}"))

        if not wb_runs:
            raise ValueError(f"No W&B runs found in {wb_entity}/{wb_project}")

        print(f"Found {len(wb_runs)} W&B run(s) to convert.")

        for i, wb_run in enumerate(wb_runs):
            print(f"[{i + 1}/{len(wb_runs)}] Converting W&B run: {wb_run.name or wb_run.id}")
            swanlab.login(relogin=True)
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
            print(f"  ✓ Finished converting run: {wb_run.name or wb_run.id}")

        print(f"All {len(wb_runs)} W&B run(s) converted successfully.")

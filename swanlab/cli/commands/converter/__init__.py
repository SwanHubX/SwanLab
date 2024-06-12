from .tfb import TFBConverter
from .wb import WandbConverter
import click


@click.command()
@click.option(
    "--type",
    "-t",
    default="tensorboard",
    type=click.Choice(["tensorboard", "wandb"]),
    help="The type of the experiment tracking tool you want to convert to.",
)
@click.option(
    "--project",
    "-p",
    default=None,
    type=str,
    help="SwanLab project name.",
)
@click.option(
    "--workspace",
    "-w",
    default=None,
    type=str,
    help="swanlab.init workspace parameter.",
)
@click.option(
    "--cloud",
    default=True,
    type=bool,
    help="swanlab.init cloud parameter.",
)
@click.option(
    "--logdir",
    "-l",
    type=str,
    help="The directory where the swanlab log files are stored.",
)
@click.option(
    "--tb_logdir",
    type=str,
    help="The directory where the tensorboard log files are stored.",
)
@click.option(
    "--wb-project",
    type=str,
    help="The project name of the wandb runs.",
)
@click.option(
    "--wb-entity",
    type=str,
    help="The entity name of the wandb runs.",
)
@click.option(
    "--wb-runid",
    type=str,
    help="The run_id of the wandb run.",
)
def convert(
        type: str,
        project: str,
        cloud: bool,
        workspace: str,
        logdir: str,
        tb_logdir: str,
        wb_project: str,
        wb_entity: str,
        wb_runid: str,
        **kwargs,
):
    """Convert the log files of other experiment tracking tools to SwanLab."""
    if type == "tensorboard":

        tfb_converter = TFBConverter(
            convert_dir=tb_logdir,
            project=project,
            workspace=workspace,
            cloud=cloud,
            logdir=logdir,
        )
        tfb_converter.run()

    elif type == "wandb":

        print(wb_project, wb_entity, wb_runid)

        wb_converter = WandbConverter(
            project=project,
            workspace=workspace,
            cloud=cloud,
            logdir=logdir,
        )
        wb_converter.run(
            wb_project=wb_project,
            wb_entity=wb_entity,
            wb_run_id=wb_runid,
        )

    else:
        raise ValueError("The type of the experiment tracking tool you want to convert to is not supported.")

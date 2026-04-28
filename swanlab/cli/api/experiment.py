import click
import orjson

from swanlab.api import Api
from swanlab.cli.api.helper import format_output, save_output


@click.group("run")
def experiment_cli():
    """Experiment(Run) management commands."""
    pass


@experiment_cli.command("info")
@click.argument("path", required=True)
@click.option(
    "--save",
    "-s",
    "name",
    is_flag=False,
    flag_value=".",
    default=None,
    help="Save output as JSON to current directory.",
)
def get_experiment(path: str, name):
    """Get Experiment(Run) info by path (username/project/run_id)."""
    api = Api()
    resp = api.run(path).wrapper()
    format_output(resp)
    if resp.ok and name is not None:
        save_output(orjson.dumps(resp.json(), option=orjson.OPT_INDENT_2), name=name)

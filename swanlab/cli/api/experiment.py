import click
import orjson

from swanlab.api.typings.common import ApiResponseType

from .helper import format_output


@click.group("run")
def experiment_cli():
    """Experiment(Run) management commands."""
    pass


@experiment_cli.command("info")
@click.argument("path", required=True)
def get_info(path: str):
    """Get run(experiment) info by path (username/project/run_id)."""
    from swanlab.api import Api

    try:
        api = Api()
        resp = api.run(path)
    except Exception as e:
        resp = ApiResponseType(ok=False, errmsg=str(e))
    click.echo(orjson.dumps(resp.to_json_dict(), option=orjson.OPT_INDENT_2).decode())

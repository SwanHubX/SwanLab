import click
import orjson

from swanlab.api.typings.common import ApiResponseType


@click.command("run")
@click.argument("path")
def get_run(path: str):
    """Get run(experiment) info by path (username/project/run_id)."""
    from swanlab.api import Api

    try:
        api = Api()
        resp = api.run(path)
    except Exception as e:
        resp = ApiResponseType(ok=False, errmsg=str(e))
    click.echo(orjson.dumps(resp.to_json_dict(), option=orjson.OPT_INDENT_2).decode())

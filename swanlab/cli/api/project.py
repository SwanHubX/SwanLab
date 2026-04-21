import click
import orjson

from swanlab.api.typings.common import ApiResponseType


@click.command("project")
@click.argument("path")
def get_project(path: str):
    """Get project info by path (username/project)."""
    from swanlab.api import Api

    try:
        api = Api()
        resp = api.project(path)
    except Exception as e:
        resp = ApiResponseType(ok=False, errmsg=str(e))
    click.echo(orjson.dumps(resp.to_json_dict(), option=orjson.OPT_INDENT_2).decode())

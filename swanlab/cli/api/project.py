import click
import orjson

from swanlab.api.typings.common import ApiResponseType


@click.group("project")
def project_cli():
    """Project management commands."""
    pass


@project_cli.command("info")
@click.argument("path", required=True)
def get_info(path: str):
    """Get project info by path (username/project)."""
    from swanlab.api import Api

    try:
        api = Api()
        resp = api.project(path)
    except Exception as e:
        resp = ApiResponseType(ok=False, errmsg=str(e))
    click.echo(orjson.dumps(resp.to_json_dict(), option=orjson.OPT_INDENT_2).decode())

import click
import orjson


@click.command("project")
@click.argument("path")
def get_project(path: str):
    """Get project info by path (username/project)."""
    from swanlab.api import Api

    api = Api()
    project = api.project(path)
    click.echo(orjson.dumps(project.to_dict(), option=orjson.OPT_INDENT_2).decode())

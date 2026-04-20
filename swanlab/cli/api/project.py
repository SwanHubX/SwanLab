import json

import click


@click.command("project")
@click.argument("path")
def get_project(path: str):
    """Get project info by path (username/project)."""
    from swanlab.api import Api

    api = Api()
    project = api.project(path)
    click.echo(json.dumps(project.to_dict(), indent=2, ensure_ascii=False))

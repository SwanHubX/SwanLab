import click


@click.command("project")
@click.argument("path")
def get_project(path: str):
    """Get project info by path (username/project)."""
    raise NotImplementedError("cli api project")

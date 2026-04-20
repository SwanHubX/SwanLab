import click


@click.command("run")
@click.argument("path")
def get_run(path: str):
    """Get run(experiment) info by path (username/project/run_id)."""
    raise NotImplementedError("cli api run")

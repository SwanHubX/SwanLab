import click
import orjson

from swanlab.sdk.internal.pkg import safe


@click.command("run")
@click.argument("path")
def get_run(path: str):
    """Get run(experiment) info by path (username/project/run_id)."""
    from swanlab.api import Api

    api = Api()
    with safe.block(message="Failed to get run info"):
        run = api.run(path)
        click.echo(orjson.dumps(run.to_dict(), option=orjson.OPT_INDENT_2).decode())

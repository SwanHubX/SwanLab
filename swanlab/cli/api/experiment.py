import json

import click


@click.command("run")
@click.argument("path")
def get_run(path: str):
    """Get run(experiment) info by path (username/project/run_id)."""
    from swanlab.api import Api

    api = Api()
    run = api.run(path)
    click.echo(json.dumps(run.to_dict(), indent=2, ensure_ascii=False))

import click
import orjson

from swanlab.api.typings.common import ApiResponseType

from .helper import format_output


@click.group("selfhosted")
def selfhosted_cli():
    """Self-hosted deployment management commands."""
    pass


@selfhosted_cli.command("info")
@click.option("-u", "--username", default=None, help="Target username (defaults to current user).")
def get_info(username):
    """Get self-hosted info for a user."""
    from swanlab.api import Api

    try:
        api = Api()
        resp = api.user(username=username)
    except Exception as e:
        resp = ApiResponseType(ok=False, errmsg=str(e))
    format_output(resp)


@selfhosted_cli.command("keys")
@click.option("-u", "--username", default=None, help="Target username (defaults to current user).")
def get_keys(username):
    """List API keys for a user (only available for the logged-in user)."""
    from swanlab.api import Api

    try:
        api = Api()
        user_resp = api.user(username=username)
        if not user_resp.ok:
            format_output(user_resp)
            return
        user = user_resp.data
        keys = user.api_keys
        format_output(ApiResponseType(ok=True, data=keys))
    except Exception as e:
        format_output(ApiResponseType(ok=False, errmsg=str(e)))


@selfhosted_cli.command("create-user")
@click.argument("new_username", required=True)
@click.argument("password", required=True)
def create_user(new_username: str, password: str):
    """Create a new user in self-hosted deployment (root only)."""
    from swanlab.api import Api

    try:
        api = Api()
        user_resp = api.user()
        if not user_resp.ok:
            format_output(user_resp)
            return
        user = user_resp.data
        success = user.create(new_username, password)
        format_output(ApiResponseType(ok=success, data=None, errmsg="" if success else "Failed to create user"))
    except Exception as e:
        format_output(ApiResponseType(ok=False, errmsg=str(e)))

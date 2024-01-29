from swanlab.db import connect
from swanlab.server.app import app

connect()


def main():
    return app

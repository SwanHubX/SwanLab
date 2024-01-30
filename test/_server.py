# from swanlab.env import ROOT
# import os

# os.environ[ROOT] = "/Users/likang/Code/swanhub-dev/swanlab/test/experiments"
from swanlab.db import connect
from swanlab.server.app import app

connect()


def main():
    return app

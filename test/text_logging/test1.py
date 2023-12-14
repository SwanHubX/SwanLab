import logging

logging.basicConfig(format="%(asctime)s %(module)s %(message)s", datefmt="%m/%d/%Y %I:%M:%S %p")
logger = logging.getLogger(__name__)


def test_fun1():
    logger.warning("Test1")


test_fun1()

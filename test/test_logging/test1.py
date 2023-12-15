from swanlab.log import Swanlog
import logging

sl = Swanlog(__name__, log_level=logging.WARNING)
test_string = "test1"


def test1():
    sl.debug(test_string)
    sl.info(test_string)
    sl.warning(test_string)
    sl.error(test_string)
    sl.critical(test_string)

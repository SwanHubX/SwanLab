from swanlab.log import swanlog as sl
import logging

test_string = "test1"


def test1():
    # sl.setLevel("error")
    sl.debug(test_string)
    sl.info(test_string)
    sl.warning(test_string)
    sl.error(test_string)
    sl.critical(test_string)

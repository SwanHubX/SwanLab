from test1 import test1
from swanlab.log import swanlog as sl
import logging

# sl.setLevel("error")
# sl.setOutput()
sl.debug("Watch out!")
sl.info("I told you so")
sl.warning("I told you so")
sl.error("I told you so")
sl.critical("I told you so")
test1()

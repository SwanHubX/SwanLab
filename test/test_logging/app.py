from test1 import test1
from swanlab.log import Swanlog
import logging

sl = Swanlog(__name__)


sl.debug("Watch out!")
sl.info("I told you so")
sl.warning("I told you so")
sl.error("I told you so")
sl.critical("I told you so")
test1()

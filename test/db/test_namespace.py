from swanlab.db import Namespace
from playhouse.shortcuts import model_to_dict

Namespace.create_namespace("123", experiment_id=1)

# print(model_to_dict(Namespace.select().first()))

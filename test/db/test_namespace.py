from swanlab.db import Namespace
from swanlab.utils.time import create_time


# Namespace.create_namespace("123", experiment_id=1)

# Namespace.delete().where(Namespace.id == 1).execute()


Namespace.create(
    name="name",
    index=1,
    experiment_id=1,
    project_id=None,
)

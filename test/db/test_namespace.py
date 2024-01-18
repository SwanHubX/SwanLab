from swanlab.db import Namespace
from swanlab.utils.time import create_time


# Namespace.create_namespace("123", experiment_id=1)

# Namespace.delete().where(Namespace.id == 1).execute()


Namespace.get_or_create(
    name="name",
    index=1,
    create_time=create_time(),
    update_time=create_time(),
)

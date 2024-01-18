from swanlab.db import Experiment

exp = Experiment.create("test", run_id="test", description="test", more={"test": "test"})

print(exp.create_time)

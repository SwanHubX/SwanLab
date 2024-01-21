from swanlab.db import Experiment

exp = Experiment.create("test-1", run_id="test-1", description="test", more={"test": "test"})

a = exp.__dict__()


print(exp.create_time)

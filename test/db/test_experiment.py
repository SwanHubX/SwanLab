from swanlab.db import Experiment

print(Experiment.create_experiment(run_id="1234").run_id)

# print(Experiment.check_run_id("123123123"))
# print(Experiment.get_experiment_by_runid("1"))
# print(Experiment.get_experiment("1"))

# print(Experiment.delete_experiment("1"))

# Experiment.update_info(2, "experiment-2222", "none")

# Experiment.update_updatetime(2)

# Experiment.update_status(2, 1)

# print(Experiment.get_experiments())

# print(Experiment.get_tags(1))

# print(Experiment.get_experiments())

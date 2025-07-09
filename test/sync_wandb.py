import wandb
import random
import swanlab
import numpy as np

swanlab.sync_wandb(wandb_run=False)

wandb.init(
  project="test",
  config={"a": 1, "b": 2},
  name="test",
  notes="test_wandb_sync",
  )

wandb.config.update({"c": 3, "d": 4})
print(swanlab.config.get("c"))

wandb.log({"im": wandb.Image(np.random.randint(0, 255, (3, 100, 100)))})

epochs = 10
offset = random.random() / 5
for epoch in range(2, epochs):
  acc = 1 - 2 ** -epoch - random.random() / epoch - offset
  loss = 2 ** -epoch + random.random() / epoch + offset

  wandb.log({"acc": acc, "loss": loss})
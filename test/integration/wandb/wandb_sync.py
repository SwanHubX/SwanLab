import wandb
import random
import swanlab

swanlab.sync_wandb()

wandb.init(
  project="test",
  config={"a": 1, "b": 2},
  name="test",
  )

epochs = 10
offset = random.random() / 5
for epoch in range(2, epochs):
  acc = 1 - 2 ** -epoch - random.random() / epoch - offset
  loss = 2 ** -epoch + random.random() / epoch + offset

  wandb.log({"acc": acc, "loss": loss})
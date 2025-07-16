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
  tags=["test", "wandb_sync"],
  )

wandb.config.update({"c": 3, "d": 4})
print(swanlab.config.get("c"))

wb_image = wandb.Image(np.random.randint(low=0, high=256, size=(100, 100, 3), dtype=np.uint8), caption="test_image")
wb_image_list = [wb_image] * 2

wandb.log({"im": wb_image_list}, step=10)

epochs = 10
offset = random.random() / 5
for epoch in range(2, epochs):
  acc = 1 - 2 ** -epoch - random.random() / epoch - offset
  loss = 2 ** -epoch + random.random() / epoch + offset

  wandb.log({"acc": acc, "loss": loss})
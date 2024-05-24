from tutils import open_dev_mode
import swanlab

swanlab.login(open_dev_mode())

import os
from lightning import Trainer
from swanlab.integration.pytorch_lightning import SwanLabLogger
from tutils import open_dev_mode
from lightning_base import BoringModel, RandomDataset
from torch.utils.data import DataLoader


def main():
    print("User process PID:", os.getpid())

    # Set up data
    num_samples = 100000
    train = DataLoader(RandomDataset(32, num_samples), batch_size=32)
    val = DataLoader(RandomDataset(32, num_samples), batch_size=32)
    test = DataLoader(RandomDataset(32, num_samples), batch_size=32)

    model = BoringModel()

    swanlab_logger = SwanLabLogger(
        project="LightningTest",
        config={"num_samples": num_samples},
    )

    # Initialize a trainer
    trainer = Trainer(
        max_epochs=2,
        logger=swanlab_logger,
    )

    # Train the model
    trainer.fit(model, train, val)
    trainer.test(dataloaders=test)


if __name__ == "__main__":
    # swanlab.login(open_dev_mode())
    main()

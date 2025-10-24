try:
    import pytorch_lightning as pl
    from pytorch_lightning import demos # Need Import Again
except ImportError:
    pl = None

import swanlab
from swanlab.data.run import SwanLabRunState
from swanlab.integration.pytorch_lightning import SwanLabLogger


def pl_strategy_for_swanlab(raw_fn, swanlab_logger: SwanLabLogger):
    def wrapper(exception: BaseException):
        if isinstance(exception, KeyboardInterrupt):
            swanlab_logger.experiment.finish(SwanLabRunState.CRASHED, "KeyboardInterrupt by user", interrupt=True)
        raw_fn(exception)
    return wrapper


def test_for_pytorch_lighting():
    # Also See SwanLab/test/integration/lightning
    if pl is None:
        Warning("Need PyTorch Lightning To Test")
        return

    model = pl.demos.boring_classes.BoringModel()
    datamodule = pl.demos.boring_classes.BoringDataModule()
    logger = SwanLabLogger(project="test-interrupt")

    # log_every_n_steps != None for enable Logger
    trainer = pl.Trainer(max_epochs=10, logger=logger, log_every_n_steps=1)

    # 加上可使pl被KeyboardInterrupt中断时，Web上的状态为 "Ctrl+C" 而不是 "中断"
    # Add support so that when pl is interrupted by a KeyboardInterrupt, the status displayed on the web is "Ctrl+C" instead of "Crashed"
    # trainer.strategy.on_exception = pl_strategy_for_swanlab(trainer.strategy.on_exception, logger) # 在pl的策略补充即可

    trainer.fit(model, datamodule)

def test_for_normal_interrupt():
    from tqdm import tqdm
    # import time

    swanlab.init(project="test-interrupt")

    for i in tqdm(range(10000)):
        # time.sleep(0.01)

        for _ in range(100):
            for _ in range(100):
                for _ in range(100):
                    pass

        if i % 10 == 0:
            swanlab.log({"progress": i}, print_to_console=True)


if __name__ == '__main__':
    test_for_pytorch_lighting()

    # test_for_normal_interrupt()


import swanlab
import time
import random
import numpy as np

epochs = 50
lr = 0.01
offset = random.random() / 5

swanlab.init(
    log_level="debug",
    config={
        "epochs": epochs,
        "learning_rate": lr,
        "test": 1,
        "debug": "这是一串" + "很长" * 100 + "的字符串",
        "verbose": 1,
    },
    logggings=True,
)
for epoch in range(2, epochs):
    if epoch % 10 == 0:
        # 测试audio
        sample_rate = 44100
        test_audio_arr = np.random.randn(2, 100000)
        swanlab.log(
            {
                "test/audio": [swanlab.Audio(test_audio_arr, sample_rate, caption="test")] * (epoch // 10),
            },
            step=epoch,
        )
        # 测试image
        test_image = np.random.randint(0, 255, (100, 100, 3))
        swanlab.log(
            {
                "test/image": swanlab.Image(test_image, caption="test"),
            },
            step=epoch,
        )
    acc = 1 - 2**-epoch - random.random() / epoch - offset
    loss = 2**-epoch + random.random() / epoch + offset
    loss2 = 3**-epoch + random.random() / epoch + offset * 3
    print(f"epoch={epoch}, accuracy={acc}, loss={loss}")
    swanlab.log({"t/accuracy": acc, "loss": loss, "loss2": loss2})
    time.sleep(0.2)

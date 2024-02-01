import swanlab
import time
import random
import numpy as np

epochs = 100
lr = 0.01
offset = random.random() / 5

swanlab.init(
    log_level="info",
    config={
        "epochs": epochs,
        "learning_rate": lr,
        "test": 1,
        "verbose": 1,
    },
)

for epoch in range(2, epochs):
    if epoch % 10 == 0:
        # audio
        sample_rate = 44100
        test_audio_arr = np.random.randn(2, 100000)
        swanlab.log({"WhiteNoise": swanlab.Audio(test_audio_arr, sample_rate, caption=123)}, step=epoch)
        swanlab.log({"Music": swanlab.Audio("/Users/zeyilin/Desktop/Coding/swanlab/test/assets/慢冷.mp3")}, step=epoch)
    acc = 1 - 2**-epoch - random.random() / epoch - offset
    loss = 2**-epoch + random.random() / epoch + offset
    loss2 = 3**-epoch + random.random() / epoch + offset * 3
    print(f"epoch={epoch}, accuracy={acc}, loss={loss}")
    swanlab.log({"accuracy": acc, "loss": loss, "loss2": loss2})
    time.sleep(0.2)

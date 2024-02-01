import swanlab
import time

epochs = 10

swanlab.init(
    experiment_name="TestTextChart",
    description="测试文本图表",
    log_level="info",
    config={
        "epochs": epochs,
    },
)

for epoch in range(1, epochs):
    if epoch == 1:
        swanlab.log({"text1": swanlab.Text("测试文本001", caption="测试文本001")})
    elif epoch == 2:
        swanlab.log({"text2": swanlab.Text("测试文本002")})
    else:
        swanlab.log({"text3": swanlab.Text("测试文本003")})
    time.sleep(0.2)

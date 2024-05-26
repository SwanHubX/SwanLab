from tutils import open_dev_mode
import swanlab

swanlab.login(open_dev_mode())

from ultralytics import YOLO
from swanlab.integration.ultralytics import add_swanlab_callback


def main():
    model = YOLO("yolov8n.pt")
    add_swanlab_callback(model)
    model.train(data="coco128.yaml", epochs=2, imgsz=640, batch=64)


if __name__ == "__main__":
    main()

from ultralytics import YOLO
from swanlab.integration.ultralytics import add_swanlab_callback


def main():
    model = YOLO("yolov8n-cls.pt")
    add_swanlab_callback(model)
    model.train(data="mnist160", epochs=1, imgsz=64)


if __name__ == "__main__":
    main()

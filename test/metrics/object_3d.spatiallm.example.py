"""
Example SwanLab code for visualizing point clouds and bounding boxes from the SpatialLM-Testset.

This code demonstrates how to load .ply files and scene description files from the
[SpatialLM-Testset](https://huggingface.co/datasets/manycore-research/SpatialLM-Testset) dataset,
process the data, and visualize it using SwanLab's 3D object logging capabilities.

Before running this code, make sure to:
1.  Install the required packages: `pip install numpy plyfile swanlab`
2.  Download the necessary .ply and .txt scene files from the SpatialLM-Testset dataset
    at [https://huggingface.co/datasets/manycore-research/SpatialLM-Testset](https://huggingface.co/datasets/manycore-research/SpatialLM-Testset).
    Specifically, you'll need `scene0000_00.ply` and `scene0000_00.txt` (or similar files).
3.  Update the `ply_file` and `scene_file` variables with the correct paths to your downloaded files.
"""

import re
from typing import List, Tuple

import numpy as np
import swanlab
from plyfile import PlyData


ply_file = "scene0000_00.ply"  # 替换为你的 .ply 文件路径
scene_file = "scene0000_00.txt"  # 替换为你的场景描述文件路径

swanlab.init(project="SpatialLM-Example-PointCloud", public=True)


def load_ply_data(ply_file_path: str) -> np.ndarray | None:
    """加载 .ply 文件并转换为 NumPy 数组。

    Args:
        ply_file_path (str): .ply 文件路径。

    Returns:
        np.ndarray: 包含顶点坐标 (x, y, z) 和颜色信息 (red, green, blue) 的 NumPy 数组，shape 为 (N, 6)，如果加载失败则返回 None。
    """
    try:
        ply_data = PlyData.read(ply_file_path)
        vertices = ply_data["vertex"].data
        vertices_array = np.vstack([vertices["x"], vertices["y"], vertices["z"]]).T
        colors_array = np.vstack(
            [vertices["red"], vertices["green"], vertices["blue"]]
        ).T
        return np.concatenate((vertices_array, colors_array), axis=1)
    except Exception as e:
        print(f"Error loading or processing .ply file: {e}")
        return None


def parse_bbox_line(line: str) -> Tuple[int, str, List[float]] | None:
    """从文本行解析识别框信息。

    Args:
        line (str): 包含识别框信息的文本行，格式为 "bbox_id=object_type(param1, param2, ...)"。

    Returns:
        Tuple[int, str, List[float]]: 包含边界框 ID、对象类型和参数的元组，如果解析失败则返回 None。
    """
    match = re.match(r"bbox_(\d+)=(.*)\((.*)\)", line)
    if not match:
        return None

    bbox_id = int(match.group(1))
    object_type = match.group(2).strip()
    params_str = match.group(3)

    try:
        params = [float(param_str.strip()) for param_str in params_str.split(",")]
        return bbox_id, object_type, params
    except ValueError as e:
        print(f"Warning: Could not convert parameter to float: {e}")
        return None


def create_bounding_box(params: List[float], object_type: str) -> dict:
    """根据参数创建边界框字典。

    Args:
        params (List[float]): 边界框参数，包括中心点坐标 (x, y, z)、绕 Z 轴的旋转角度 (弧度)、宽度、深度和高度。
        object_type (str): 对象类型。

    Returns:
        dict: 包含边界框信息的字典，包括颜色、角点、标签和置信度。
    """
    center_x, center_y, center_z = params[:3]  # 中心点坐标
    rotation_z = params[3]  # 绕 Z 轴的旋转角度 (弧度)
    width, depth, height = params[4:7]  # 宽度、深度、高度

    # 计算旋转矩阵
    rotation_matrix = np.array(
        [
            [np.cos(rotation_z), -np.sin(rotation_z), 0],
            [np.sin(rotation_z), np.cos(rotation_z), 0],
            [0, 0, 1],
        ]
    )

    # 定义未旋转的角点
    corners_unrotated = np.array(
        [
            [-width / 2, -depth / 2, -height / 2],
            [width / 2, -depth / 2, -height / 2],
            [width / 2, depth / 2, -height / 2],
            [-width / 2, depth / 2, -height / 2],
            [-width / 2, -depth / 2, height / 2],
            [width / 2, -depth / 2, height / 2],
            [width / 2, depth / 2, height / 2],
            [-width / 2, depth / 2, height / 2],
        ]
    )

    # 旋转角点
    corners_rotated = np.dot(corners_unrotated, rotation_matrix.T)

    # 平移角点到最终位置
    corners = corners_rotated + np.array([center_x, center_y, center_z])

    # 随机生成一个颜色
    color = np.random.randint(0, 256, size=3).tolist()

    return {
        "color": color,
        "corners": corners.tolist(),
        "label": object_type,
        "score": 0.9,  # 假设一个置信度
    }


def load_bounding_boxes(scene_file_path: str) -> List[dict] | None:
    """加载场景描述文件并提取识别框信息。

    Args:
        scene_file_path (str): 场景描述文件路径。

    Returns:
        List[dict]: 包含边界框信息的字典列表，如果加载失败则返回 None。
    """
    try:
        with open(scene_file_path, "r") as f:
            bounding_boxes = [
                create_bounding_box(params, object_type)
                for line in f
                if line.startswith("bbox_")
                and (bbox_info := parse_bbox_line(line))
                and (bbox_id, object_type, params := bbox_info)  # noqa: F841
            ]  # noqa: F821
        return bounding_boxes
    except FileNotFoundError:
        print(f"Error: Scene file not found at {scene_file_path}")
        return None
    except Exception as e:
        print(f"Error reading or parsing scene file: {e}")
        return None


def main(ply_file_path: str, scene_file_path: str) -> dict | None:
    """主函数，加载 .ply 文件和识别框信息。

    Args:
        ply_file_path (str): .ply 文件路径。
        scene_file_path (str): 场景描述文件路径。

    Returns:
        dict: 包含点云数据和边界框信息的字典，如果加载失败则返回 None。
    """
    points = load_ply_data(ply_file_path)
    boxes = load_bounding_boxes(scene_file_path)

    if points is None or boxes is None:
        print("Failed to load data.")
        return None

    print(f"Successfully loaded {points.shape[0]} points from {ply_file_path}")
    print(f"Successfully loaded {len(boxes)} bounding boxes from {scene_file_path}")

    return {"points": points, "boxes": boxes}


if __name__ == "__main__":
    data = main(ply_file, scene_file)

    if data is None:
        print("Data loading failed.")
        exit()

    print("Data loading and processing complete.")
    print("Points shape:", data["points"].shape)
    print("Number of bounding boxes:", len(data["boxes"]))

    pc = swanlab.Object3D({"points": data["points"], "boxes": data["boxes"]})
    swanlab.log({"example": pc})

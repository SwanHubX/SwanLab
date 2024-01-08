import torch
import torchvision


def format_size(size):
    # Define thresholds for units
    K, M, B = 1e3, 1e6, 1e9
    if size == 0:
        return "0"
    elif size < M:
        return f"{size / K:.1f}K"
    elif size < B:
        return f"{size / M:.1f}M"
    else:
        return f"{size / B:.1f}B"


def get_pytorch_model_info(model: torch.nn.Module) -> (str, list):
    """
    输入一个PyTorch Model对象，返回模型的总参数量（格式化为易读格式）以及每一层的名称、尺寸、精度、参数量、是否可训练。

    :param model: PyTorch Model
    :return: (总参数量dict[包括总参数量、可训练的总参数量、不可被训练的总参数量], 参数列表[包括每层的名称、尺寸和数据类型])
    """

    params_list = []
    total_params = 0
    total_params_non_trainable = 0

    for name, param in model.named_parameters():
        params_cout = param.numel()
        trainable = param.requires_grad
        params_list.append(
            [
                {
                    "tensor": name,
                    "shape": str(list(param.size())),
                    "prediction": str(param.dtype).split(".")[-1],
                    "params_count": str(params_cout),
                    "trainable": str(trainable),
                }
            ]
        )
        total_params += params_cout
        if not trainable:
            total_params_non_trainable += params_cout

    total_params_trainable = total_params - total_params_non_trainable

    total_params_info = {
        "total_params": format_size(total_params),
        "total_params_trainable": format_size(total_params_trainable),
        "total_params_non_trainable": format_size(total_params_non_trainable),
    }

    return total_params_info, params_list


if __name__ == "__main__":
    model = torchvision.models.resnet50(weights=None)
    in_features = model.fc.in_features
    model.fc = torch.nn.Linear(in_features, 2)

    formatted_params, params_list = get_pytorch_model_info(model)

    print(formatted_params)
    for item in params_list:
        print(item)

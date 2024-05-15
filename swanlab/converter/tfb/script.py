from utils import find_tfevents, get_tf_events_tags_type, get_tf_events_tags_data
import os

# TFEVENT日志文件路径
LOGDIR = "/Users/zeyilin/Desktop/Coding/wandb/Transfer-TF-WANDB/logs"
# 测试空路径
# LOGDIR = "/Users/zeyilin/Desktop/Coding/wandb/Transfer-TF-WANDB/script"

# 找到所有TFEvent文件
path_dict = find_tfevents(LOGDIR)
if path_dict:
    print("找到的文件路径路径字典:", path_dict)

# 打印字典，了解结构
"""
Dir: a; Path:['/Users/zeyilin/Desktop/Coding/wandb/Transfer-TF-WANDB/logs/a/events.out.tfevents.1715134571.qingdeMBP']
Dir: b; Path:['/Users/zeyilin/Desktop/Coding/wandb/Transfer-TF-WANDB/logs/b/events.out.tfevents.1715134202.qingdeMBP']
"""
for dir, paths in path_dict.items():
    print("=====>>>> Dir: {};\tPath:{}".format(dir, paths))
    for path in paths:
        filename = os.path.basename(path)
        # 获取所有的tag及其类型
        """
        {'training_loss': 'scalar', 'fake_image': 'image', 'fake audio': 'audio', 'fake text/text_summary': 'text'}
        """
        tags = get_tf_events_tags_type(path)
        print("=====>>>> Tags: {} \n".format(tags))

        # 根据tag提取数据
        """
        step, value, wall_time
        {
            'training_loss': [
                (0, 0.0, '2024-05-08 10:16:11'),
                (1, 0.019999999552965164, '2024-05-08 10:16:11'),
                (2, 0.03999999910593033, '2024-05-08 10:16:11')
                ]
        }
        """
        if tags:
            import swanlab

            run = swanlab.init(
                project="wandb-transfer-tf",
                experiment_name=f"{dir}/{filename}",
                config={"tfevent_path": path},
            )

            data_by_tag = get_tf_events_tags_data(path, tags)
            # print("=====>>>> Data by tag: {} \n".format(data_by_tag))

            if data_by_tag:
                # 打印并转换数据到SwanLab
                for tag, data in data_by_tag.items():
                    print(f"Data for {tag} ({tags[tag]}):")
                    for step, value, time in data:
                        print(f"Step: {step}, Time: {time}, Value: {value}")
                        # 如果是标量
                        if tags[tag] == "scalar":
                            swanlab.log({tag: value}, step=step)
                        # 如果是图片
                        elif tags[tag] == "image":
                            swanlab.log({tag: swanlab.Image(value)}, step=step)
                        # 如果是音频
                        elif tags[tag] == "audio":
                            swanlab.log(
                                {tag: swanlab.Audio(value[0], sample_rate=value[1])},
                                step=step,
                            )
                        # 如果是文本
                        elif tags[tag] == "text":
                            swanlab.log({tag: swanlab.Text(value)}, step=step)

            run.finish()

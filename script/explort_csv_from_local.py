################ 使用方法
# 用于从本地文件中提取指标数据，并保存为csv
# 使用方法
# ```python explort_csv_from_local.py --exp <实验名称> --indicator <指标名称>```
# 参数说明：
#     --root swanlab本地日志的路径，默认为./swanlog
#     --exp 实验名称，由于swanlab会默认在实验名后增加时间戳，需要完整填写
#     --indicator 需要打印的指标名称，默认为loss
#     --output 导出的csv路径，默认为.output.csv

################


import os
from pathlib import Path
import sqlite3
import json
import csv
import argparse


parser = argparse.ArgumentParser(description="Process some experiment log data.")
parser.add_argument(
    "--exp",
    type=str,
    required=True,
    help="Name of the experiment, including the timestamp suffix.",
)
parser.add_argument(
    "--root",
    type=str,
    default="swanlog",
    help="Path to the root directory where logs are saved.",
)
parser.add_argument(
    "--indicator",
    type=str,
    default="loss",
    help="Name of the indicator to be processed.",
)
parser.add_argument(
    "--output",
    type=str,
    default="output.csv",
    help="Path to save the output CSV file.",
)
args = parser.parse_args()
print(f"log root: {args.root}")
print(f"exp name: {args.exp}")
print(f"indicator: {args.indicator}")
print(f"output_path: {args.output}")
swanlog_db_path = os.path.join(args.root, "runs.swanlab")
conn = sqlite3.connect(swanlog_db_path)
cursor = conn.cursor()
# cursor.execute(f"SELECT id FROM experiment WHERE name = ?", ("{args.exp}",))
try:
    cursor.execute(f"""SELECT id FROM experiment WHERE name =="{args.exp}";""")
    exp_id = cursor.fetchone()[0]
except Exception as e:
    print("未查找到实验，注意swanlab会默认给实验名后面加时间后缀，需要完整填写")
    print(e)
    raise
cursor.execute(f"""SELECT run_id FROM experiment WHERE name =="{args.exp}";""")
run_id = cursor.fetchone()[0]
# cursor.execute(
#     f'''SELECT folder FROM tag WHERE name =="{args.indicator}" and experiment_id=="{exp_id}"'''
# )
cursor.execute(f"""SELECT folder FROM tag WHERE name = "{args.indicator}" AND experiment_id = "{exp_id}";""")
folder_id = cursor.fetchone()[0]
data_root = os.path.join(args.root, run_id, "logs", folder_id)
log_path_list = list(Path(data_root).glob("*.log"))
data_list = []
for fp in log_path_list:
    with open(fp, "r") as file:
        for line in file:
            data = json.loads(line)
            data_list.append(data)
with open(args.output, "w", newline="") as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=data_list[0].keys())
    writer.writeheader()
    writer.writerows(data_list)
print(f"Data has been saved to {args.output}")

# 实用脚本说明文档

本文件夹下用于存放一些临时性的过度脚本，或者一些实用的批处理方法

## explort_csv_from_local.py

用于从本地文件中提取指标数据，并保存为csv

**使用方法：**

```bash
python explort_csv_from_local.py --exp <实验名称> --indicator <指标名称>
```

**参数说明：**

```bash
--root swanlab本地日志的路径，默认为./swanlog

--exp 实验名称，由于swanlab会默认在实验名后增加时间戳，需要完整填写

--indicator 需要打印的指标名称，默认为loss

--output 导出的csv路径，默认为.output.csv
```

## transfer_logfile_0.1.4.py

Swanlab update script for upgrading from v0.1.4 to v0.1.5, transitioning from file system to database.
This script depends on swanlab v0.1.5, so please upgrade to swanlab v0.1.5 and backup the original log data before using this script—even though this script will not delete the original log data.

The following is the usage method:

```bash
python transfer_logfile_0.1.4.py -i /path/to/swanlog_old -o /path/to/swanlog_new
```
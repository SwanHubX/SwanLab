#!/bin/bash
source .env
read -p "请输入 run_dir: " run_dir

while true; do
    swanlab sync "$run_dir" --id "auto"
    sleep 5
done
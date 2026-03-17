#!/bin/bash
# Description: 递归清除指定目录下的 __pycache__ 和 .pyc 文件
# 如果删除 __pycache__ 后导致其父目录为空，则连同由于删除产生的空父目录一起递归向上删除

# 目标目录，默认当前目录
TARGET_DIR="${1:-.}"

echo "开始清理 $TARGET_DIR 目录下的 Python 缓存文件..."

# 1. 寻找所有的 __pycache__ 文件夹，并循环处理
# 注意：find 结果通过 while 读取，处理带有空格等特殊符号的路径
find "$TARGET_DIR" -type d -name "__pycache__" | while read -r cache_dir; do
    # 获取 __pycache__ 所在的父目录
    parent_dir="$(dirname "$cache_dir")"
    
    # 删除该 __pycache__ 目录及里面的缓存文件
    rm -rf "$cache_dir"
    
    echo "已清理: $cache_dir"
    
    # 尝试递归往上删除空的父目录
    # rmdir 的特性是：只有在此目录没有任何文件/文件夹时才会删除并返回成功
    # 因此我们循环往上 rmdir，如果哪个父层级除了 __pycache__ 还有其它文件，自然会失败并跳出循环
    while [ "$parent_dir" != "." ] && [ "$parent_dir" != "/" ] && [ "$parent_dir" != "$TARGET_DIR" ]; do
        if rmdir "$parent_dir" 2>/dev/null; then
            echo "已删除原本仅含 __pycache__ 的空父目录: $parent_dir"
            # 指向更上层父目录继续检查
            parent_dir="$(dirname "$parent_dir")"
        else
            # 遇到了非空目录（包含了其它文件），终止向上删除
            break
        fi
    done
done

# 2. 清理散落在其它地方的独立 .pyc 后缀缓存文件
find "$TARGET_DIR" -type f -name "*.pyc" -delete 2>/dev/null

echo "清理完成！"

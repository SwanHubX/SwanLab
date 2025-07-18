"""
@author: cunyue
@file: transformers_fake_train.py
@time: 2025/7/18 16:33
@description: 模拟 transformers 训练过程，使用 Qwen2 模型进行简单的训练和评估
NOTE: 该脚本运行于项目根目录，在项目根目录下执行: `python test/integration/transformers/transformers_fake_train.py`
"""

import argparse
import os.path

import dotenv
import torch
from datasets import Dataset
from transformers import Qwen2Config, Qwen2ForCausalLM, TrainingArguments, Trainer

parser = argparse.ArgumentParser(description="Qwen2模型模拟训练脚本")
parser.add_argument("--num_samples", type=int, default=5000, help="生成的样本数量（默认为5000）")
parser.add_argument("--seq_length", type=int, default=32, help="每个序列的长度（默认为32）")
args = parser.parse_args()

dotenv.load_dotenv()

config = Qwen2Config(
    vocab_size=24,  # GPT-2 标准词汇表大小
    hidden_size=16,
    intermediate_size=32,
    num_hidden_layers=3,
    num_attention_heads=4,
    num_key_value_heads=4,
)

model = Qwen2ForCausalLM(config)


seq = torch.randint(0, 24, (args.num_samples, args.seq_length)).tolist()
ds = Dataset.from_dict({"sequence": seq})
ds = ds.train_test_split(0.2)


train_args = TrainingArguments(
    output_dir=os.path.join(os.path.dirname(__file__), "test_transformers_output"),
    max_steps=120,
    per_device_train_batch_size=2,
    report_to="swanlab",
    logging_steps=1,
    eval_strategy="steps",
    eval_steps=0.1,
    remove_unused_columns=False,
)


def collator(examples):
    data = [example["sequence"] for example in examples]
    data = torch.tensor(data, dtype=torch.int)
    return {
        "input_ids": data,
        "attention_mask": torch.ones_like(data),
        "labels": data,
    }


trainer = Trainer(model, train_args, data_collator=collator, train_dataset=ds["train"], eval_dataset=ds["test"])

trainer.train()

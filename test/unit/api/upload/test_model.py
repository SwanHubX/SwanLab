#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/6/10 20:08
@File: test_model.py
@IDE: pycharm
@Description:
    测试上传的模型类型
"""

import random
import time

from nanoid import generate

from swanlab.api.upload.model import FileModel


class TestFileModel:

    def test_create(self):
        """
        测试FileModel的create函数
        """
        # 生成一系列的FileModel
        file_models = []
        lr = None  # 最后一个requirements, str
        lm = None  # 最后一个metadata, dict
        lc = None  # 最后一个config, dict
        lo = None  # 最后一个conda str
        for i in range(10):
            lr = generate() if random.random() > 0.5 else lr
            lm = {generate(): generate()} if random.random() > 0.5 else lm
            lc = {generate(): generate()} if random.random() > 0.5 else lc
            lo = generate() if random.random() > 0.5 else lo
            file_models.append(FileModel(lr, lm, lc, lo))
            time.sleep(0.01)
        # 验证
        # file_models打乱顺序
        random.shuffle(file_models)
        file_model = FileModel.create(file_models)
        assert file_model.requirements == lr
        assert file_model.metadata == lm
        assert file_model.config == lc
        assert file_model.conda == lo

    def test_empty(self):
        assert FileModel().empty is True
        assert FileModel("1", None, None, None).empty is False
        assert FileModel(None, {}, None, None).empty is False

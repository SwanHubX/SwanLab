#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2025/5/1 17:30
@File: project.py
@IDE: pycharm
@Description:
    项目相关的开放API
"""

from swanlab.api.http import HTTP
from swanlab.api.openapi.base import ApiBase, ApiHTTP


class ProjectAPI(ApiBase):
    def __init__(self, http: ApiHTTP):
        super().__init__(http)


    def list_projects(self, username: str, detail: bool):
        """
        列出一个 workspace 下的所有项目

        Args:
            username (str): 工作空间名, 默认为用户个人空间
            detail (bool): 是否返回实验详细信息，默认传入 True
        """
        base_url = f"/project/{username}?detail={detail}"
        
        result = {
            'total': 0,
            'list': []
        }

        raw_resp = self.http.get(base_url)
        if not isinstance(raw_resp, dict):
            return result
        
        result['total'] = raw_resp.get('total', 0)
        
        page_count = raw_resp.get('pages', 1)
        # 分页获取 
        for page in range(1, page_count + 1):
            page_url = f"{base_url}&page={page}"
            page_resp = self.http.get(page_url)
            
            # 移除avatar字段
            for item in page_resp.get('list', []):
                if item.get('group') and 'avatar' in item['group']:
                    del item['group']['avatar']
                result['list'].append(item)
        
        return result
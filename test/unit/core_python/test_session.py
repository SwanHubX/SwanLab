"""
@author: cunyue
@file: test_session.py
@time: 2025/9/9 15:12
@description: $END$
"""

import responses
from responses import registries

from swanlab.core_python import create_session


@responses.activate(registry=registries.OrderedRegistry)
def test_retry():
    """
    测试重试机制
    """
    from swanlab.package import get_host_api

    url = get_host_api() + "/retry"

    [responses.add(responses.GET, url, body="Error", status=500) for _ in range(5)]
    responses.add(responses.GET, url, body="Success", status=200)
    s = create_session()
    resp = s.get(url)
    assert resp.text == "Success"
    assert len(responses.calls) == 6

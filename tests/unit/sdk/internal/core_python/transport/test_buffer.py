from swanlab.sdk.internal.core_python.transport.buffer import RecordBuffer


def test_record_buffer_has_no_records_property():
    """不应暴露内部 records 引用。"""
    buf = RecordBuffer()
    assert not hasattr(buf, "records")


def test_record_buffer_extend_dedups_by_num(make_scalar_record):
    """extend() 先收集后批量追加，并按 num 去重。"""
    buf = RecordBuffer()
    first = make_scalar_record(step=1)
    second = make_scalar_record(step=2)
    first.num = 1
    second.num = 1

    accepted = buf.extend([first, second])

    assert accepted == 1
    assert buf.drain() == [first]


def test_record_buffer_prepend_keeps_order_and_dedups(make_scalar_record):
    """prepend() 保持传入顺序，并跳过已存在 num。"""
    buf = RecordBuffer()
    existing = make_scalar_record(step=1)
    existing.num = 1
    buf.extend([existing])

    first = make_scalar_record(step=2)
    second = make_scalar_record(step=3)
    duplicate = make_scalar_record(step=4)
    first.num = 2
    second.num = 3
    duplicate.num = 1

    accepted = buf.prepend([first, second, duplicate])

    assert accepted == 2
    assert [record.num for record in buf.drain()] == [2, 3, 1]


def test_record_buffer_drain_clears_contents(make_scalar_record):
    """drain() 取出全部记录并清空自身。"""
    buf = RecordBuffer()
    record = make_scalar_record(step=1)
    record.num = 7
    buf.extend([record])

    drained = buf.drain()

    assert drained == [record]
    assert not buf
    assert len(buf) == 0

import pytest

from swanlab.sdk.internal.core_python.transport.buffer import RecordBuffer


def test_record_buffer_has_no_records_property():
    """不应暴露内部 records 引用。"""
    buf = RecordBuffer()
    assert not hasattr(buf, "records")


def test_record_buffer_extend_replaces_by_num(make_scalar_record):
    """extend() 用后到记录覆盖同编号内容。"""
    buf = RecordBuffer()
    first = make_scalar_record(step=1)
    second = make_scalar_record(step=2)
    first.num = 1
    second.num = 1

    accepted = buf.extend([first, second])

    assert accepted == 1
    assert buf.drain() == [second]


def test_record_buffer_prepend_keeps_order_and_dedups(make_scalar_record):
    """prepend() 保持传入顺序，并跳过已存在 num。"""
    buf = RecordBuffer()
    existing = make_scalar_record(step=1)
    existing.num = 1
    buf.extend([existing])

    first = make_scalar_record(step=2)
    second = make_scalar_record(step=3)
    duplicate = make_scalar_record(step=4)
    replacement = make_scalar_record(step=5)
    first.num = 2
    second.num = 3
    duplicate.num = 1
    replacement.num = 2

    accepted = buf.prepend([first, second, duplicate, replacement])

    assert accepted == 2
    assert buf.drain() == [replacement, second, existing]


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


@pytest.mark.parametrize("separate_calls", [False, True])
@pytest.mark.parametrize("record_kind", ["scalar", "config"])
def test_record_buffer_keeps_latest_record_in_place(
    make_config_record, make_scalar_record, separate_calls, record_kind
):
    buf = RecordBuffer()
    make_record = make_config_record if record_kind == "config" else make_scalar_record
    first = make_record()
    latest = make_record()
    first.num = latest.num = 7
    first.timestamp.seconds = 1
    latest.timestamp.seconds = 2
    scalar = make_scalar_record()
    scalar.num = 1

    if separate_calls:
        assert buf.extend([first, scalar]) == 2
        assert buf.extend([latest]) == 0
    else:
        assert buf.extend([first, scalar, latest]) == 2

    assert len(buf) == 2
    assert buf.drain() == [latest, scalar]
    assert first.timestamp.seconds == 1
    assert buf.extend([latest]) == 1


def test_record_buffer_prepend_preserves_buffered_config(make_config_record, make_scalar_record):
    buf = RecordBuffer()
    old = make_config_record()
    latest = make_config_record()
    old.num = latest.num = -3
    old.save.payload = b"old"
    latest.save.payload = b"latest"
    scalar = make_scalar_record()
    scalar.num = 1
    buf.extend([latest])

    assert buf.prepend([scalar, old, scalar]) == 1
    assert buf.drain() == [scalar, latest]

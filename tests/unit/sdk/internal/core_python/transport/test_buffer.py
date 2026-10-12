import pytest

from swanlab.sdk.internal.core_python.transport.worker.buffer import RecordBuffer



def test_record_buffer_drain_clears_contents(make_scalar_record):
    """drain() 取出记录并清空缓冲区。"""
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

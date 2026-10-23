import os
from unittest.mock import Mock

import pytest

from swanlab.sdk.internal.bus import EmitterProtocol
from swanlab.sdk.internal.run.components.terminal import TerminalProxy
from swanlab.sdk.internal.run.components.terminal.emulator import Char, TerminalEmulator


@pytest.fixture
def proxy():
    return TerminalProxy(Mock(spec=EmitterProtocol), "all", max_log_length=1024, init_pid=os.getpid())


def emitted_lines(proxy):
    return [call.args[0].line for call in proxy._emitter.emit.call_args_list]


@pytest.mark.parametrize("stream", ["stdout", "stderr"])
def test_large_write_flushes_before_allocating_rows(proxy, monkeypatch, stream):
    monkeypatch.setattr(TerminalEmulator, "_MAX_LINES", 8)
    emulator = getattr(proxy, f"_{stream}_emulator")
    original_copy = Char.copy

    def checked_copy(char, **kwargs):
        # Check during allocation, not after the final read() has trimmed rows.
        assert emulator.num_lines <= emulator._MAX_LINES
        assert len(emulator._buffer) <= emulator._MAX_LINES
        return original_copy(char, **kwargs)

    monkeypatch.setattr(Char, "copy", checked_copy)
    lines = [f"line {i}" for i in range(100)]
    getattr(proxy, f"_{stream}_queue").put("\n".join(lines) + "\n")
    proxy._stopped.set()
    proxy._worker_loop()

    assert emitted_lines(proxy) == lines
    assert emulator.num_lines <= emulator._MAX_LINES
    proxy._flush()
    assert emitted_lines(proxy) == lines


def test_flush_keeps_pending_progress_until_finalize(proxy, monkeypatch):
    monkeypatch.setattr(TerminalEmulator, "_MAX_LINES", 4)
    lines = [f"line {i}" for i in range(12)]
    proxy._stdout_emulator.write("\n".join(lines) + "\npending\r\033[2Kdone")
    proxy._flush()
    assert emitted_lines(proxy) == lines

    proxy._installed = True
    proxy.uninstall()
    assert emitted_lines(proxy) == lines + ["done"]


def test_repeated_writes_preserve_ansi_and_empty_lines(proxy, monkeypatch):
    monkeypatch.setattr(TerminalEmulator, "_MAX_LINES", 4)
    lines = []
    for i in range(20):
        proxy._stdout_emulator.write(f"\033[31mitem {i}\033[0m\r\033[2K你好 {i}\n\n")
        lines.extend([f"你好 {i}", ""])
    proxy._flush()
    assert emitted_lines(proxy) == lines


def test_streams_keep_independent_pending_lines(proxy, monkeypatch):
    monkeypatch.setattr(TerminalEmulator, "_MAX_LINES", 4)
    proxy._stderr_emulator.write("error pending")
    lines = [f"out {i}" for i in range(12)]
    proxy._stdout_emulator.write("\n".join(lines) + "\n")
    proxy._flush()
    assert emitted_lines(proxy) == lines

    proxy._stderr_emulator.write("\r\033[2Kerror complete\n")
    proxy._flush()
    assert emitted_lines(proxy) == lines + ["error complete"]

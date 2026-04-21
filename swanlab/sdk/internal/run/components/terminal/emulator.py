"""
@author: cunyue
@file: emulator.py
@time: 2026/4/20
@description: 终端模拟器 — 2D 字符缓冲区 + ANSI CSI 解析

基于 wandb TerminalEmulator 的设计，做了以下调整：
- 移除 numpy 依赖，纯 Python 实现，因为SwanLab不需要重建ANSI序列
- read() 返回结构化 [(纯文本行, is_new_line), ...] 而非带 \r 的原始字符串
- _get_plain_line() 剥离 ANSI 码，ConsoleEvent.line 为纯文本
"""

from __future__ import annotations

import itertools
import re
from collections import defaultdict
from functools import cached_property
from typing import Callable, List, Tuple

from swanlab.sdk.internal.pkg import safe

# ==========================================
# ANSI 常量
# ==========================================

ANSI_CSI_RE = re.compile(r"\001?\033\[((?:\d|;)*)([a-zA-Z])\002?")
ANSI_OSC_RE = re.compile(r"\001?\033]([^\a]*)(\a)\002?")

SEP_RE = re.compile(r"\r|\n|" + "|".join([chr(i) for i in range(2**8) if repr(chr(i)).startswith("'\\x")]))

ANSI_FG = list(map(str, itertools.chain(range(30, 40), range(90, 98))))
ANSI_BG = list(map(str, itertools.chain(range(40, 50), range(100, 108))))

ANSI_FG_DEFAULT = "39"
ANSI_BG_DEFAULT = "49"
ANSI_RESET = "0"

ANSI_STYLES = {
    "1": "bold",
    "2": "/bold",
    "3": "italics",
    "4": "underscore",
    "5": "blink",
    "7": "reverse",
    "9": "strikethrough",
    "22": "/bold",
    "23": "/italics",
    "24": "/underscore",
    "25": "/blink",
    "27": "/reverse",
    "29": "/strikethrough",
}

# ==========================================
# Char / Cursor
# ==========================================


class Char:
    """单字符单元，包含前景/背景色和文本样式属性。"""

    __slots__ = (
        "data",
        "fg",
        "bg",
        "bold",
        "italics",
        "underscore",
        "blink",
        "strikethrough",
        "reverse",
    )

    def __init__(
        self,
        data: str = " ",
        fg: str = ANSI_FG_DEFAULT,
        bg: str = ANSI_BG_DEFAULT,
        bold: bool = False,
        italics: bool = False,
        underscore: bool = False,
        blink: bool = False,
        strikethrough: bool = False,
        reverse: bool = False,
    ) -> None:
        self.data = data
        self.fg = fg
        self.bg = bg
        self.bold = bold
        self.italics = italics
        self.underscore = underscore
        self.blink = blink
        self.strikethrough = strikethrough
        self.reverse = reverse

    def reset(self) -> None:
        """重置样式属性为默认值，保留 data。"""
        default = self.__class__()
        for k in self.__slots__[1:]:
            setattr(self, k, getattr(default, k))

    def copy(self, **kwargs: object) -> Char:
        """浅拷贝，可选覆盖部分字段。"""
        attrs = {k: kwargs.get(k, getattr(self, k)) for k in self.__slots__}
        return self.__class__(**attrs)  # type: ignore

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Char):
            return NotImplemented
        return all(getattr(self, k) == getattr(other, k) for k in self.__slots__)


_default_char = Char()


class Cursor:
    """2D 光标，跟踪位置和继承的字符样式。"""

    __slots__ = ("x", "y", "char")

    def __init__(self, x: int = 0, y: int = 0, char: Char | None = None) -> None:
        if char is None:
            char = Char()
        self.x = x
        self.y = y
        self.char = char


# ==========================================
# TerminalEmulator
# ==========================================


class TerminalEmulator:
    """终端模拟器 — 正常终端是将输出写入到tty，这里是将输出写入到事件发射器
    字符存储在 2D 矩阵中，由光标索引

    支持的 ANSI CSI 序列：
    - 光标移动：A (上), B (下), C (右), D (左), H/f (定位)
    - 行擦除：K
    - 屏幕擦除：J
    - 插入行：L
    - 颜色/样式：m (SGR)

    控制字符：
    - \\r → 归位 (carriage return)
    - \\n → 换行 (linefeed)
    - \\b → 退格 (cursor left)
    """

    _MAX_LINES = 1024

    def __init__(self) -> None:
        self._buffer: defaultdict[int, defaultdict[int, Char]] = defaultdict(lambda: defaultdict(lambda: _default_char))
        self._cursor = Cursor()
        self._num_lines_cache: int | None = None

        # Diff 追踪
        self._prev_num_lines: int = 0
        self._prev_last_line: str = ""
        # 已提交的行数（光标已通过 \n 离开的行）
        # 光标当前所在行是"进行中"的，不算 committed
        self._committed_lines: int = 0

    def finalize(self) -> None:
        """将所有行标记为 committed，使 read() 返回 pending 行。
        尾部空行不提交，避免发射无意义的空行事件。
        """
        while self.num_lines > self._committed_lines:
            last = self._get_plain_line(self.num_lines - 1)
            if last:
                self._committed_lines = self.num_lines
                break
            # 尾部空行，回退
            if self.num_lines - 1 in self._buffer:
                del self._buffer[self.num_lines - 1]
            self._cursor.y = min(self._cursor.y, self.num_lines - 2)
            self._num_lines_cache = None

    # ----------------------------------
    # 光标操作
    # ----------------------------------

    def cursor_up(self, n: int = 1) -> None:
        n = min(n, self._cursor.y)
        self._cursor.y -= n

    def cursor_down(self, n: int = 1) -> None:
        self._cursor.y += n

    def cursor_left(self, n: int = 1) -> None:
        n = min(n, self._cursor.x)
        self._cursor.x -= n

    def cursor_right(self, n: int = 1) -> None:
        self._cursor.x += n

    def carriage_return(self) -> None:
        self._cursor.x = 0

    def cursor_position(self, line: int, column: int) -> None:
        self._cursor.x = max(column, 1) - 1
        self._cursor.y = max(line, 1) - 1

    def linefeed(self) -> None:
        # 当前行已通过 \n 完成，更新 committed 边界
        self._committed_lines = max(self._committed_lines, self._cursor.y + 1)
        self.cursor_down()
        self.carriage_return()

    @cached_property
    def _cursor_fns(self) -> dict[str, Callable[[int], None]]:
        return {
            "A": self.cursor_up,
            "B": self.cursor_down,
            "C": self.cursor_right,
            "D": self.cursor_left,
        }

    # ----------------------------------
    # 擦除
    # ----------------------------------

    def erase_line(self, mode: int = 0) -> None:
        curr_line = self._buffer[self._cursor.y]
        if mode == 0:
            for i in range(self._cursor.x, self._get_line_len(self._cursor.y)):
                if i in curr_line:
                    del curr_line[i]
        elif mode == 1:
            for i in range(self._cursor.x + 1):
                if i in curr_line:
                    del curr_line[i]
        else:
            curr_line.clear()

    def erase_screen(self, mode: int = 0) -> None:
        if mode == 0:
            for i in range(self._cursor.y + 1, self.num_lines):
                if i in self._buffer:
                    del self._buffer[i]
            self.erase_line(mode)
        elif mode == 1:
            for i in range(self._cursor.y):
                if i in self._buffer:
                    del self._buffer[i]
            self.erase_line(mode)
        elif mode in (2, 3):
            self._buffer.clear()

    def insert_lines(self, n: int = 1) -> None:
        for i in range(self.num_lines - 1, self._cursor.y, -1):
            self._buffer[i + n] = self._buffer[i]
        for i in range(self._cursor.y + 1, self._cursor.y + 1 + n):
            if i in self._buffer:
                del self._buffer[i]

    # ----------------------------------
    # 写入
    # ----------------------------------

    def write(self, data: str) -> None:
        """处理一块终端输出（纯文本 + ANSI 序列）。"""
        self._num_lines_cache = None
        data = self._remove_osc(data)
        prev_end = 0
        for match in ANSI_CSI_RE.finditer(data):
            start, end = match.span()
            text = data[prev_end:start]
            prev_end = end
            self._write_text(text)
            self._handle_csi(*match.groups())
        self._write_text(data[prev_end:])

    def _handle_csi(self, params: str, command: str) -> None:
        """分发 CSI 命令。单个 CSI 解析异常不影响后续序列处理。"""
        with safe.block(message="Terminal emulator CSI parse error", write_to_tty=False):
            if command == "m":
                # 只取首个参数；SGR 多参数（如 \033[31;1m）暂不逐个遍历，
                # 因为 read() 只消费 Char.data 纯文本，样式属性目前无消费者
                p = params.split(";")[0]
                if not p:
                    p = "0"
                if p in ANSI_FG:
                    self._cursor.char.fg = p
                elif p in ANSI_BG:
                    self._cursor.char.bg = p
                elif p == ANSI_RESET:
                    self._cursor.char.reset()
                elif p in ANSI_STYLES:
                    style = ANSI_STYLES[p]
                    off = style.startswith("/")
                    if off:
                        style = style[1:]
                    setattr(self._cursor.char, style, not off)
            else:
                cursor_fn = self._cursor_fns.get(command)
                if cursor_fn:
                    cursor_fn(int(params.split(";")[0]) if params else 1)
                elif command == "J":
                    p = int(params.split(";")[0]) if params else 0
                    self.erase_screen(p)
                elif command == "K":
                    p = int(params.split(";")[0]) if params else 0
                    self.erase_line(p)
                elif command == "L":
                    p = int(params.split(";")[0]) if params else 1
                    self.insert_lines(p)
                elif command in "Hf":
                    p = params.split(";")
                    if len(p) == 2:
                        self.cursor_position(int(p[0]), int(p[1]))
                    elif len(p) == 1 and p[0]:
                        self.cursor_position(int(p[0]), 1)
                    else:
                        self.cursor_position(1, 1)

    def _write_text(self, text: str) -> None:
        """处理包含 \\r/\\n/\\b 分隔符的文本。"""
        prev_end = 0
        for match in SEP_RE.finditer(text):
            start, end = match.span()
            self._write_plain_text(text[prev_end:start])
            prev_end = end
            c = match.group()
            if c == "\n":
                self.linefeed()
            elif c == "\r":
                self.carriage_return()
            elif c == "\b":
                self.cursor_left()
            else:
                continue
        self._write_plain_text(text[prev_end:])

    def _write_plain_text(self, plain_text: str) -> None:
        """在当前光标位置写入纯文本。"""
        self._buffer[self._cursor.y].update(
            [(self._cursor.x + i, self._cursor.char.copy(data=c)) for i, c in enumerate(plain_text)]
        )
        self._cursor.x += len(plain_text)

    @staticmethod
    def _remove_osc(text: str) -> str:
        """移除 OSC (Operating System Command) 序列。"""
        return re.sub(ANSI_OSC_RE, "", text)

    # ----------------------------------
    # 读取
    # ----------------------------------

    @property
    def num_lines(self) -> int:
        """
        当前有多少行（含空行）
        """
        if self._num_lines_cache is not None:
            return self._num_lines_cache
        # 考虑已提交的行、光标所在的行以及缓冲区中已有的行
        max_idx = max(self._buffer.keys()) if self._buffer else -1
        ret = max(max_idx + 1, self._committed_lines, self._cursor.y + 1)
        self._num_lines_cache = ret
        return ret

    def _get_line_len(self, n: int) -> int:
        if n not in self._buffer:
            return 0
        line = self._buffer[n]
        if not line:
            return 0
        max_idx = max(line.keys())
        for i in range(max_idx, -1, -1):
            if line[i] != _default_char:
                return i + 1
        return 0

    def _get_plain_line(self, n: int) -> str:
        """返回第 n 行的纯文本（ANSI 码已剥离，尾部空白已修剪）。"""
        if n not in self._buffer:
            return ""
        line = self._buffer[n]
        line_len = self._get_line_len(n)
        return "".join(line[i].data for i in range(line_len))

    def read(self) -> List[Tuple[str, bool]]:
        """返回自上次 read() 以来的增量 diff。

        返回 [(line_text, is_new_line), ...] 列表：
        - is_new_line=True: 已完成的行（光标已通过 \n 离开）
        - is_new_line=False: 进行中的行（光标还在该行，可能被 \r 覆盖）

        读取后，若总行数超过 _MAX_LINES，自动滚动缓冲区。
        """
        num_lines = self.num_lines
        result: List[Tuple[str, bool]] = []

        # 计算新增的已完成行
        for i in range(self._prev_num_lines, min(self._committed_lines, num_lines)):
            result.append((self._get_plain_line(i), True))

        # 进行中的行（光标当前所在行，可能被 \r 覆盖）
        # 只有当进行中的行内容发生了变化时才返回
        if num_lines > self._committed_lines:
            pending_line = self._get_plain_line(num_lines - 1)
            if pending_line != self._prev_last_line:
                result.append((pending_line, False))

        # 滚动缓冲区
        if num_lines > self._MAX_LINES:
            shift = num_lines - self._MAX_LINES
            for i in range(shift, num_lines):
                self._buffer[i - shift] = self._buffer[i]
            for i in range(self._MAX_LINES, max(self._buffer.keys()) + 1):
                if i in self._buffer:
                    del self._buffer[i]
            self._cursor.y -= min(self._cursor.y, shift)
            self._committed_lines -= min(self._committed_lines, shift)
            self._num_lines_cache = self._MAX_LINES
            num_lines = self._MAX_LINES

        self._prev_num_lines = min(self._committed_lines, num_lines)
        self._prev_last_line = self._get_plain_line(num_lines - 1) if num_lines > 0 else ""
        return result

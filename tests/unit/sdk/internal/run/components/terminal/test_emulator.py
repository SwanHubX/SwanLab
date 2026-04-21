"""
@author: cunyue
@file: test_emulator.py
@time: 2026/4/21
@description: 终端模拟器单元测试

测试策略：
    - 优先通过公共接口 write() + read() 验证行为
    - 光标移动等 read() 无法捕获的场景，通过 _get_plain_line / _cursor 辅助验证
    - 不直接依赖 _buffer / _committed_lines 等内部数据结构
"""

from __future__ import annotations

from swanlab.sdk.internal.run.components.terminal.emulator import TerminalEmulator

# ==========================================
# 辅助
# ==========================================


def _feed(em: TerminalEmulator, data: str) -> list[tuple[str, bool]]:
    """写入数据后立即读取 diff。"""
    em.write(data)
    return em.read()


def _line(em: TerminalEmulator, n: int) -> str:
    """获取第 n 行纯文本（辅助验证，非公共 API）。"""
    return em._get_plain_line(n)  # noqa


# ==========================================
# 基础写入 & read() 语义
# ==========================================


class TestWriteRead:
    def test_single_line(self):
        assert _feed(TerminalEmulator(), "Hello\n") == [("Hello", True)]

    def test_multiple_lines(self):
        assert _feed(TerminalEmulator(), "A\nB\nC\n") == [
            ("A", True),
            ("B", True),
            ("C", True),
        ]

    def test_no_trailing_newline_is_pending(self):
        """没有 \n 的行是进行中（is_new_line=False）。"""
        assert _feed(TerminalEmulator(), "Hello") == [("Hello", False)]

    def test_empty_lines_preserved(self):
        """空行保留，不跳过。"""
        assert _feed(TerminalEmulator(), "A\n\nB\n") == [("A", True), ("", True), ("B", True)]

    def test_empty_write(self):
        assert _feed(TerminalEmulator(), "") == []

    def test_incremental_read(self):
        """多次写入 + 读取，diff 互不影响。"""
        em = TerminalEmulator()
        assert _feed(em, "A\n") == [("A", True)]
        assert _feed(em, "B\n") == [("B", True)]

    def test_read_without_write_is_empty(self):
        em = TerminalEmulator()
        _feed(em, "A\n")
        assert em.read() == []

    def test_crlf(self):
        """Windows 风格 \r\n：\r 回行首，\n 换行。"""
        assert _feed(TerminalEmulator(), "A\r\nB\r\n") == [("A", True), ("B", True)]

    def test_pending_then_committed(self):
        """进行中行 → 补 \n → 变为已完成。"""
        em = TerminalEmulator()
        assert _feed(em, "A") == [("A", False)]
        # 补 \n 后 A 从 pending 变为 committed
        assert _feed(em, "\n") == [("A", True), ("", False)]

    def test_mixed_committed_and_pending(self):
        assert _feed(TerminalEmulator(), "A\nB\nC") == [
            ("A", True),
            ("B", True),
            ("C", False),
        ]

    def test_long_line(self):
        """超长单行（无 \n）仍然正确。"""
        text = "X" * 10000
        assert _feed(TerminalEmulator(), text) == [(text, False)]

    def test_unicode(self):
        """中文 / emoji 不影响行分割。"""
        assert _feed(TerminalEmulator(), "你好\n🌍\n") == [
            ("你好", True),
            ("🌍", True),
        ]


# ==========================================
# \r — 回车
# ==========================================


class TestCarriageReturn:
    def test_overwrite_prefix(self):
        """\r 回行首后覆盖前缀，尾部保留。"""
        assert _feed(TerminalEmulator(), "abcdef\rxy\n") == [("xycdef", True)]

    def test_full_overwrite(self):
        """\r 回行首后覆盖整行。"""
        assert _feed(TerminalEmulator(), "abc\rxyz\n") == [("xyz", True)]

    def test_r_then_newline(self):
        assert _feed(TerminalEmulator(), "foo\rbar\n") == [("bar", True)]

    def test_consecutive_cr(self):
        """连续多次 \r 不带内容，再写入。"""
        assert _feed(TerminalEmulator(), "foo\r\r\rbar\n") == [("bar", True)]

    def test_cr_clear_then_newline(self):
        """\r + 擦除 + \n → 空行 committed。"""
        em = TerminalEmulator()
        em.write("some text\r\033[K\n")
        assert em.read() == [("", True)]

    def test_cr_without_overwrite_keeps_tail(self):
        """\r 回行首，写入更短内容，原尾部保留。"""
        assert _feed(TerminalEmulator(), "abcdef\rxy\n") == [("xycdef", True)]


# ==========================================
# \b — 退格
# ==========================================


class TestBackspace:
    def test_overwrite_previous_char(self):
        assert _feed(TerminalEmulator(), "abc\bX\n") == [("abX", True)]

    def test_at_line_start_noop(self):
        """行首退格不移动光标。"""
        assert _feed(TerminalEmulator(), "\bA\n") == [("A", True)]

    def test_multiple_backspace(self):
        assert _feed(TerminalEmulator(), "abcde\b\bXY\n") == [("abcXY", True)]

    def test_excessive_backspace_clamped(self):
        """退格超过行首被钳制，多余退格无效。"""
        em = TerminalEmulator()
        em.write("AB\b\b\b\bX\n")
        # AB 后光标在 x=2，退格到 x=0，后续退格无效 → X 覆盖 A → XB
        assert em.read() == [("XB", True)]


# ==========================================
# CSI 光标移动
# ==========================================


class TestCursorMovement:
    def test_cursor_up(self):
        """CSI A — 光标上移后写入。"""
        em = TerminalEmulator()
        em.write("Line1\nLine2\n")
        em.read()
        # 光标在 y=2,x=0（linefeed 后），上移1行到 y=1,x=0
        em.write("\033[A**")
        # 在 Line2 行首覆盖写入 ** → **ne2
        assert _line(em, 1) == "**ne2"

    def test_cursor_up_n(self):
        """CSI nA — 上移 n 行。"""
        em = TerminalEmulator()
        em.write("A\nB\nC\n")
        em.read()
        # 光标在 y=3,x=0，上移2行到 y=1,x=0
        em.write("\033[2A*")
        assert em._cursor.y == 1
        # x=0 处写入 *，覆盖 B → "*"
        assert _line(em, 1) == "*"

    def test_cursor_up_clamped_at_top(self):
        """上移不超过第 0 行。"""
        em = TerminalEmulator()
        em.write("A\nB\n")
        em.read()
        em.write("\033[99A*")
        assert em._cursor.y == 0

    def test_cursor_down(self):
        """CSI B — 光标下移，x 保持不变。"""
        em = TerminalEmulator()
        em.write("AB")
        em.read()
        # cursor at (0, 2)，下移后写入
        em.write("\033[BX\n")
        assert _line(em, 1) == "  X"  # x=2 位置写入 X

    def test_cursor_right(self):
        """CSI C — 光标右移跳过字符。"""
        assert _feed(TerminalEmulator(), "AB\033[2CX\n") == [("AB  X", True)]

    def test_cursor_left(self):
        """CSI D — 光标左移后覆盖。"""
        assert _feed(TerminalEmulator(), "ABCDE\033[2DXY\n") == [("ABCXY", True)]

    def test_cursor_left_default_one(self):
        """CSI D 无参数默认左移 1 格。"""
        assert _feed(TerminalEmulator(), "AB\033[DX\n") == [("AX", True)]

    def test_cursor_left_clamped(self):
        """左移不超过行首。"""
        assert _feed(TerminalEmulator(), "\033[99DX\n") == [("X", True)]


# ==========================================
# CSI H/f — 光标定位
# ==========================================


class TestCursorPosition:
    def test_row_and_column(self):
        """CSI row;colH 定位。"""
        em = TerminalEmulator()
        em.write("\033[2;5H*")
        assert em._cursor.y == 1
        assert em._cursor.x == 5

    def test_row_only(self):
        """CSI rowH — 列默认 1。"""
        em = TerminalEmulator()
        em.write("\033[3H*")
        assert em._cursor.y == 2
        assert em._cursor.x == 1

    def test_home(self):
        """CSI H — 无参数回到 (1,1)。"""
        em = TerminalEmulator()
        em.write("Some text\033[H*")
        assert em._cursor.y == 0
        assert em._cursor.x == 1

    def test_f_equals_h(self):
        """CSI f 等同于 CSI H。"""
        em = TerminalEmulator()
        em.write("\033[2;3f*")
        assert em._cursor.y == 1
        assert em._cursor.x == 3

    def test_zero_params_clamped_to_one(self):
        """参数 0 被钳制为 1。"""
        em = TerminalEmulator()
        em.write("\033[0;0H*")
        assert em._cursor.y == 0
        assert em._cursor.x == 1


# ==========================================
# CSI K — 行擦除
# ==========================================


class TestEraseLine:
    def test_erase_to_end_default(self):
        """CSI K 默认 mode=0：擦除光标到行尾。"""
        em = TerminalEmulator()
        em.write("ABCDE\033[2D\033[K\n")
        # 光标左移2格到位置3(C处)，\033[K 擦除位置3~行尾 → 保留 ABC
        assert em.read() == [("ABC", True)]

    def test_erase_to_end_explicit(self):
        """CSI 0K 显式 mode=0。"""
        em = TerminalEmulator()
        em.write("ABCDE\033[2D\033[0K\n")
        assert em.read() == [("ABC", True)]

    def test_erase_from_start(self):
        """CSI 1K — 擦除行首到光标（含光标位）。"""
        em = TerminalEmulator()
        em.write("ABCDE\033[2D\033[1K\n")
        # 光标在位置 3，擦除 0~3(含) → 保留位置 4 的 E，前导空格
        assert em.read() == [("    E", True)]

    def test_erase_full_line(self):
        """CSI 2K — 擦除整行。"""
        em = TerminalEmulator()
        em.write("ABCDE\033[2K\n")
        assert em.read() == [("", True)]


# ==========================================
# CSI J — 屏幕擦除
# ==========================================


class TestEraseScreen:
    def test_erase_below(self):
        """CSI 0J — 擦除光标下方。"""
        em = TerminalEmulator()
        em.write("A\nB\nC\n")
        em.read()
        em.write("\033[2A\033[J")
        # 光标移到 y=1，擦除 y=1 之后 → C 被清除
        assert _line(em, 0) == "A"
        assert _line(em, 2) == ""

    def test_erase_above(self):
        """CSI 1J — 擦除光标上方。"""
        em = TerminalEmulator()
        em.write("A\nB\nC\n")
        em.read()
        em.write("\033[1J")
        # 光标在 y=3，擦除 y=0~2
        assert _line(em, 0) == ""

    def test_erase_all(self):
        """CSI 2J — 擦除全部。"""
        em = TerminalEmulator()
        em.write("A\nB\nC\n")
        assert em.num_lines == 4  # A, B, C + 空行
        em.read()
        em.write("\033[2J")
        # buffer 已清，但是历史行数还是4行
        assert em.num_lines == 4
        em.write("D")
        assert em.num_lines == 4
        assert em.read() == [("D", False)]


# ==========================================
# CSI L — 插入行
# ==========================================


class TestInsertLines:
    def test_insert_pushes_content_down(self):
        """CSI L — 插入空行，下方内容下推。"""
        em = TerminalEmulator()
        em.write("A\nB\nC\n")
        em.read()
        # 光标在 y=3，上移到 y=1
        em.write("\033[2A")
        assert em._cursor.y == 1
        em.write("\033[L")
        # y=1 内容保留(B)，y=2 被清空，原 C 推到 y=3
        assert _line(em, 1) == "B"
        assert _line(em, 2) == ""
        assert _line(em, 3) == "C"

    def test_insert_multiple(self):
        """CSI nL — 插入多行。"""
        em = TerminalEmulator()
        em.write("A\nB\nC\n")
        em.read()
        em.write("\033[3A")  # 上移到 y=0
        em.write("\033[2L")
        # y=0 内容保留(A)，y=1~2 清空，原 B/C 推到 y=3~4
        assert _line(em, 0) == "A"
        assert _line(em, 1) == ""
        assert _line(em, 2) == ""
        assert _line(em, 3) == "B"
        assert _line(em, 4) == "C"


# ==========================================
# CSI m — 颜色/样式（SGR）
# ==========================================


class TestSGR:
    def test_foreground_color_ignored_in_plain_text(self):
        """read() 返回纯文本，颜色信息不影响。"""
        assert _feed(TerminalEmulator(), "\033[31mRed\033[0m\n") == [("Red", True)]

    def test_bold_ignored_in_plain_text(self):
        assert _feed(TerminalEmulator(), "\033[1mBold\033[0m\n") == [("Bold", True)]

    def test_reset(self):
        assert _feed(TerminalEmulator(), "\033[1;31mX\033[0mY\n") == [("XY", True)]

    def test_no_params_defaults_to_reset(self):
        assert _feed(TerminalEmulator(), "\033[1mX\033[mY\n") == [("XY", True)]

    def test_sgr_does_not_break_line_splitting(self):
        """样式切换不会打断行内容。"""
        assert _feed(TerminalEmulator(), "\033[31mA\033[32mB\n") == [("AB", True)]


# ==========================================
# OSC 移除
# ==========================================


class TestOSCRemoval:
    def test_osc_stripped(self):
        assert _feed(TerminalEmulator(), "\033]0;title\007Hello\n") == [("Hello", True)]

    def test_osc_in_middle(self):
        assert _feed(TerminalEmulator(), "A\033]0;x\007B\n") == [("AB", True)]

    def test_multiple_osc(self):
        assert _feed(TerminalEmulator(), "\033]0;a\007X\033]1;b\007Y\n") == [("XY", True)]


# ==========================================
# 进度条场景（tqdm 风格）— 核心业务场景
# ==========================================


class TestProgressBar:
    def test_r_overwrite_progress_bar(self):
        """\r 反复覆盖，只有最终 \n 的行被 committed。"""
        em = TerminalEmulator()
        em.write("Epoch 1/10\r")
        assert em.read() == [("Epoch 1/10", False)]

        em.write("Epoch 2/10\r")
        assert em.read() == [("Epoch 2/10", False)]

        em.write("Epoch 3/10\n")
        assert em.read() == [("Epoch 3/10", True), ("", False)]

    def test_tqdm_style_full_cycle(self):
        """完整 tqdm 周期：多次 \r 更新 → 最终 \n。"""
        em = TerminalEmulator()
        em.write("  0%|          | 0/100 it/s\r")
        r1 = em.read()
        assert all(not is_new for _, is_new in r1)

        em.write("100%|##########| 100/100 done\n")
        r2 = em.read()
        assert all(is_new for _, is_new in r2[:-1])  # 最后一行永远是空行
        assert r2[0] == ("100%|##########| 100/100 done", True)

    def test_progress_bar_then_print(self):
        """进度条完成后紧接 print 输出。"""
        em = TerminalEmulator()
        em.write("Training started\n")
        assert em.read() == [("Training started", True)]

        em.write("  0%|  | 0/10\r")
        em.read()

        em.write("100%|##| 10/10\n")
        result = em.read()
        assert result[0] == ("100%|##| 10/10", True)

    def test_identical_pending_not_reduplicated(self):
        """连续写入相同的进行中行，不重复返回。"""
        em = TerminalEmulator()
        em.write("same\r")
        assert em.read() == [("same", False)]

        em.write("same\r")
        assert em.read() == []

    def test_pending_cleared_to_empty(self):
        """进行中行被 \r + 擦除清空后换行 → 空行 committed。"""
        em = TerminalEmulator()
        em.write("some text\r")
        em.read()

        em.write("\033[K\n")
        # 空行 committed + 新空行 pending
        assert em.read() == [("", True), ("", False)]

    def test_flush_skips_pending_lines(self):
        """_flush 场景：只有 is_new_line=True 的行被发射。"""
        em = TerminalEmulator()
        em.write("pending line\r")
        diff = em.read()
        assert [t for t, is_new in diff if is_new] == []
        assert [t for t, is_new in diff if not is_new] == ["pending line"]

        em.write("committed line\n")
        diff = em.read()
        assert [t for t, is_new in diff if is_new] == ["committed line"]


# ==========================================
# 缓冲区滚动
# ==========================================


class TestScrollBuffer:
    def test_overflow_truncates_oldest(self):
        em = TerminalEmulator()
        for i in range(em._MAX_LINES + 10):
            em.write(f"Line {i}\n")
        em.read()
        assert em.num_lines <= em._MAX_LINES

    def test_preserves_recent_lines(self):
        em = TerminalEmulator()
        total = em._MAX_LINES + 5
        for i in range(total):
            em.write(f"Line {i}\n")
        em.read()
        # 滚动后最后一行有内容（可能不是最后一行，因为空行也在）
        assert _line(em, em.num_lines - 2) == f"Line {total - 1}"

    def test_cursor_adjusted_after_scroll(self):
        em = TerminalEmulator()
        for i in range(em._MAX_LINES + 10):
            em.write(f"Line {i}\n")
        em.read()
        assert em._cursor.y >= 0


# ==========================================
# 边界 & 异常
# ==========================================


class TestEdgeCases:
    def test_malformed_csi_no_crash(self):
        """损坏的 CSI 序列不应崩溃，最终正常文本仍能产出。"""
        em = TerminalEmulator()
        em.write("\033[99999m")  # 超大参数 → try/except 静默跳过
        em.write("ok\n")
        result = em.read()
        # ok 被正确产出
        assert ("ok", True) in result

    def test_null_byte(self):
        """空字节不崩溃。"""
        em = TerminalEmulator()
        em.write("A\x00B\n")
        # \x00 被 SEP_RE 匹配为控制字符，被 skip
        result = em.read()
        assert len(result) > 0

    def test_only_control_chars(self):
        """纯控制字符输入。"""
        em = TerminalEmulator()
        em.write("\r\n\r\n")
        # 两次 \r\n 产生两个空行 committed + 一个空行 pending
        result = em.read()
        assert all(line == "" for line, _ in result)

    def test_consecutive_newlines(self):
        em = TerminalEmulator()
        em.write("\n\n\n")
        # 3个 \n 产生3个空行 committed + 1个空行 pending
        result = em.read()
        assert all(line == "" for line, _ in result)
        assert sum(1 for _, is_new in result if is_new) == 3

    def test_large_write(self):
        """大批量写入不崩溃。"""
        em = TerminalEmulator()
        for i in range(200):
            em.write(f"Line {i}\n")
        result = em.read()
        assert len(result) > 0

    def test_unicode_with_ansi(self):
        """中文 + ANSI 颜色。"""
        assert _feed(TerminalEmulator(), "\033[31m你好\033[0m\n") == [("你好", True)]


# ==========================================
# 组合场景
# ==========================================


class TestCombinations:
    def test_color_then_overwrite(self):
        assert _feed(TerminalEmulator(), "\033[31mRed\033[0m\rGreen\n") == [("Green", True)]

    def test_mixed_control_chars(self):
        """\b + \r 组合。"""
        assert _feed(TerminalEmulator(), "ABC\bD\rEF\n") == [("EFD", True)]

    def test_write_after_clear(self):
        """擦除屏幕后继续写入。"""
        em = TerminalEmulator()
        em.write("A\nB\n")
        em.read()
        em.write("\033[2J")
        assert em.num_lines == 3  # A, B + 空行(cursor)
        em.write("New\n")
        # 清屏后写 New，New 在新行 committed
        result = em.read()
        assert ("New", True) in result


class TestFinalize:
    def test_pending_to_committed(self):
        """finalize 将 pending 行变为 committed。"""
        em = TerminalEmulator()
        em.write("pending line")
        assert em.read() == [("pending line", False)]
        em.finalize()
        assert em.read() == [("pending line", True)]

    def test_only_pending_emitted(self):
        """finalize 只发射有内容的 pending 行，忽略尾部空行。"""
        em = TerminalEmulator()
        em.write("line1\nline2\npending")
        assert em.read() == [("line1", True), ("line2", True), ("pending", False)]
        em.finalize()
        assert em.read() == [("pending", True)]
        # 尾部空行场景
        em2 = TerminalEmulator()
        em2.write("A\nB\n\nC")
        em2.read()
        em2.finalize()
        assert em2.read() == [("C", True)]

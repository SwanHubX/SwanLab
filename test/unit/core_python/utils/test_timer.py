"""
@author: cunyue
@file: test_timer.py
@time: 2025/12/30 23:36
@description: 测试定时相关功能
"""

import threading
import time
from unittest.mock import MagicMock

import pytest

# 假设你的 Timer 类定义在 swanlab.utils.timer 模块中
# 请根据实际情况修改 import 路径
from swanlab.core_python.utils.timer import Timer


# ==========================================
# Fixtures (测试夹具)
# ==========================================


@pytest.fixture
def mock_swanlog(mocker):
    """
    Mock 掉 swanlog，既为了防止真实写日志，
    也用于断言日志方法是否被正确调用。
    注意：请将 'swanlab.utils.timer.swanlog' 替换为你实际代码中引入 swanlog 的路径
    """
    return mocker.patch("swanlab.core_python.utils.timer.swanlog")


@pytest.fixture
def timer_manager():
    """
    Timer 管理器 fixture。
    作用：收集测试中创建的所有 Timer 实例，在测试结束（teardown）时统一 cancel 和 join。
    防止因为测试失败导致后台线程残留，引起后续测试卡死或报错。
    """
    timers = []

    def _create_timer(*args, **kwargs):
        t = Timer(*args, **kwargs)
        timers.append(t)
        return t

    yield _create_timer

    # Teardown: 清理所有创建的 timer
    for t in timers:
        try:
            t.cancel()
            t.join(timeout=1.0)  # 设置超时防止死锁
        except Exception:
            pass


# ==========================================
# Test Cases (测试用例)
# ==========================================


def test_basic_interval(timer_manager):
    """测试：固定间隔的基本执行功能"""
    counter = 0
    lock = threading.Lock()

    def task():
        nonlocal counter
        with lock:
            counter += 1

    # 创建并启动 Timer，间隔 0.1s
    timer = timer_manager(task, interval=0.1, immediate=False)
    timer.run()

    # 等待 0.25s
    # 理论执行时间点：T+0.1s, T+0.2s
    time.sleep(0.25)

    # 手动停止
    timer.cancel()
    timer.join()

    # 断言至少执行了 2 次 (多线程环境时间不绝对精确，用 >= 比较稳妥)
    assert counter >= 2


def test_immediate_execution(timer_manager):
    """测试：immediate=True 是否立即执行"""
    counter = 0

    def task():
        nonlocal counter
        counter += 1

    # 设置一个很长的间隔（10秒），确保在测试的短时间内，只有 immediate 那一次触发
    timer = timer_manager(task, interval=10.0, immediate=True)
    timer.run()

    # 稍微给一点时间让线程启动并执行第一次任务
    time.sleep(0.1)

    # 断言确实执行了一次
    assert counter == 1


def test_dynamic_interval(timer_manager):
    """测试：传入函数作为 interval，实现动态间隔"""
    intervals = []

    def interval_strategy(count):
        intervals.append(count)
        # count 是从 0 开始的
        # 第 0 次 (Wait) -> 返回 0.1s
        # 第 1 次 (Wait) -> 返回 0.1s
        # 第 2 次 (Wait) -> 返回 0.5s (变慢)
        if count < 2:
            return 0.1
        return 0.5

    # 任务本身不做事
    timer = timer_manager(lambda: None, interval=interval_strategy)
    timer.run()

    # 我们等待 0.35s：
    # T+0:   启动 -> wait(0.1) -> 此时 count=0
    # T+0.1: 执行 task, count变为1 -> wait(0.1) -> 此时 count=1
    # T+0.2: 执行 task, count变为2 -> wait(0.5) -> 此时 count=2
    time.sleep(0.35)

    timer.cancel()
    timer.join()

    # 验证 interval 策略函数是否被正确调用，并且传入了正确的 count
    assert 0 in intervals
    assert 1 in intervals
    assert 2 in intervals

    # 因为第三次 wait 是 0.5s，所以在 0.35s 内不应该完成第三次等待
    # Timer 的 _count 应该在 2 左右 (0和1执行完了，2正在wait)
    assert timer._count >= 2


def test_cancel_interrupts_sleep(timer_manager):
    """测试：cancel() 是否能立即打断长时间的 wait"""
    start_time = time.time()

    # 设置 10 秒间隔
    timer = timer_manager(lambda: None, interval=10.0)
    timer.run()

    # 确保线程已经进入 wait 状态
    time.sleep(0.1)

    # 调用 cancel，代码中 _stop_event.wait() 应该立即返回
    timer.cancel()
    timer.join()

    duration = time.time() - start_time

    # 如果 cancel 没有打断 wait，这里会耗时 10s 以上
    # 只要小于 1s 都说明打断成功了
    assert duration < 1.0


def test_error_resilience(timer_manager, mock_swanlog):
    """测试：任务抛出异常时，定时器不崩溃，且计数器继续增加"""

    # 模拟一个只会报错的任务
    mock_task = MagicMock(side_effect=ValueError("Test Error"))

    timer = timer_manager(mock_task, interval=0.1)
    timer.run()

    # 等待执行几次
    time.sleep(0.35)

    # 1. 验证任务被尝试执行了多次（说明线程没因为异常挂掉）
    assert mock_task.call_count >= 3

    # 2. 验证是否调用了 swanlog.error
    assert mock_swanlog.error.called
    assert "Error executing task" in mock_swanlog.error.call_args[0][0]

    # 3. 验证 _count 是否增加 (你的代码中 finally 块处理了计数)
    assert timer._count >= 3


def test_double_run_warning(timer_manager, mock_swanlog):
    """测试：重复调用 run() 是否发出警告"""
    timer = timer_manager(lambda: None, interval=1.0)
    timer.run()

    # 再次启动
    timer.run()

    # 验证 swanlog.warning 被调用
    mock_swanlog.warning.assert_called_with("Timer already running")

    # 验证活跃线程数没有异常增加 (应该只有1个Timer线程 + 主线程)
    # 这里不做严格断言，因为 pytest 环境可能有其他线程，但在逻辑上确保没开新线程


def test_restart_capability(timer_manager):
    """测试：停止后是否可以重启"""
    counter = 0

    def task():
        nonlocal counter
        counter += 1

    timer = timer_manager(task, interval=0.1)

    # 第一次运行
    timer.run()
    time.sleep(0.15)
    timer.cancel()
    timer.join()
    run_count_1 = counter
    assert run_count_1 >= 1

    # 重启
    timer.run()
    time.sleep(0.15)
    timer.cancel()
    timer.join()

    # 计数器应该继续增加
    assert counter > run_count_1

"""
swanlab public API type stubs.

This file is the single source of truth for all type signatures exposed at the
top-level `swanlab` namespace.  Add new public symbols here when they are added
to swanlab/__init__.py.
"""

from concurrent.futures import Future
from typing import Any, Callable, List, Mapping, Optional, Union

from . import utils
from .sdk import Audio, Callback, Image, Run, Settings, Text, Video, config, logout, verify
from .sdk.cmd.init import ConfigLike
from .sdk.typings.run import AsyncLogType, FinishType, ModeType, ResumeType
from .sdk.typings.run.column import ScalarXAxisType
from .sdk.typings.run.transforms import CaptionsType
from .sdk.typings.run.transforms.audio import AudioDatasType, AudioRatesType
from .sdk.typings.run.transforms.image import ImageDatasType, ImageFilesType, ImageModesType, ImageSizesType
from .sdk.typings.run.transforms.text import TextDatasType
from .sdk.typings.run.transforms.video import VideoDatasType

__version__: str

__all__ = [
    # cmd
    "merge_settings",
    "init",
    "finish",
    "login",
    "logout",
    "verify",
    "log",
    "log_text",
    "log_image",
    "log_audio",
    "log_video",
    "define_scalar",
    "async_log",
    # run
    "run",  # type: ignore [no-redef]
    "Run",
    "has_run",
    "get_run",
    # config
    "config",
    # data
    "Text",
    "Audio",
    "Image",
    "Video",
    # utils
    "utils",
    "Settings",
    "Callback",
]

# ── lifecycle ──────────────────────────────────────────────────────────────────

def init(
    *,
    reinit: Optional[bool] = None,
    logdir: Optional[str] = None,
    mode: Optional[ModeType] = None,
    workspace: Optional[str] = None,
    project: Optional[str] = None,
    public: Optional[bool] = None,
    name: Optional[str] = None,
    color: Optional[str] = None,
    description: Optional[str] = None,
    job_type: Optional[str] = None,
    group: Optional[str] = None,
    tags: Optional[List[str]] = None,
    id: Optional[str] = None,
    resume: Optional[Union[ResumeType, bool]] = None,
    config: Optional[ConfigLike] = None,
    settings: Optional[Settings] = None,
    callbacks: Optional[List[Callback]] = None,
    **kwargs: Any,
) -> Run:
    """Initialize a new SwanLab run to track experiments.

    This function starts a new run for logging metrics, artifacts, and metadata.
    After calling this, use `swanlab.log()` to log data and `swanlab.finish()` to
    close the run. SwanLab automatically finishes runs at program exit.

    :param reinit: If True, finish the current run before starting a new one. Defaults to False.
    :param logdir: Directory to store logs. Defaults to "./swanlog".
    :param mode: Run mode. Options: "cloud" (sync to cloud), "local" (local only),
        "offline" (save locally for later sync), "disabled" (no logging). Defaults to "cloud".
    :param workspace: Workspace or organization name. Defaults to current user.
    :param project: Project name. Defaults to current directory name.
    :param public: Make project publicly visible (cloud mode only). Defaults to False.
    :param name: Experiment name. Auto-generated if not provided.
    :param color: Experiment color for visualization. Auto-generated if not provided.
    :param description: Experiment description.
    :param job_type: Job type label (e.g., "train", "eval").
    :param group: Group name for organizing related experiments.
    :param tags: List of tags for categorizing experiments.
    :param id: Run ID for resuming a previous run (cloud mode only).
    :param resume: Resume behavior. Options: "must" (must resume), "allow" (resume if exists),
        "never" (always create new). Defaults to "never".
    :param config: Experiment configuration dict or path to config file (JSON/YAML).
    :param settings: Custom Settings object for advanced configuration.
    :param callbacks: List of callback functions triggered on run events.
    :return: The initialized Run object.
    :raises RuntimeError: If a run is already active and reinit=False.

    Examples:

        Basic local run:

        >>> import swanlab
        >>> swanlab.init(mode="local", project="my_project")
        >>> swanlab.log({"loss": 0.5})
        >>> swanlab.finish()

        Cloud run with configuration:

        >>> import swanlab
        >>> swanlab.login(api_key="your_key")
        >>> swanlab.init(
        ...     mode="cloud",
        ...     project="image_classification",
        ...     name="resnet50_experiment",
        ...     config={"lr": 0.001, "batch_size": 32}
        ... )
        >>> swanlab.log({"accuracy": 0.95})
        >>> swanlab.finish()
    """
    ...

def finish(state: FinishType = "success", error: Optional[str] = None) -> None:
    """Finish the current run and close the experiment.

    This function safely closes the current run and waits for all logs to be flushed.
    SwanLab automatically calls this function at program exit, but you can call it
    manually to mark the experiment as completed with a specific state.

    :param state: Final state of the run. Must be one of: "success", "crashed", "aborted".
        Defaults to "success".
    :param error: Error message if state is "crashed". Required when state="crashed".
    :raises RuntimeError: If called without an active run.

    Examples:

        Finish a successful run:

        >>> import swanlab
        >>> swanlab.init(mode="local")
        >>> swanlab.log({"loss": 0.5})
        >>> swanlab.finish()

        Mark run as crashed with error message:

        >>> import swanlab
        >>> swanlab.init(mode="local")
        >>> try:
        ...     raise ValueError("Training failed")
        ... except Exception as e:
        ...     swanlab.finish(state="crashed", error=str(e))
    """
    ...

def login(
    api_key: Optional[str] = None,
    relogin: bool = False,
    host: Optional[str] = None,
    save: bool = False,
    timeout: int = 10,
) -> bool:
    """Authenticate with SwanLab Cloud.

    This function authenticates your environment with SwanLab. If already logged in
    and `relogin` is False, this function does nothing. Call this before `swanlab.init()`
    to use cloud features.

    :param api_key: Your SwanLab API key. If not provided, will attempt to read from
        environment or prompt for input.
    :param relogin: If True, forces re-authentication and overwrites existing credentials.
        Defaults to False.
    :param host: Custom API host URL. If not provided, uses the default SwanLab cloud host.
    :param save: Whether to save the API key locally for future sessions. Defaults to False.
    :param timeout: Network request timeout in seconds. Defaults to 10.
    :return: True if login was successful, False otherwise.
    :raises RuntimeError: If called while a run is active.
    :raises AuthenticationError: If login fails due to invalid credentials or network issues.

    Examples:

        Login with an API key:

        >>> import swanlab
        >>> swanlab.login(api_key="your_api_key_here")
        >>> swanlab.init(mode="cloud")

        Force re-login and save credentials:

        >>> import swanlab
        >>> swanlab.login(api_key="new_api_key", relogin=True, save=True)
    """
    ...

def merge_settings(settings: Union[Settings, dict]) -> None:
    """Merge custom settings into the global SwanLab configuration.

    This function allows you to customize SwanLab's behavior before initializing a run.
    It must be called before `swanlab.init()`.

    :param settings: Custom settings to merge. Can be either a Settings object or a dict.
    :raises RuntimeError: If called while a run is active.

    Examples:

        >>> import swanlab
        >>> swanlab.merge_settings({"mode": "local", "logdir": "./my_logs"})
        >>> swanlab.init()
    """
    ...

# ── run access ─────────────────────────────────────────────────────────────────

run: Optional[Run]
"""The current active SwanLab run, or None if no run is active."""

def has_run() -> bool:
    """Check if there is an active SwanLab run.

    :return: True if a run is currently active, False otherwise.

    Examples:

        >>> import swanlab
        >>> if swanlab.has_run():
        ...     swanlab.log({"metric": 1.0})
        ... else:
        ...     print("No active run")
    """
    ...

def get_run() -> Run:
    """Get the current active SwanLab run.

    :return: The active Run instance.
    :raises RuntimeError: If no run is currently active.

    Examples:

        >>> import swanlab
        >>> swanlab.init(mode="local")
        >>> run = swanlab.get_run()
        >>> print(run.id)
        >>> swanlab.finish()
    """
    ...

# ── logging ────────────────────────────────────────────────────────────────────

def log(data: Mapping[str, Any], step: Optional[int] = None) -> None:
    """Log metrics and data to the current run.

    :param data: Dictionary of metric names and values to log.
    :param step: Optional step number. If not provided, auto-increments.
    :raises RuntimeError: If called without an active run.

    Examples:

        Log multiple metrics:

        >>> import swanlab
        >>> swanlab.init(mode="local")
        >>> swanlab.log({"loss": 0.5, "accuracy": 0.95})
        >>> swanlab.finish()

        Log with explicit step:

        >>> import swanlab
        >>> swanlab.init(mode="local")
        >>> swanlab.log({"loss": 0.5}, step=10)
        >>> swanlab.finish()
    """
    ...

def log_text(
    *,
    key: str,
    data: TextDatasType,
    caption: CaptionsType = None,
    step: Optional[int] = None,
) -> None:
    """A syntactic sugar for logging text data.

    :param key: The key for the text data.
    :param data: The text data itself or a Text object.
    :param caption: Optional caption for the text data.
    :param step: Optional step number. If not provided, auto-increments.
    :raises RuntimeError: If called without an active run.

    Examples:

        Log simple text:

        >>> import swanlab
        >>> swanlab.init(mode="local")
        >>> swanlab.log_text(key="output", data="Training started")
        >>> swanlab.finish()

        Log text with caption:

        >>> import swanlab
        >>> swanlab.init(mode="local")
        >>> swanlab.log_text(key="prediction", data="cat", caption="Model output")
        >>> swanlab.finish()
    """
    ...

def log_image(
    *,
    key: str,
    data: ImageDatasType,
    mode: ImageModesType = None,
    caption: CaptionsType = None,
    file_type: ImageFilesType = None,
    size: ImageSizesType = None,
    step: Optional[int] = None,
) -> None:
    """A syntactic sugar for logging image data.

    :param key: The key for the image data.
    :param data: The image data itself or an Image object.
    :param mode: PIL mode applied when converting to PIL.Image (e.g. 'RGB', 'L').
    :param caption: Optional caption for the image data.
    :param file_type: Output file format. One of ['png', 'jpg', 'jpeg', 'bmp']. Defaults to 'png'.
    :param size: Resize policy (int for max side, (w, h) tuple for exact size).
    :param step: Optional step number. If not provided, auto-increments.
    :raises RuntimeError: If called without an active run.

    Examples:

        >>> import swanlab, numpy as np
        >>> swanlab.init(mode="local")
        >>> img = np.zeros((64, 64, 3), dtype=np.uint8)
        >>> swanlab.log_image(key="sample", data=img, caption="blank image")
        >>> swanlab.finish()
    """
    ...

def log_audio(
    *,
    key: str,
    data: AudioDatasType,
    sample_rate: AudioRatesType = 44100,
    caption: CaptionsType = None,
    step: Optional[int] = None,
) -> None:
    """A syntactic sugar for logging audio data.

    :param key: The key for the audio data.
    :param data: The audio data itself or an Audio object.
    :param sample_rate: Sample rate of the audio (used when data is raw numpy array).
    :param caption: Optional caption for the audio data.
    :param step: Optional step number. If not provided, auto-increments.
    :raises RuntimeError: If called without an active run.

    Examples:

        >>> import swanlab, numpy as np
        >>> swanlab.init(mode="local")
        >>> audio = np.zeros((1, 44100), dtype=np.float32)
        >>> swanlab.log_audio(key="sound", data=audio, sample_rate=44100)
        >>> swanlab.finish()
    """
    ...

def log_video(
    *,
    key: str,
    data: VideoDatasType,
    caption: CaptionsType = None,
    step: Optional[int] = None,
) -> None:
    """A syntactic sugar for logging video data.

    Currently supported formats: GIF.

    :param key: The key for the video data.
    :param data: The video data itself or a Video object.
    :param caption: Optional caption for the video data.
    :param step: Optional step number. If not provided, auto-increments.
    :raises RuntimeError: If called without an active run.

    Examples:

        >>> import swanlab
        >>> swanlab.init(mode="local")
        >>> with open("animation.gif", "rb") as f:
        ...     swanlab.log_video(key="rollout", data=f.read())
        >>> swanlab.finish()
    """
    ...

def async_log(
    func: Callable,
    *args: Any,
    step: Optional[int] = None,
    mode: AsyncLogType = "threading",
    **kwargs: Any,
) -> Future:
    """Asynchronously execute a function and automatically log its return value.

    ``func`` is submitted to a background thread, process, or the asyncio event loop
    (depending on *mode*).  When it completes, its return value — a ``dict`` — is
    passed to :func:`log` automatically.  The call returns a
    :class:`~concurrent.futures.Future` immediately.

    :func:`finish` waits for all outstanding ``async_log`` tasks before flushing,
    so no data is lost.

    :param func: A callable returning a ``dict`` suitable for :func:`log`.
    :param args: Positional arguments forwarded to *func*.
    :param step: Optional step index.  If ``None``, auto-incremented when the task
        **completes** (not when submitted).  Pass an explicit value if step ordering matters.
    :param mode: Execution mode:

        - ``"asyncio"`` — schedule on the running asyncio event loop.  *func* must be a
          coroutine (``async def``).  Raises :exc:`RuntimeError` if no loop is running.

        - ``"threading"`` (default) — background thread.  *func* can access
          ``swanlab.config`` and return media objects (:class:`Image`, :class:`Audio`, etc.).
          Subject to the GIL.

        - ``"spawn"`` — new child process (``mp_context=spawn``).  Bypasses the GIL, ideal
          for CPU-bound work.  *func*, its arguments, and its return value **must be
          pickle-serializable** (no :class:`Image`, ``torch.Tensor``, etc.).

        - ``"fork"`` — **reserved**.  Will be enabled after ``swanlab-core`` ships.

    :param kwargs: Keyword arguments forwarded to *func*.
    :return: A :class:`~concurrent.futures.Future`.
    :raises RuntimeError: No active Run, or no asyncio event loop (``"asyncio"`` mode only).

    Examples:

        Asyncio mode — coroutine function for IO-bound work:

        >>> import swanlab
        >>> run = swanlab.init()
        >>> async def slow_compute():
        ...     import asyncio
        ...     await asyncio.sleep(2)
        ...     return {"score": 0.95}
        >>> future = swanlab.async_log(slow_compute, step=1, mode="asyncio")

        Threading mode (default) — IO-bound or returning media objects:

        >>> def fetch_score():
        ...     import time
        ...     time.sleep(2)
        ...     return {"score": 0.95}
        >>> future = swanlab.async_log(fetch_score, step=1)

        Spawn mode — CPU-bound, pickle-safe return values:

        >>> def compute_loss():
        ...     return {"loss": 0.123, "acc": 0.95}
        >>> future = swanlab.async_log(compute_loss, step=2, mode="spawn")
    """
    ...

def define_scalar(
    *,
    key: str,
    name: Optional[str] = None,
    color: Optional[str] = None,
    x_axis: Optional[ScalarXAxisType] = None,
    chart_name: Optional[str] = None,
) -> None:
    """Explicitly define a scalar column.

    Call this before logging to customize how a scalar metric is displayed,
    such as setting a display name, color, or x-axis type.

    :param key: The key for the scalar column. Supports glob patterns (e.g. "train/*") to match multiple columns at once.
    :param name: Optional display name for the scalar column.
    :param color: Optional hex color for the scalar line in charts.
    :param x_axis: Optional x-axis type. One of "_step", "_relative_time", or a custom key.
    :param chart_name: Optional name for the chart group this column belongs to.
    :raises RuntimeError: If called without an active run.

    Examples:

        >>> import swanlab
        >>> swanlab.init(mode="local")
        >>> swanlab.define_scalar("loss", color="#FF5733", x_axis="_step")
        >>> swanlab.log({"loss": 0.5})
        >>> swanlab.finish()
    """
    ...

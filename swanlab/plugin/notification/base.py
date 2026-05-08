from abc import abstractmethod
from pathlib import Path
from typing import Dict, Optional

from swanlab.sdk.internal.core_python.pkg.executor import SafeThreadPoolExecutor
from swanlab.sdk.internal.pkg.helper import fmt_run_path
from swanlab.sdk.internal.settings import Settings
from swanlab.sdk.protocol import Callback

# ---------------------------------------------------------------------------
# Message templates
# ---------------------------------------------------------------------------

_TEMPLATES: Dict[str, Dict[str, str]] = {
    "en": {
        "title": "SwanLab Message Notification\n",
        "msg_success": "SwanLab | Your experiment completed successfully\n",
        "msg_error": "Your SwanLab experiment encountered an error: {error}\n",
        "link_text": (
            "Project: {project}\n"
            "Workspace: {workspace}\n"
            "Name: {exp_name}\n"
            "Description: {description}\n"
            "Experiment Link: {link}"
        ),
    },
    "zh": {
        "title": "SwanLab 消息通知\n",
        "msg_success": "SwanLab | 您的实验已成功完成\n",
        "msg_error": "您的 SwanLab 实验遇到错误: {error}\n",
        "link_text": (
            "项目: {project}\n工作区: {workspace}\n实验名: {exp_name}\n描述: {description}\n实验链接: {link}"
        ),
    },
}


# ---------------------------------------------------------------------------
# URL builder
# ---------------------------------------------------------------------------


def build_run_url(settings: Settings, path: Optional[str]) -> Optional[str]:
    if path and settings.mode in ("online", "offline"):
        return f"{settings.web_host}{fmt_run_path(path)}"
    return None


# ---------------------------------------------------------------------------
# Base class
# ---------------------------------------------------------------------------


class NotificationCallback(Callback):
    """Base class for notification callbacks.

    Uses a shared ``SafeThreadPoolExecutor`` so that ``on_run_finished``
    returns immediately — the actual HTTP/SMTP call happens in a background
    thread.  During Python shutdown ``run()`` degrades to synchronous
    execution so no notification is lost.
    """

    _executor = SafeThreadPoolExecutor(max_workers=1, thread_name_prefix="swanlab-plugin-notification")

    def __init__(self, language: str = "zh"):
        self.language = language
        self._settings: Optional[Settings] = None
        self._path: Optional[str] = None

    @property
    def name(self) -> str:
        return self.__class__.__name__

    # -- Callback hooks -----------------------------------------------------

    def on_run_initialized(self, run_dir: Path, path: str, *args, **kwargs) -> None:
        settings = kwargs.get("settings")
        if isinstance(settings, Settings):
            self._settings = settings
            self._path = path

    def on_run_finished(self, state: str, error: Optional[str] = None, **kwargs) -> None:
        if self._settings is not None:
            self._executor.run(self._send_notification, state, error)

    # -- Template helpers ---------------------------------------------------

    def _get_url(self) -> Optional[str]:
        if self._settings is None:
            return None
        return build_run_url(self._settings, self._path)

    def _build_content(self, state: str, error: Optional[str]) -> str:
        tpl = _TEMPLATES.get(self.language, _TEMPLATES["zh"])
        content: str = tpl["title"]
        if error:
            content += tpl["msg_error"].format(error=error)
        else:
            content += tpl["msg_success"]
        url = self._get_url()
        if url and self._settings:
            content += tpl["link_text"].format(
                project=self._settings.project.name or "",
                workspace=self._settings.project.workspace or "",
                exp_name=self._settings.experiment.name or "",
                description=self._settings.experiment.description or "",
                link=url,
            )
        return content

    # -- Extension point ----------------------------------------------------

    @abstractmethod
    def _send_notification(self, state: str, error: Optional[str]) -> None:
        """Send the notification.  Runs inside the shared thread pool."""
        ...

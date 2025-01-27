from swankit.core import SwanLabSharedSettings

from swanlab.api import get_http
from swanlab.env import get_mode


class SwanlabCloudConfig:
    """
    public data for the SwanLab project when running in cloud mode.
    """

    def __init__(self):
        self.__http = None
        if get_mode() == "cloud":
            try:
                self.__http = get_http()
            except ValueError:
                pass
        self.__available = self.__http is not None

    def __get_property_from_http(self, name: str):
        """
        Get the property from the http object.
        if the http object is None, it will be initialized.
        if initialization fails, it will return None.
        """
        if self.available:
            return getattr(self.__http, name)
        return None

    @property
    def available(self):
        """
        Whether the SwanLab is running in cloud mode.
        """
        return self.__available

    @property
    def project_name(self):
        """
        The name of the project. Equal to `run.public.project_name`.
        If swanlab is not running in cloud mode, it will return None.
        """
        return self.__get_property_from_http("projname")

    @property
    def project_url(self):
        """
        The url of the project. It is the url of the project page on the SwanLab.
        If swanlab is not running in cloud mode, it will return None.
        """
        if not self.available:
            return None
        return self.__get_property_from_http("web_proj_url")

    @property
    def experiment_name(self):
        """
        The name of the experiment. It may be different from the name of swanboard.
        """
        return self.__get_property_from_http("expname")

    @property
    def experiment_url(self):
        """
        The url of the experiment. It is the url of the experiment page on the SwanLab.
        """
        if not self.available:
            return None
        return self.__get_property_from_http("web_exp_url")


class SwanLabPublicConfig:
    """
    Public data for the SwanLab project.
    """

    def __init__(self, project_name: str, settings: SwanLabSharedSettings):
        self.__project_name = project_name
        self.__cloud = SwanlabCloudConfig()
        self.__settings = settings

    def json(self):
        """
        Return a dict of the public config.
        This method is used to serialize the public config to json.
        """
        return {
            "project_name": self.project_name,
            "version": self.version,
            "run_id": self.run_id,
            "swanlog_dir": self.swanlog_dir,
            "run_dir": self.run_dir,
            "cloud": {
                "project_name": self.cloud.project_name,
                "project_url": self.cloud.project_url,
                "experiment_name": self.cloud.experiment_name,
                "experiment_url": self.cloud.experiment_url,
            },
        }

    @property
    def cloud(self):
        """
        The cloud configuration.
        """
        return self.__cloud

    @property
    def project_name(self):
        """
        The name of the project. Equal to `run.public.project_name`.
        """
        return self.__project_name

    # ---------------------------------- 继承settings的属性 ----------------------------------

    @property
    def version(self) -> str:
        """
        The version of the SwanLab.
        """
        return self.__settings.version

    @property
    def run_id(self) -> str:
        """
        The id of the run.
        """
        return self.__settings.run_id

    @property
    def swanlog_dir(self) -> str:
        """
        The directory of the SwanLab log.
        """
        return self.__settings.swanlog_dir

    @property
    def run_dir(self) -> str:
        """
        The directory of the run.
        """
        return self.__settings.run_dir

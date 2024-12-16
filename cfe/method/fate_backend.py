from abc import ABC, abstractmethod

from .._logging import logger


# Backend: abstract class, used for subsequent specific implementation such as "DockerBackend" class and "FunctionBackend" class
class Backend(ABC):
    def __init__(self, method_name):
        pass

    @abstractmethod
    def load_backend(self):
        pass


# DockerBackend: specific implementation of abstract Backend class using Docker.
class DockerBackend(Backend):
    # def __init__(self, method_name):
    #     logger.debug("DockerBackend __init__")
    #     pass

    def load_backend(self):
        # load docker image
        # check docker image exists
        pass

    def create_ti_method_container(self):
        logger.debug("create_ti_method_container")


# FunctionBackend: specific implementation of abstract Backend class using Python functions.
class FunctionBackend(Backend):
    # def __init__(self, method_name):
    #     logger.debug("FunctionBackend __init__")

    def load_backend(self):
        # load docker image
        # check pypi package exists
        pass

    def install_pipy_package(self):
        logger.debug("install_pipy_package")

from .._logging import logger
from .fate_backend import DockerBackEnd


# CfeDockerBackend: specific implementation of abstract Backend class using CFE Docker..
class CFEDockerBackend(DockerBackEnd):
    def __init__(self, method_name):
        logger.debug("CFEDockerBackend __init__")

    def load_backend(self):
        # load docker image
        pass

    def run(self):
        pass

    def _load_definition(self):
        """
        Load definition yml file during load_backend
        """
        pass

    def install_pipy_package(self):
        logger.debug("install_pipy_package")

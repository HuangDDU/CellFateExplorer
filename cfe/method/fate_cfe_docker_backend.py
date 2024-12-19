from .._logging import logger
from .fate_backend import Backend


# CfeDockerBackend: specific implementation of abstract Backend class using CFE Docker..
class CFEDockerBackend(Backend):
    def __init__(self, method_name):
        logger.debug("CFEDockerBackend __init__")

    def load_backend(self):
        # load docker image
        pass

    def preprocess(self, inputs, tmp_wd):
        pass

    def execute(self, tmp_wd):
        pass

    def postprocess(self):
        pass

    def install_pipy_package(self):
        logger.debug("install_pipy_package")

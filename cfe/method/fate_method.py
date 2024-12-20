import os.path
import yaml

from .._logging import logger
from .._settings import settings

from ..data.fate_anndata import FateAnnData

from .fate_dynverse_docker_backend import DynverseDockerBackend
from .fate_function_backend import FunctionBackend
from .fate_cfe_docker_backend import CFEDockerBackend


class FateMethod():

    def __init__(self, method_name="paga", backend=None):
        # logger.debug("FateMethod __init__")
        self.method_name = method_name
        self.choose_backend(backend)

    def choose_backend(self, backend=None):
        """
        ref: pydynverse.methods.method_choose_backend.method_choose_backend
        """
        # logger.debug("FateMethod choose_backend")
        backend = settings.backend if backend is None else backend
        if backend is None:
            # backend in function parameteres and setting file are both None, choose backend
            # input msg for choosing backend
            answer = input("""
                    You can run this method as an Python function (1), CFE Docker containter(2), Dynverse Docker container(3, default)
                    Which do you want to use?
                    1: Python function 
                    2: CFE Docker Container
                    else: Dynverse Docker Container[default]
            """)

            if answer == "1":
                backend = "python_function"
            elif answer == "2":
                backend = "cfe_docker"
            else:
                backend = "dynverse_docker"
            settings["backend"] = backend  # update default backend in setting

        with open(os.path.join(os.path.dirname(__file__), "method_backend.yml"), 'r') as file:
            method_backend_dict = yaml.safe_load(file)

        if backend == "python_function":
            function_name = method_backend_dict[self.method_name]["python_function"]
            self.method_backend = FunctionBackend(function_name)
            logger.info(f"backend: Python Function")
        elif backend == "cfe_docker":
            # TODO: file is same as python_function
            image_id = method_backend_dict[self.method_name]["cfe_docker"]
            self.method_backend = CFEDockerBackend(image_id)
            logger.info(f"backend: Python Function")
        else:
            # backend == "dynverse_docker"
            image_id = method_backend_dict[self.method_name]["dynverse_docker"]
            self.method_backend = DynverseDockerBackend(image_id=image_id)
            logger.info(f"backend: Dynverse Docker")

        self.backend = backend

    def infer_trajectory(
        self,
        fadata: FateAnnData = None,
        parameters: dict = None,
    ):
        """
        ref: 
            - pydynverse/wrap/method_execute._method_execute

        """
        # logger.debug("FateMethod infer_trajectory")
        self.method_backend.run(fadata, parameters)

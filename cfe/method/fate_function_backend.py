import os
import importlib.util
import yaml

from .._logging import logger
from ..data import FateAnnData
from .fate_backend import Backend, Definition


class FunctionBackend(Backend):
    """Specific implementation of abstract Backend class using Python functions.
    """

    def __init__(self, function_name="cf_paga"):
        """Initialize the FunctionBackend class.

        Args:
            function_name (str, optional):  name of the function backend .
        """
        logger.debug("FunctionBackend __init__")

        self.function_name = function_name
        self.load_backend()

    def load_backend(self) -> None:
        """load backend from python funct6ion
        """
        function_file_path = f"{os.path.dirname(__file__)}/function/{self.function_name}.py"
        # Load the module
        spec = importlib.util.spec_from_file_location(self.function_name, function_file_path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        # Get the function from the module
        self.function = getattr(module, self.function_name)
        logger.info(f"Loaded function: {self.function} from {function_file_path}")

        self._load_definition()

    def run(self, fadata: FateAnnData, parameters: dict) -> None:
        """call the python function to get trajectory dict

        Args:
            fadata (FateAnnData): the input FateAnnData object for trajectory inference method
            parameters (dict):  the  parameters for trajectory inference method
        """

        prior_information = self._extract_prior_information(fadata, self.definition.get_inputs_df())  # check prior information and add to fadata
        default_parameters = self.definition.get_parameters()
        if parameters is not None:
            default_parameters.update(parameters)
        parameters = default_parameters

        trajectory_dict = self.function(fadata, prior_information, parameters)

        fadata.add_trajectory_by_type(trajectory_dict)

    def _load_definition(self) -> None:
        """load definition from yaml file and ceate Definition object
        """
        definition_file_path = f"{os.path.dirname(__file__)}/definition/{self.function_name}.yml"
        with open(definition_file_path, 'r') as file:
            definition_raw = yaml.safe_load(file)

        definition = Definition(definition_raw)
        definition["run"] = {"backend": "python_function", "function_name": self.function_name}
        self.definition = definition

    def install_pipy_package(self):
        # TODO: install the relevant package from pipy
        logger.debug("install_pipy_package")

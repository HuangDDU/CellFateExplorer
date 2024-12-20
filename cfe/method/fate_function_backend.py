import os
import importlib.util
import yaml

from .._logging import logger
from .fate_backend import Backend, Definition


# FunctionBackend: specific implementation of abstract Backend class using Python functions.
class FunctionBackend(Backend):
    def __init__(self, function_name="cf_paga"):
        logger.debug("FunctionBackend __init__")

        self.function_name = function_name
        self.load_backend()

    def load_backend(self):
        function_file_path = f"{os.path.dirname(__file__)}/function/{self.function_name}.py"
        # Load the module
        spec = importlib.util.spec_from_file_location(self.function_name, function_file_path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        # Get the function from the module
        self.function = getattr(module, self.function_name)
        logger.info(f"Loaded function: {self.function} from {function_file_path}")

        self._load_definition()

    def run(self, fadata, parameters):

        priors = self._extract_priors(fadata, self.definition.get_inputs_df())  # only simple check prior information, don't use the priors
        default_parameters = self.definition.get_parameters()
        if not parameters is None:
            default_parameters.update(parameters)
        parameters = default_parameters

        self.function(fadata, parameters)  # call fadata.add_trajectory

    def _load_definition(self):
        definition_file_path = f"{os.path.dirname(__file__)}/definition/{self.function_name}.yml"
        with open(definition_file_path, 'r') as file:
            definition_raw = yaml.safe_load(file)

        definition = Definition(definition_raw)
        definition["run"] = {"backend": "python_function", "function_name": self.function_name}
        self.definition = definition

    def install_pipy_package(self):
        logger.debug("install_pipy_package")

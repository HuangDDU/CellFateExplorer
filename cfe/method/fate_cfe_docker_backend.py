import tempfile

from .._logging import logger
from .fate_backend import DockerBackEnd


# CfeDockerBackend: specific implementation of abstract Backend class using CFE Docker..
class CFEDockerBackend(DockerBackEnd):
    def __init__(self, image_id):
        logger.debug("CFEDockerBackend __init__")

        self.image_id = image_id
        self.load_backend()

    def load_backend(self):
        # load docker image
        pass

    def preprocess(self, tmp_wd):
        # TODO: save adata h5ad , prior information and parameters json file in tmp_wd dir
        tmp_wd

    def execute(self, tmp_wd):
        # TODO: CFE Docker run, save dict.pkl in tmp_wd dir, return trajectory_dict
        tmp_wd
        trajectory_dict = {}
        return trajectory_dict

    def postprocess(self, fadata, trajectory_dict):
        # TODO: save trajectory_dict
        fadata.add_trajectory_by_type(trajectory_dict)

    def run(self, fadata, parameters):
        # TODO:
        self.preprocess()

        with tempfile.TemporaryDirectory() as tmp_wd:
            self.preprocess()

            trajectory_dict = self.execute(tmp_wd)

            self.postprocess(fadata, trajectory_dict)

    def _load_definition(self):
        # TODO:
        """
        Load definition yml file during load_backend
        """
        pass

    def install_pipy_package(self):
        logger.debug("install_pipy_package")

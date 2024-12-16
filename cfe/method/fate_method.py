
from .._logging import logger
from .._settings import settings

from ..data.fate_anndata import FateAnnData


class FateMethod():

    def __init__(self, method_name="PAGA", backend_str=settings.backend):
        logger.debug("FateMethod __init__")
        self.method_name = method_name

        backend = self.choose_backend(backend_str)

    def choose_backend(self, backend_str):
        logger.debug("FateMethod choose_backend")
        return None

    def infer_trajectory(self, fadata: FateAnnData = None):
        logger.debug("FateMethod infer_trajectory")
        # fadata.add_trajectory(self)

    #     self.execute()

    # def execute(self):
    #     logger.debug("execute")
    #     self.process_definition()
    #     self.extract_args()

    # def process_definition(self):
    #     logger.debug("process_definition")

    # def extract_args(self):
    #     logger.debug("extract_args")

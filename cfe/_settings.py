from ._logging import logger


class CellFateExplorerConfig:
    def __init__(self):
        # plot settings
        self.plot_format = "pdf"

        # backend settings
        self.backend = "python_function"  # ["python_function", "cfe_docker", "dynverse_docker"]

        # check if rpy2 is available
        try:
            import rpy2
            self.r_available = True
            logger.debug("R is available. You can use dynverse backend.")
        except:
            self.r_available = False
            logger.warning("R not available. You cannot use dynverse backend.")
        # manually settings r_available=False for testing
        self.r_available = False

    def __getitem__(self, key):
        if hasattr(self, key):
            return getattr(self, key)
        else:
            return None

    def __setitem__(self, key, value):
        setattr(self, key, value)


settings = CellFateExplorerConfig()

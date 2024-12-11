class CellFateExplorerConfig:
    def __init__(self):
        # plot settings
        self.plot_format = "pdf"

        # backend settings
        self.backend = None

    def __getitem__(self, key):
        if hasattr(self, key):
            return getattr(self, key)
        else:
            return None

    def __setitem__(self, key, value):
        setattr(self, key, value)


settings = CellFateExplorerConfig()

import pytest
import cfe

import pandas as pd


class TestDynverseDockerBackend:
    def setup_method(self):
        self.dyn_docker = cfe.method.DynverseDockerBackend("dynverse/ti_slingshot:v1.0.3")

    def test_load_backend(self):
        # no progression bar show here, user permission may not be sufficient
        definition_object = self.dyn_docker.load_backend()

    def test_process_definition(self):
        pass

    def test_run(self):
        pass


if __name__ == "__main__":
    pytest.main(["-v", __file__])

import pytest
import cfe


class TestCFEDockerBackend():
    def setup_method(self):
        cfe

    def test_load_backend(self):
        pass


if __name__ == "__main__":
    pytest.main(["-v", __file__])

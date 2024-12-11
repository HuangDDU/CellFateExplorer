import pytest
import cfe

def test_settings():
    backend = cfe.settings["backend"]
    assert backend in [None, "function", "container"]


if __name__ == "__main__":
    pytest.main(["-v", __file__])
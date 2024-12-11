import pytest
import cfe

import anndata as ad

def test_fate_anndata():
    fadata = cfe.data.FateAnnData()
    assert isinstance(fadata, ad.AnnData)


if __name__ == "__main__":
    pytest.main(["-v", __file__])

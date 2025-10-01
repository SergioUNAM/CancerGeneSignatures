import pandas as pd
import pytest

from app.core.fold_change import compute_fold_change


@pytest.fixture()
def controles_muestras(qpcr_long_df):
    long_df = qpcr_long_df
    controles = long_df[long_df["test"].str.startswith("CTRL")].copy()
    muestras = long_df[long_df["test"].str.startswith("TUM")].copy()
    assert not controles.empty and not muestras.empty
    return controles, muestras


def test_compute_fold_change_returns_reference_gene(controles_muestras):
    controles, muestras = controles_muestras

    result = compute_fold_change(controles, muestras)

    assert result.reference_gene == "GAPDH"
    assert {"target", "fold_change"}.issubset(result.by_means.columns)
    assert {"target", "fold_change"}.issubset(result.by_refgene.columns)

    myc_fc = result.by_refgene.loc[result.by_refgene["target"] == "MYC", "fold_change"].iloc[0]
    assert myc_fc > 1.0  # ΔΔCt negativo → sobreexpresión relativa en MYC


def test_compute_fold_change_requires_shared_genes(controles_muestras):
    controles, muestras = controles_muestras
    muestras_sin_overlap = muestras.assign(target=muestras["target"].apply(lambda g: f"{g}_ALT"))

    with pytest.raises(ValueError, match="No hay genes en común"):
        compute_fold_change(controles, muestras_sin_overlap)

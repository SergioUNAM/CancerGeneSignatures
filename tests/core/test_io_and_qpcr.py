import pandas as pd


def test_parse_qpcr_wide_handles_duplicates_and_tokens(qpcr_load_result):
    load_result = qpcr_load_result

    assert list(load_result.df.columns[:2]) == ["Well", "Target Name"]
    assert load_result.meta["sample_names"] == ["CTRL-01", "CTRL-01_2", "TUM-01", "TUM-01_2"]

    # Undetermined token becomes NaN and meta captures the record
    myc_row = load_result.df[load_result.df["Target Name"] == "MYC"].iloc[0]
    assert pd.isna(myc_row["CTRL-01"])
    undetermined = load_result.meta["undetermined_records"]
    assert undetermined and undetermined[0]["target"] == "MYC"


def test_melt_wide_to_long_produces_numeric_ct(qpcr_long_df):
    long_df = qpcr_long_df

    # Expected width: 3 genes * 4 tests = 12 rows
    assert len(long_df) == 12
    assert set(long_df.columns) == {"Well", "target", "test", "ct"}

    # Ensure numeric conversion succeeded and NaN preserved for undetermined values
    assert long_df["ct"].dtype.kind in {"f", "i"}
    myc_ctrl = long_df[(long_df["target"] == "MYC") & (long_df["test"] == "CTRL-01")]["ct"]
    assert myc_ctrl.isna().all()

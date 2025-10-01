import sys
from io import BytesIO
from pathlib import Path
from typing import Iterator

import pandas as pd
import pytest


@pytest.fixture(scope="session", autouse=True)
def _ensure_src_on_path() -> Iterator[Path]:
    """Guarantee that `src/` is importable during the test session."""
    root = Path(__file__).resolve().parents[1]
    src_path = root / "src"
    if str(src_path) not in sys.path:
        sys.path.insert(0, str(src_path))
    yield root
    # No cleanup required; pytest finalises the process afterwards.


@pytest.fixture(scope="session")
def synthetic_qpcr_workbook() -> BytesIO:
    """Create an in-memory Excel workbook emulating the qPCR template layout."""
    raw = pd.DataFrame(
        [
            ["Tipo de cancer", "Gastrico", "", "", "", ""],
            ["Metodo", "qPCR", "", "", "", ""],
            ["Contexto", "Metastasis", "", "", "", ""],
            ["Prefijo controles", "CTRL", "", "", "", ""],
            ["Prefijo muestras", "TUM", "", "", "", ""],
            ["", "", "", "", "", ""],
            ["", "", "CTRL-01", "CTRL-01", "TUM-01", "TUM-01"],
            ["", "", "CT", "CT", "CT", "CT"],
            ["Well", "Target Name", "CTRL-01", "CTRL-02", "TUM-01", "TUM-02"],
            ["A1", "GAPDH", 20.1, 20.3, 25.4, 25.2],
            ["A2", "ACTB", 22.2, 22.3, 26.0, 26.4],
            ["A3", "MYC", "Undetermined", 28.1, 31.5, 31.7],
        ]
    )
    buffer = BytesIO()
    with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
        raw.to_excel(writer, index=False, header=False)
    buffer.seek(0)
    buffer.name = "synthetic_qpcr.xlsx"
    return buffer


@pytest.fixture()
def qpcr_load_result(synthetic_qpcr_workbook):  # type: ignore[override]
    from src.core.io import parse_qpcr_wide

    workbook = BytesIO(synthetic_qpcr_workbook.getvalue())
    workbook.name = synthetic_qpcr_workbook.name
    return parse_qpcr_wide(workbook, sheet_name=0)


@pytest.fixture()
def qpcr_long_df(qpcr_load_result):  # type: ignore[override]
    from src.core.qpcr import melt_wide_to_long

    return melt_wide_to_long(qpcr_load_result.df)

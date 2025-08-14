import pytest


def test_parsetopos():
    from grins.gen_params import parse_topos
    import pandas as pd
    import pandas.testing as pdt

    # Parsing the topo file
    ts_df = parse_topos("./resources/TS.topo")
    # Reading the parquet file to crosscheck
    ts_chk_df = pd.read_parquet("./resources/TS_topoparse.parquet")

    pdt.assert_frame_equal(ts_df, ts_chk_df)

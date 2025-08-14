def test_parsetopos():
    from grins.gen_params import parse_topos
    import pandas as pd
    import pandas.testing as pdt
    import os

    # getting the current path
    curr_path = os.path.dirname(__file__)

    # Parsing the topo file
    # ts_df = parse_topos("/test/resources/TS.topo")
    ts_df = parse_topos(os.path.join(curr_path, "resources", "TS.topo"))
    # Reading the parquet file to crosscheck
    ts_chk_df = pd.read_parquet(
        os.path.join(curr_path, "resources", "TS_topoparse.parquet")
    )

    pdt.assert_frame_equal(ts_df, ts_chk_df)

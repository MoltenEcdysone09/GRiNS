import pytest


def test_racipesim():
    from grins.racipe_run import topo_simulate
    import pandas as pd
    import pandas.testing as pdt

    # Read the parameter file and initial condition file
    prs_df = pd.read_csv("./resources/TS_params.csv")
    ini_df = pd.read_csv("./resources/TS_initconds.csv")

    # Simulating the topo file
    res_df = topo_simulate(
        topo_file="./resources/TS.topo",
        replicate_dir="./resources/",
        initial_conditions=ini_df,
        parameters=prs_df,
        ode_term_dir="./resources",
    )

    # Reading the pregenerated solution file
    res_chk_df = pd.read_csv("./resources/TS_solution.csv")

    # Dropping the time column
    res_df.drop(columns=["Time"], inplace=True)
    res_chk_df.drop(columns=["Time"], inplace=True)

    # As res_df is float32 cause of GPU convert into float64
    res_df = res_df.astype("float64")

    # Check if the values of the simulation match
    pdt.assert_frame_equal(res_df, res_chk_df, rtol=1e-2, atol=1e-2)

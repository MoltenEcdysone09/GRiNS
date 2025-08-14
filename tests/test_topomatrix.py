def test_topomatrix():
    from grins.ising_bool import parse_topo_to_matrix
    import numpy as np
    import numpy.testing as npt
    import os

    curr_dir = os.path.dirname(__file__)

    # Read the TS matrix
    ts_mat, node_list = parse_topo_to_matrix(
        os.path.join(curr_dir, "resources", "TS.topo")
    )

    # Check the adjacnecny matrix against the TS adjmat
    ts_chk_mat = np.array([[0, -1], [-1, 0]])

    npt.assert_array_equal(ts_mat, ts_chk_mat)

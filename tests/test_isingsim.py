import pytest


def test_isingsim():
    from grins.ising_bool import parse_topo_to_matrix, simulate_sync_trajectory
    import numpy as np
    import numpy.testing as npt
    import jax.numpy as jnp

    # Getting the adjmat of the TS topo
    ts_adj, ts_nodes = parse_topo_to_matrix("./resources/TS.topo")

    # Simulating the topo file to get the steady state
    ts_sim = simulate_sync_trajectory(
        jnp.array([-1, 1], dtype=jnp.int16), ts_adj, jnp.array([-1, 1]), jnp.arange(3)
    )

    npt.assert_array_equal(np.array([0, 1]), np.array(ts_sim[-1, 1:]))

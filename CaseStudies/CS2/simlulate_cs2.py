from grins import racipe_run as racipe
from grins import ising_bool as ising
import os
import jax.numpy as jnp


# Save Directory
save_dir = "CS2_SimulationResults"
os.makedirs(save_dir, exist_ok=True)

# Simulation Parameters
num_replicates = 3
num_params = 10000
num_init_conds = 100

# Generating the Parameters and Intial Conditions
racipe.gen_topo_param_files(
    "EMT22N.topo",
    save_dir,
    num_replicates,
    num_params,
    num_init_conds,
    sampling_method="Sobol",
)

# Simulation
racipe.run_all_replicates(
    "EMT22N.topo",
    save_dir,
    batch_size=1000,
)

# Running Ising Simulations
ising.run_all_replicates_ising(
    "EMT22N.topo",
    num_initial_conditions=1000,
    max_steps=1000,
    batch_size=100000,
    mode="sync",
    packbits=False,
    replacement_values=jnp.array([-1, 1]),
)

ising.run_all_replicates_ising(
    "EMT22N.topo",
    num_initial_conditions=1000,
    max_steps=1000,
    batch_size=100000,
    mode="async",
    packbits=False,
    replacement_values=jnp.array([-1, 1]),
)

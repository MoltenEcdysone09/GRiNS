from grins import racipe_run as racipe
import os
import glob

# Save Directory
save_dir = "CS1_SimulationResults"
os.makedirs(save_dir, exist_ok=True)

# Specifying the folder with the topo files
topos_folder = "TOPOS"
# Listing all the topo files
topo_file_list = sorted(glob.glob(os.path.join(topos_folder, "*.topo")))
print(f"Number of GRNs to simulate: {len(topo_file_list)}")

# Simulation Parameters
num_replicates = 3
num_params = 10000
num_init_conds = 1000


for topo_file in topo_file_list:
    # Generating the Parameters and Intial Conditions
    racipe.gen_topo_param_files(
        topo_file,
        save_dir,
        num_replicates,
        num_params,
        num_init_conds,
        sampling_method="Sobol",
    )
    # Simulation
    racipe.run_all_replicates(
        topo_file, save_dir, batch_size=10000, normalize=True, discretize=True
    )

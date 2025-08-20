from grins import racipe_run as racipe
import os
import glob
from itertools import product
import time
import pandas as pd

# Save Directory
save_dir = "Benchmark_SimulationResults"
os.makedirs(save_dir, exist_ok=True)

# Specifying the folder with the topo files
topos_folder = "Benchmark_TOPOS"
# Listing all the topo files
topo_file_list = sorted(glob.glob(os.path.join(topos_folder, "*.topo")))
print(f"Number of GRNs to simulate: {len(topo_file_list)}")

# Simulation Parameters
num_replicates = 1

# Generating the combinations of the simulation parameters
num_param_inicond_comb = list(product([10, 100, 1000, 10000], [10, 100, 1000, 10000]))


for rpl in range(1, 4):
    # List to store the time elapsed
    time_df = []
    # STart a file to log teh times
    with open(f"benchmark_time_{rpl}.csv", "w") as bnch_file:
        bnch_file.write(
            ",".join(["NetName", "NumParams", "NumInitConds", "Time"]) + "\n"
        )
        for num_prs, num_ics in num_param_inicond_comb:
            # Make a save directory
            prsic_save_dir = f"{num_prs:06}P_{num_ics:06}I"
            # Loop through the topo files and simulate
            for topo_file in topo_file_list:
                # Topo Name extraction
                topo_name = topo_file.split("/")[-1].replace(".topo", "")
                # Batch Size
                if topo_name == "EMT22N":
                    batch_size = 1000
                else:
                    batch_size = 10000
                # Result folder name
                res_dir = os.path.join(
                    save_dir, prsic_save_dir, f"{num_prs:06}P_{num_ics:06}I_{topo_name}"
                )
                # Start timer
                start_time = time.time()
                # Generating the Parameters and Intial Conditions
                racipe.gen_topo_param_files(
                    topo_file,
                    res_dir,
                    num_replicates,
                    num_prs,
                    num_ics,
                    sampling_method="Sobol",
                )
                # Simulation
                racipe.run_all_replicates(
                    topo_file,
                    res_dir,
                    batch_size=batch_size,
                    normalize=False,
                    discretize=False,
                )
                # Append the time to the time dataframe
                time_df.append([topo_name, num_prs, num_ics, time.time() - start_time])
                # Add a line to the csv
                bnch_file.write(
                    ",".join(
                        [
                            topo_name,
                            str(num_prs),
                            str(num_ics),
                            str(time.time() - start_time),
                        ]
                    )
                    + "\n"
                )
        #     break
        # break

    # Converting the time list into pandas dataframe
    time_df = pd.DataFrame(
        time_df, columns=["NetName", "NumParams", "NumInitConds", "Time"]
    )
    time_df.to_csv(f"{save_dir}/Benchmark_{rpl}.csv", index=False)
    print(time_df)

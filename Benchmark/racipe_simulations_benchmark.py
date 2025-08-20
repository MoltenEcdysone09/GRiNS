import os
import glob
from itertools import product
import time
import pandas as pd
import shutil
import subprocess

# Get the current working directry
cwd = os.getcwd()

# Save Directory
save_dir = "Benchmark_RACIPE_SimulationResults"
os.makedirs(save_dir, exist_ok=True)

# Specifying the folder with the topo files
topos_folder = "Benchmark_TOPOS"
# Listing all the topo files
topo_file_list = sorted(glob.glob(os.path.join(topos_folder, "*.topo")))
print(f"Number of GRNs to simulate: {len(topo_file_list)}")

# Simulation Parameters
num_replicates = 1

# Generating the combinations of the simulation parameters
num_param_inicond_comb = list(product([10, 100, 1000], [10, 100, 1000]))
# num_param_inicond_comb = list(product([10], [10]))

for rpl in range(1, 4):
    # List to store the time elapsed
    time_df = []
    # STart a file to log teh times
    with open(f"racipe_benchmark_time_{rpl}.csv", "w") as bnch_file:
        bnch_file.write(
            ",".join(["NetName", "NumParams", "NumInitConds", "Time"]) + "\n"
        )
        for num_prs, num_ics in num_param_inicond_comb:
            # Make a save directory
            prsic_save_dir = f"{num_prs:06}P_{num_ics:06}I_RACIPE"
            # Loop through the topo files and simulate
            for topo_file in topo_file_list:
                # Topo Name extraction
                topo_name = topo_file.split("/")[-1].replace(".topo", "")
                # Result folder name
                res_dir = os.path.join(
                    save_dir,
                    prsic_save_dir,
                    f"{topo_name}_RACIPE",
                )
                # Making the directory
                os.makedirs(res_dir, exist_ok=True)
                # Copy the topo file to the result directory
                shutil.copy(topo_file, res_dir)
                # # Change the directory to the directory
                # os.chdir(res_dir)
                # print(os.getcwd())
                # Start timer
                start_time = time.time()
                # Running the RACIPE command
                subprocess.run(
                    [
                        "zsh",
                        "-i",
                        "-c",
                        f"cd '/home/csb/Pradyumna/Benchmark/{res_dir}' && OLD_RACIPE '{topo_name}.topo' -num_paras {num_prs} -num_ode {num_ics} && wait",
                    ]
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
    # Converting the time list into pandas dataframe
    time_df = pd.DataFrame(
        time_df, columns=["NetName", "NumParams", "NumInitConds", "Time"]
    )
    time_df.to_csv(f"{save_dir}/RACIPE_Benchmark_{rpl}.csv", index=False)
    print(time_df)
    # break

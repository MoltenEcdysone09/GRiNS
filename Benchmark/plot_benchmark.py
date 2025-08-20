import pandas as pd
import os
import glob
import matplotlib

matplotlib.use("Qt5Agg")
import seaborn as sns
from matplotlib import pyplot as plt

# # Listing out the benchmark files from the benchmarks directory
# benchmark_dir = "./GRiNS_Benchmark/"
# benchmark_files = glob.glob(f"{benchmark_dir}/*.csv")
#
# print(benchmark_files)
#
# # Combined benchmark dataframe
# bench_df = []
#
# # Lopping through the benchmark files
# for bfl in benchmark_files:
#     tmp_df = pd.read_csv(bfl)
#     # If racipe in bfl name add RACIPE as the simulation method
#     if "RACIPE" in bfl:
#         tmp_df["SimulMethod"] = "RACIPE"
#     else:
#         tmp_df["SimulMethod"] = "GRiNS"
#     # Getting the replicate number
#     rep_num = bfl.split("_")[-1].rstrip(".csv")
#     # Adding a coulmn for replicate number
#     tmp_df["Replicate"] = rep_num
#     # Remove all the rows which have more than 1000 params or initconds as racipe does not have them
#     tmp_df = tmp_df[(tmp_df["NumParams"] <= 1000) & (tmp_df["NumInitConds"] <= 1000)]
#     # COmbined combinations
#     tmp_df["NumParamInitCondCombs"] = tmp_df["NumParams"] * tmp_df["NumInitConds"]
#     # Append to the list
#     bench_df.append(tmp_df)
#
# # Merge all the dataframes
# bench_df = pd.concat(bench_df)
# print(bench_df)
# # Saving the combined dataframe
# bench_df.to_csv("./benchmark_combined.csv", index=False)

# Reading the combined benchmark file
bench_df = pd.read_csv("./benchmark_combined.csv")
print(bench_df)

# PLotting the combined benchmarks for each network
for net in bench_df["NetName"].unique():
    print(net)
    tmpdf = bench_df[bench_df["NetName"] == net]
    print(tmpdf)
    # Plot with standard deviation as error bars
    sns.lineplot(
        data=tmpdf,
        x="NumParamInitCondCombs",
        y="Time",
        hue="SimulMethod",
        errorbar="sd",
        err_style="bars",
        marker="o",
        err_kws={"capsize": 2.0},
    )
    plt.xscale("log")  # optional: makes sense for large ranges
    plt.yscale("log")  # optional: time likely scales with complexity
    plt.xlabel("Number of Parameter - Initial Condition Combinations")
    plt.ylabel("Simulation Time (s)")  # adjust units if needed
    plt.title(f"Simulation Time vs Combinations for {net}")
    plt.legend(
        title="Simulation Method",  # Legend title
        bbox_to_anchor=(1.05, 1),  # Position outside the axes (to the right)
        loc="upper left",  # Anchor point of the legend
    )
    plt.tight_layout()

    plt.savefig(f"./Plots/{net}_lineplot.png", dpi=450)
    plt.savefig(f"./Plots/{net}_lineplot.svg", dpi=300)
    plt.show()
    plt.close()

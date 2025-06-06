import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

# Set font and plot styles in Matplotlib
plt.rcParams.update(
    {
        "font.family": "Roboto",
        # "font.weight": "medium",
        "legend.fontsize": 13,
        "legend.title_fontsize": 14,
        "font.size": 16,
        "axes.titlesize": 16,
        "axes.labelsize": 16,
        "xtick.labelsize": 16,
        "ytick.labelsize": 16,
        "axes.edgecolor": "black",  # For black tick borders
        "axes.linewidth": 1.5,  # Make axis lines more visible
        # "axes.spines.left": True,  # Display the left spine
        # "axes.spines.bottom": True,  # Display the bottom spine
        # "axes.spines.right": True,  # Display the right spine
        # "axes.spines.top": True,  # Display the top spine
        # "lines.linewidth": 0.8,  # Set the default line width for plots
        "xtick.major.width": 1.5,  # X-axis major ticks
        "ytick.major.width": 1.5,  # Y-axis major ticks
        "lines.linewidth": 1.5,  # Lines in line plots
        "patch.linewidth": 1.5,  # For bar edge lines (used by histplot)
        "legend.edgecolor": "black",  # frame color for the legend
        "savefig.dpi": 300,  # Figure DPI while saving
    }
)

# Apply the "ticks" style manually by removing the grid
plt.style.use("seaborn-v0_8-deep")


new_order = [
    "miR141",
    "miR101",
    "miR34a",
    "miR200a",
    "miR200c",
    "miR200b",
    "GSC",
    "TWIST1",
    "TWIST2",
    "FOXC2",
    "TGFbeta",
    "SNAI1",
    "SNAI2",
    "ZEB1",
    "ZEB2",
]

# Specifying the results directory
result_dir = "IsingSimulResults"

# Specifying the topo name
topo_name = "EMT22N"

# Lists of df to store replcaite and a/sync resutls
res_df = []

# Loop through the updates
for updatever in ["sync", "async"]:
    # Loop through the replicates
    for replicate in range(1, 4):
        ising_df = pd.read_parquet(
            f"./{result_dir}/{topo_name}/{replicate}/{topo_name}_{updatever}_ising_results.parquet"
        )
        # Get the max step value, as things start from 0 it will be the max value of step column + 1
        print(ising_df["Step"].max() + 1)
        # Adding intial condition numbers
        ising_df["InitNum"] = ising_df.index // ising_df["Step"].max() + 1
        ising_df["State"] = ising_df[new_order].astype(str).agg("".join, axis=1)
        ising_df = ising_df[["Step", "State", "InitNum"]]

        # Subsetting the to keep the end states
        stg_df = pd.merge(
            ising_df[ising_df["Step"] == 999],
            ising_df[ising_df["Step"] == 1000],
            on="InitNum",
            suffixes=("_Prev", "_Last"),
        )

        stg_df = stg_df[["InitNum", "State_Prev", "State_Last"]]
        stg_df = stg_df[stg_df["State_Prev"] == stg_df["State_Last"]]
        stg_df = stg_df[["InitNum", "State_Last"]]
        stg_df["SteadyState"] = "'" + stg_df["State_Last"] + "'"
        stg_df["Replicate"] = replicate
        stg_df = stg_df[["InitNum", "SteadyState", "Replicate"]]
        stg_df["Mode"] = updatever.capitalize()
        print(stg_df)
        stg_df.to_csv(
            f"./{result_dir}/{topo_name}/{replicate}/{topo_name}_{updatever}_SteadyStates.csv",
            index=False,
        )

        state_counts = stg_df["SteadyState"].value_counts(normalize=True).reset_index()
        state_counts.columns = ["SteadyState", "Fraction"]
        state_counts["Replicate"] = replicate
        state_counts["Mode"] = updatever.capitalize()
        print(state_counts)
        # Appendind to thh results dataframe
        res_df.append(state_counts)
        state_counts.to_csv(
            f"./{result_dir}/{topo_name}/{replicate}/{topo_name}_{updatever}_StateCounts.csv",
            index=False,
        )

# Concatenate all the state count dataframes and save
res_df = pd.concat(res_df, axis=0)
res_df.to_csv(
    f"./{result_dir}/{topo_name}/{topo_name}_StateCounts_Main.csv", index=False
)

res_df = pd.read_csv(f"./{result_dir}/{topo_name}/{topo_name}_StateCounts_Main.csv")

# Filter the counts to include steady states which have a frequency higher than 1%
res_df = res_df[res_df["Fraction"] > 0.02]

# Remove the quote from strings
res_df["SteadyState"] = res_df["SteadyState"].str.strip("'")
print(res_df)

plt.figure(figsize=(5.5, 6))
ax = sns.barplot(
    data=res_df,
    x="SteadyState",
    y="Fraction",
    hue="Mode",
    errorbar="sd",  # Use standard deviation as the error bar
    capsize=0.3,
    edgecolor="black",
    linewidth=1.5,
    err_kws={"color": "black", "linewidth": 1.5},
)
plt.xlabel("Steady State")
plt.ylabel("Fraction of States")
sns.move_legend(ax, loc="upper right", title="Update Type")

# Extending the y axis to avoid text overlap
ymin, ymax = ax.get_ylim()
ax.set_ylim(ymin, ymax * 1.03)

# Add value labels above bars
for container in ax.containers:
    for bar in container:
        height = bar.get_height()
        if height > 0:  # Avoid labeling bars with 0 height
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                height + 0.02,
                f"{height:.2f}",
                ha="center",
                va="bottom",
                fontsize=13,
                color="black",
            )

# plt.xticks(rotation=90)  # Rotate x-labels for readability
# plt.title("")
plt.tight_layout()
plt.savefig("./IsingSimulResults/EMT22N/async_scync_freqside.png", dpi=300)
plt.savefig("./IsingSimulResults/EMT22N/async_scync_freqside.svg", dpi=300)

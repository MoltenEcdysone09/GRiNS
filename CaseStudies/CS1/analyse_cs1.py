import pandas as pd
import seaborn as sns
import os
import glob
from matplotlib import pyplot as plt
# import numpy as np
# import matplotlib.patches as patches
# from scipy.cluster.hierarchy import linkage, dendrogram, fcluster

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


# Define a function that will split the State column and add new columns based on the list of column names
def split_state_and_add_columns(df, state_column, column_names):
    # Split the state column into individual characters and create new columns with the 'S_' prefix
    split_columns = (
        df[state_column].apply(lambda x: pd.Series(list(x.strip("'")))).astype(int)
    )
    # Rename the new columns with the prefix "S_" and append them to the original dataframe
    split_columns.columns = [f"S_{col}" for col in column_names]
    # Add the new columns to the dataframe
    df = pd.concat([df, split_columns], axis=1)
    # Return the updated dataframe
    return df


results_dir = "CS1_SimulationResults"

# Listing the results folders
topo_result_dirs = sorted(glob.glob(results_dir + "/TS*/"))
print(topo_result_dirs)

# Loop through the results
for topo_result in topo_result_dirs:
    topo_name = os.path.basename(topo_result.rstrip("/"))
    print(topo_name)
    replicate_dirs = [
        path
        for path in glob.glob(topo_result + "/*/")
        if os.path.basename(path.rstrip("/")).isnumeric()
    ]
    # LIst of dataframe to concatenate at the end
    soldf = []
    steady_counts_df = []
    steady_type_df = []
    for rep_dir in replicate_dirs:
        rep_num = os.path.basename(rep_dir.rstrip("/"))
        print("Replicate Number :", rep_num)
        # Reading the solution dataframe
        sol_df = pd.read_parquet(
            os.path.join(
                rep_dir, f"{topo_name}_steadystate_solutions_{rep_num}.parquet"
            )
        )
        # Adding Rep number
        sol_df["RepNum"] = int(rep_num)
        soldf.append(sol_df)
        # Filter the steady states
        sol_df = sol_df[sol_df["SteadyStateFlag"] == 1]
        # reading the Steady state counts
        steady_counts = sol_df["State"].value_counts(normalize=True).reset_index()
        steady_counts.columns = ["State", "Fraction"]
        steady_counts["RepNum"] = int(rep_num)
        steady_counts_df.append(steady_counts)
        # Getting the unique number of steady states in for the parameter sets
        multi_stable_df = (
            sol_df.groupby("ParamNum")["State"]
            .nunique()
            .value_counts(normalize=True)
            .reset_index()
        )
        multi_stable_df.columns = ["Stability Type", "Fraction"]
        multi_stable_df["RepNum"] = int(rep_num)
        # Get the frequency of each unique State per ParamNum
        state_freq_df = (
            sol_df.groupby("ParamNum")["State"]
            .value_counts(normalize=True)
            .reset_index(name="Frequency")
        )
        # Optional: sort each group by alphabetical order of state (within each ParamNum)
        state_freq_df = state_freq_df.sort_values(["ParamNum", "State"]).reset_index(
            drop=True
        )
        # Getting the multi stable states too
        state_summary = (
            state_freq_df.groupby("ParamNum")["State"]
            .agg(lambda x: "'" + " ".join(sorted(set(x))).replace("'", "") + "'")
            .reset_index()
        )
        state_summary.columns = ["ParamNum", "MultiStableState"]
        state_freq_df = pd.merge(
            state_freq_df, state_summary, on="ParamNum", how="left"
        )
        state_freq_df["Stability Type"] = (
            state_freq_df["MultiStableState"].str.count(" ") + 1
        )
        state_freq_df["RepNum"] = int(rep_num)
        steady_type_df.append(state_freq_df)
    soldf = pd.concat(soldf, axis=0)
    print(sol_df)
    # State Counts Dataframe - Steady State
    steady_counts_df = pd.concat(steady_counts_df, axis=0)
    steady_counts_df["Motif"] = topo_name
    print(steady_counts_df)
    steady_counts_df.to_csv(
        os.path.join(topo_result, f"{topo_name}_StateCounts.csv"), index=False
    )
    # State Type and Multistable states data
    steady_type_df = pd.concat(steady_type_df, axis=0)
    steady_type_df["Motif"] = topo_name
    steady_type_df.to_csv(
        os.path.join(topo_result, f"{topo_name}_StateTypes.csv"), index=False
    )
    print(steady_type_df)

    # Node GK normalised columns
    gk_cols = [cl for cl in soldf.columns if "gk_" in cl]
    gene_cols = [cl.replace("gk_", "") for cl in soldf.columns if "gk_" in cl]
    # print(soldf[gk_cols])
    # Sampling a section of the dataframe for better plotting
    smpl_df = sol_df.sample(n=100000).sort_values(by=gene_cols)
    smpl_df_clustered = smpl_df.sort_values(by="State")
    # Step 5: plot
    plt.figure(figsize=(4, 6))
    ax = sns.heatmap(
        smpl_df_clustered[gk_cols],
        cmap=sns.dark_palette("#8fbcbb", reverse=False, as_cmap=True),
        yticklabels=False,
        xticklabels=gene_cols,
        cbar_kws={"shrink": 0.6, "aspect": 10},
    )
    # Access the colorbar
    colorbar = ax.collections[0].colorbar
    # Set colorbar ticks
    colorbar.set_ticks([0, 0.5, 1])
    # Set colorbar border color
    colorbar.outline.set_edgecolor("black")
    colorbar.outline.set_linewidth(1.5)  # Optional: make the border thicker
    plt.tight_layout()
    plt.savefig(os.path.join(topo_result, f"{topo_name}_GKValsHeatMap.png"))
    plt.clf()
    smpl_df_clustered = split_state_and_add_columns(
        smpl_df_clustered, "State", gene_cols
    )
    # Step 5: plot
    plt.figure(figsize=(4, 6))
    ax = sns.heatmap(
        smpl_df_clustered[["S_" + g for g in gene_cols]],
        cmap=sns.dark_palette("#69d", reverse=False, as_cmap=True),
        yticklabels=False,
        xticklabels=gene_cols,
        cbar_kws={"shrink": 0.6, "aspect": 10},
    )
    # Access the colorbar
    colorbar = ax.collections[0].colorbar
    # # Set colorbar ticks
    # colorbar.set_ticks([0, 0.5, 1])
    # Set colorbar border color
    colorbar.outline.set_edgecolor("black")
    colorbar.outline.set_linewidth(1.5)  # Optional: make the border thicker
    plt.tight_layout()
    # plt.savefig(os.path.join(topo_result, f"{topo_name}_BinValsHeatMap.svg"))
    plt.savefig(os.path.join(topo_result, f"{topo_name}_BinValsHeatMap.png"))
    plt.clf()
    # plt.show()
    # sns.histplot(
    #     smpl_df[gk_cols].values.flatten(),
    #     bins=100,
    # )
    # plt.show()
    # Create the barplot
    plt.figure(figsize=(5, 5))
    steady_counts_df["State"] = steady_counts_df["State"].str.strip("'")
    ax = sns.barplot(
        data=steady_counts_df,
        x="State",
        y="Fraction",
        estimator="mean",  # Use mean across replicates
        color="#DD8452",
        errorbar="sd",  # For Seaborn >=0.12; use ci="sd" for older versions
        err_kws={"color": "black"},  # errorbar color
        capsize=0.1,
        edgecolor="black",
    )
    plt.xticks(rotation=90)
    plt.title(f"Steady State Distribution of {topo_name}")
    plt.ylabel("Fraction")
    plt.xlabel("State")
    plt.tight_layout()
    # plt.savefig(os.path.join(topo_result, f"{topo_name}_StateCounts.svg"))
    plt.savefig(os.path.join(topo_result, f"{topo_name}_StateCounts.png"))
    plt.clf()
    # # plt.show()
    # print(steady_type_df.columns)
    multi_type_counts = pd.DataFrame(
        steady_type_df.drop_duplicates(subset=["ParamNum", "RepNum"], keep="first")
        .groupby("RepNum")["Stability Type"]
        .value_counts(normalize=True)
    ).reset_index()
    multi_type_counts.columns = ["RepNum", "Stability Type", "Fraction"]
    print(multi_type_counts)
    plt.figure(figsize=(5, 5))
    ax = sns.barplot(
        data=multi_type_counts,
        x="Stability Type",
        y="Fraction",
        estimator="mean",  # Use mean across replicates
        color="#64B5CD",
        errorbar="sd",  # For Seaborn >=0.12; use ci="sd" for older versions
        err_kws={"color": "black"},  # errorbar color
        capsize=0.1,
        edgecolor="black",
    )
    plt.title(f"Multi-Stability Type Distribution of {topo_name}")
    plt.ylabel("Fraction")
    plt.xlabel("Number of Stable States")
    plt.tight_layout()
    plt.savefig(os.path.join(topo_result, f"{topo_name}_StateType.png"))
    # plt.show()
    plt.clf()
    multi_type_state = pd.DataFrame(
        steady_type_df.drop_duplicates(subset=["ParamNum", "RepNum"], keep="first")
        .groupby("RepNum")["MultiStableState"]
        .value_counts(normalize=True)
    ).reset_index()
    multi_type_state.columns = ["RepNum", "MultiStableState", "Fraction"]
    multi_type_state["MultiStableState"] = multi_type_state[
        "MultiStableState"
    ].str.replace("'", "")
    # print(multi_type_state)
    plt.figure(figsize=(5, 5))
    ax = sns.barplot(
        data=multi_type_state,
        x="MultiStableState",
        y="Fraction",
        estimator="mean",  # Use mean across replicates
        color="#64B5CD",
        errorbar="sd",  # For Seaborn >=0.12; use ci="sd" for older versions
        err_kws={"color": "black"},  # errorbar color
        capsize=0.1,
        edgecolor="black",
    )
    plt.xticks(rotation=90)
    plt.title(f"State Distribution of {topo_name}")
    plt.ylabel("Fraction")
    plt.xlabel("State")
    plt.tight_layout()
    plt.savefig(os.path.join(topo_result, f"{topo_name}_StateMulti.png"))
    # plt.show()
    plt.clf()
    # break


# for motif in ["TS", "TT"]:
for motif in ["TS"]:
    print(motif)
    topo_result_dirs = sorted(glob.glob(results_dir + f"/{motif}*/"))
    statecount_df = []
    statetype_df = []
    # Loop through the results
    for topo_result in topo_result_dirs:
        topo_name = os.path.basename(topo_result.rstrip("/"))
        print(topo_name)
        stsc_df = pd.read_csv(os.path.join(topo_result, f"{topo_name}_StateCounts.csv"))
        stsc_df["State"] = stsc_df["State"].str.strip("'")
        statecount_df.append(stsc_df)
        stst_df = pd.read_csv(os.path.join(topo_result, f"{topo_name}_StateTypes.csv"))
        statetype_df.append(stst_df)
    statecount_df = pd.concat(statecount_df, axis=0)
    print(statecount_df)
    statetype_df = pd.concat(statetype_df, axis=0)
    print(statetype_df)
    # Ensure no duplicates per ParamNum before grouping, if needed
    freq_ratios = (
        statetype_df[statetype_df["Stability Type"] == 2]
        .groupby(["ParamNum", "RepNum", "Motif"])["Frequency"]
        .apply(lambda x: x.sort_values().iloc[0] / x.sort_values().iloc[1])
        .dropna()
        .reset_index(name="RelativeStability")
        .sort_values(by=["Motif", "RelativeStability"])
        .reset_index(drop=True)
    )

    def extract_extremes_and_median(group):
        min_row = group.loc[group["RelativeStability"].idxmin()]
        max_row = group.loc[group["RelativeStability"].idxmax()]
        median_val = group["RelativeStability"].median()

        # Find the row closest to the median
        median_row = group.iloc[
            (group["RelativeStability"] - median_val).abs().argsort().iloc[0]
        ]

        return pd.DataFrame([min_row, median_row, max_row])

    freq_ratios = (
        freq_ratios.groupby("Motif", group_keys=False)
        .apply(extract_extremes_and_median)
        .reset_index(drop=True)
    )
    # print(freq_ratios)

    for row in freq_ratios.itertuples(index=False):
        # print(row)
        # print(row.ParamNum)
        # print(row.RepNum)
        # print(row.Motif)
        # Reading the solution dataframe
        sol_df = pd.read_parquet(
            os.path.join(
                results_dir,
                f"{row.Motif}/{row.RepNum:>03}/{row.Motif}_steadystate_solutions_{row.RepNum:>03}.parquet",
            )
        )
        sol_df = sol_df[sol_df["ParamNum"] == row.ParamNum]
        print(sol_df.head())
        param_df = pd.read_parquet(
            os.path.join(
                results_dir,
                f"{row.Motif}/{row.RepNum:>03}/{row.Motif}_params_{row.RepNum:>03}.parquet",
            )
        )
        param_df = param_df[param_df["ParamNum"] == row.ParamNum]
        # Reading the Initial condition dataframe
        initcond_df = pd.read_parquet(
            os.path.join(
                results_dir,
                f"{row.Motif}/{row.RepNum:>03}/{row.Motif}_init_conds_{row.RepNum:>03}.parquet",
            )
        )
        # Get the node columns from the solution dataframe
        node_cols = [
            col.replace("Prod_", "") for col in param_df.columns if "Prod_" in col
        ]
        # Get the production and degradation columns
        prod_cols = [f"Prod_{node}" for node in node_cols]
        deg_cols = [f"Deg_{node}" for node in node_cols]
        # Get the gk columns
        gk_cols = [f"gk_{node}" for node in node_cols]
        # Compute the gk values
        gk_vals = param_df[prod_cols].values / param_df[deg_cols].values
        initcond_df[node_cols] = initcond_df[node_cols].values / gk_vals[0]
        initcond_df = pd.merge(
            initcond_df, sol_df[["State", "Time", "InitCondNum"]], on="InitCondNum"
        )
        initcond_df["State"] = initcond_df["State"].str.replace("'", "")
        steadystate_vals = sol_df.groupby("State", as_index=False).mean()[
            ["State"] + gk_cols
        ]
        steadystate_vals.columns = ["State"] + node_cols
        plt.figure(figsize=(7, 6))
        sns.scatterplot(data=initcond_df, x="A", y="B", hue="State", alpha=0.8)
        plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
        plt.tight_layout()
        plt.savefig(
            f"./{row.Motif}_{row.RepNum:>03}_ICScatter_{int(row.ParamNum)}_{int(row.RelativeStability * 100)}.png"
        )

        # Plot
        plt.figure(figsize=(6, 5))
        ax = sns.barplot(
            data=statecount_df,
            x="State",
            y="Fraction",
            hue="Motif",
            errorbar="sd",  # Use standard deviation as the error bar
            capsize=0.4,
            edgecolor="black",
            linewidth=1.5,
            err_kws={"color": "black", "linewidth": 1.5},
            palette=["#55A868", "#C44E52"],
            # palette=["#7287fd", "#dd7878"],
        )
        # Annotate each bar
        for bar in ax.patches:
            height = bar.get_height()
            if not pd.isna(height) and height > 0:
                x = bar.get_x() + bar.get_width() / 2
                offset = height * 0.06 if height > 0.1 else 0.06
                ax.annotate(
                    f"{height:.2f}",
                    (x, height + offset),
                    ha="center",
                    va="bottom",
                    fontsize=14,
                )
        ax.set_ylim(0, ax.get_ylim()[1] * 1.2)
        # Customize axes
        ax.set_xlabel("State")
        ax.set_ylabel("Fraction")
        plt.legend(title="Motif")
        plt.tight_layout()
        plt.savefig(f"./{topo_name}_StateCounts.png")
        # plt.savefig(f"./{topo_name}_StateCounts.svg")
        plt.clf()
        # Plot
        multi_type_counts = pd.DataFrame(
            statetype_df.drop_duplicates(
                subset=["ParamNum", "RepNum", "Motif"], keep="first"
            )
            .groupby(["RepNum", "Motif"])["Stability Type"]
            .value_counts(normalize=True)
        ).reset_index()
        multi_type_counts.columns = ["RepNum", "Motif", "Stability Type", "Fraction"]
        multi_type_counts = multi_type_counts[multi_type_counts["Fraction"] >= 0.1]
        print(multi_type_counts)
        plt.figure(figsize=(6, 5))
        ax = sns.barplot(
            data=multi_type_counts,
            x="Stability Type",
            y="Fraction",
            hue="Motif",
            errorbar="sd",  # Use standard deviation as the error bar
            capsize=0.4,
            edgecolor="black",
            linewidth=1.5,
            err_kws={"color": "black", "linewidth": 1.5},
            # palette=["#7287fd", "#dd7878"],
            palette=["#55A868", "#C44E52"],
        )
        # Annotate each bar
        for bar in ax.patches:
            height = bar.get_height()
            if not pd.isna(height) and height > 0:
                x = bar.get_x() + bar.get_width() / 2
                offset = height * 0.05 if height > 0 else 0.01
                ax.annotate(
                    f"{height:.2f}",
                    (x, height + offset),
                    ha="center",
                    va="bottom",
                    fontsize=14,
                )
        ax.set_ylim(0, ax.get_ylim()[1] * 1.2)
        # Customize axes
        ax.set_xlabel("Stability Type")
        ax.set_ylabel("Fraction")
        plt.legend(title="Motif")
        plt.tight_layout()
        plt.savefig(f"./{topo_name}_StateTypes.png")
        # plt.savefig(f"./{topo_name}_StateTypes.svg")
        plt.clf()
        # Plot
        multi_type_state = pd.DataFrame(
            statetype_df.drop_duplicates(
                subset=["ParamNum", "RepNum", "Motif"], keep="first"
            )
            .groupby(["RepNum", "Motif"])["MultiStableState"]
            .value_counts(normalize=True)
        ).reset_index()
        multi_type_state.columns = ["RepNum", "Motif", "MultiStableState", "Fraction"]
        multi_type_state["MultiStableState"] = multi_type_state[
            "MultiStableState"
        ].str.replace("'", "")
        multi_type_state = multi_type_state[multi_type_state["Fraction"] >= 0.05]
        print(multi_type_state)
        plt.figure(figsize=(7, 5))
        ax = sns.barplot(
            data=multi_type_state,
            x="MultiStableState",
            y="Fraction",
            hue="Motif",
            errorbar="sd",  # Use standard deviation as the error bar
            capsize=0.4,
            edgecolor="black",
            linewidth=1.5,
            err_kws={"color": "black", "linewidth": 1.5},
            # palette=["#7287fd", "#dd7878"],
            palette=["#55A868", "#C44E52"],
        )
        # Annotate each bar
        for bar in ax.patches:
            height = bar.get_height()
            if not pd.isna(height) and height > 0:
                x = bar.get_x() + bar.get_width() / 2
                offset = height * 0.05 if height > 0 else 0.01
                ax.annotate(
                    f"{height:.2f}",
                    (x, height + offset),
                    ha="center",
                    va="bottom",
                    fontsize=14,
                )
        ax.set_ylim(0, ax.get_ylim()[1] * 1.2)
        # Customize axes
        ax.set_xlabel("State")
        ax.set_ylabel("Fraction")
        plt.xticks(rotation=90)
        plt.legend(title="Motif")
        plt.tight_layout()
        plt.savefig(f"./{topo_name}_StateMulti.png")
        # plt.savefig(f"./{topo_name}_StateTypes.svg")
        plt.clf()

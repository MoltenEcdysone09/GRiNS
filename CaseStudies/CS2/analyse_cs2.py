import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import os
import glob
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

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


def run_and_plot_pca(sol_df, gene_columns, n_components=5):
    """
    Run PCA on selected gene columns and create:
    - A 2D PCA scatter plot with axes and explained variance
    - A scree plot showing variance explained by each component

    Parameters:
    - sol_df: pd.DataFrame with gene expression data
    - gene_columns: list of gene column names to use
    - n_components: number of PCA components to compute
    """
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
    epi_list = ["miR141", "miR101", "miR34a", "miR200a", "miR200c", "miR200b"]
    mes_list = [
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
    print(sol_df)
    print(sol_df.columns)
    # Get the non signalling node columns
    gk_cols = ["gk_" + cl for cl in new_order]
    gene_columns = gk_cols
    sol_df = sol_df[new_order + gk_cols]
    # 1. Standardize the data
    # X = sol_df[gk_cols].values
    # gene_columns = new_order
    X = sol_df[gene_columns].values
    # X = np.log2(X + 1)
    X_scaled = StandardScaler().fit_transform(X)

    # 2. Run PCA
    pca = PCA(n_components=n_components)
    X_pca = pca.fit_transform(X_scaled)
    explained_var = pca.explained_variance_ratio_ * 100  # in %

    # Getting the Scores of the states
    epi_mean = sol_df[["gk_" + cl for cl in epi_list]].mean(axis=1)
    mes_mean = sol_df[["gk_" + cl for cl in mes_list]].mean(axis=1)
    em_score = mes_mean - epi_mean
    em_score = (em_score - em_score.min()) / (em_score.max() - em_score.min())

    # 3. Plot 2D PCA with explained variance
    fig, ax = plt.subplots(figsize=(6, 6))
    scatter = ax.scatter(
        X_pca[:, 0],
        X_pca[:, 1],
        alpha=0.8,
        # color=sns.color_palette("deep")[1],
        c=em_score,
        cmap="coolwarm_r",
        s=12,
        linewidth=0.1,
        edgecolor="k",
        # edgecolor="#3b4252",
        vmin=em_score.min(),
        vmax=em_score.max(),
    )
    ax.axhline(0, linestyle="--", color="grey", linewidth=1.5)
    ax.axvline(0, linestyle="--", color="grey", linewidth=1.5)
    ax.set_xlabel(f"PC1 ({explained_var[0]:.2f}%)")
    ax.set_ylabel(f"PC2 ({explained_var[1]:.2f}%)")
    # ax.set_title("PCA: PC1 vs PC2")
    #
    # Plotting Colorbar
    cbar = plt.colorbar(
        scatter, ax=ax, fraction=0.046, pad=0.04
    )  # fraction & pad control size/spacing
    cbar.set_ticks([0, 0.5, 1])
    cbar.set_ticklabels(["0 (Epi)", "0.5", "1 (Mes)"])
    cbar.set_label("Normalized EMT Score (Mes - Epi)", fontsize=19)

    # # Plot loading vectors (axes)
    # loadings = pca.components_.T
    # for i, gene in enumerate(gene_columns):
    #     ax.arrow(
    #         0,
    #         0,
    #         loadings[i, 0] * 3,
    #         loadings[i, 1] * 3,
    #         color="black",
    #         alpha=1,
    #         head_width=0.05,
    #     )
    #     ax.text(
    #         loadings[i, 0] * 3.2,
    #         loadings[i, 1] * 3.2,
    #         gene,
    #         fontsize=14,
    #         ha="center",
    #         va="center",
    #         color="black",
    #     )
    #
    # ax.grid(True)
    plt.tight_layout()
    # plt.savefig("./CS2_SimulationResults/EMT22N/PCA.svg")
    plt.savefig("./CS2_SimulationResults/EMT22N/PCA.png", dpi=400)
    # plt.savefig("./CS2_SimulationResults/EMT22N/PCA.svg", dpi=300)
    # plt.savefig()
    # plt.show()
    plt.clf()
    #
    # 4. Plot scree plot of explained variance
    fig2, ax2 = plt.subplots(figsize=(6, 4))
    bars = ax2.bar(
        range(1, n_components + 1),
        explained_var,
        edgecolor="k",
        color=sns.color_palette("deep")[9],
    )
    # Format the labels: two decimal places and a percent sign
    ax2.bar_label(
        bars,
        labels=[
            f"{h:.1f}%" for h in explained_var
        ],  # <-- custom list of formatted labels
        fontsize=11,
    )
    # ax2.bar_label(ax2.containers[0], fontsize=12)
    ax2.set_xlabel("Principal Component")
    ax2.set_ylabel("Variance Explained (%)")
    ax2.set_ylim(0, max(explained_var) * 1.10)
    # ax2.set_title("Scree Plot")
    ax2.set_xticks(range(1, n_components + 1))
    plt.tight_layout()
    # plt.savefig("./CS2_SimulationResults/EMT22N/PCA_Skree.svg")
    plt.savefig("./CS2_SimulationResults/EMT22N/PCA_Skree.png", dpi=300)
    plt.savefig("./CS2_SimulationResults/EMT22N/PCA_Skree.svg", dpi=300)
    plt.clf()
    # plt.show()

    # 5. Plot loadings as a vertical bar plot (for PC1)
    loadings = pca.components_.T  # (n_genes, n_components)

    fig3, ax3 = plt.subplots(figsize=(4, 6))  # Tall vertical figure

    pc_idx = 0  # 0 for PC1, 1 for PC2, etc.
    genes = [g.replace("gk_", "") for g in gene_columns]
    values = loadings[:, pc_idx]

    # Prepare colors: miRs = red, others = blue
    palette = sns.color_palette("deep")
    red_color = palette[3]  # usually a nice red
    blue_color = palette[0]  # usually a nice blue
    colors = [red_color if "miR" in gene else blue_color for gene in genes]

    # Now plot
    ax3.barh(genes, values, color=colors, edgecolor="k")
    ax3.axvline(0, color="black")

    ax3.set_ylabel("Genes")
    ax3.set_xlabel(f"Loadings on PC{pc_idx + 1}")
    # ax3.set_title(f"Gene Loadings on Principal Component {pc_idx + 1}", fontsize=16)
    ax3.tick_params(axis="both", which="major")
    plt.gca().invert_yaxis()  # So that top genes are at the top
    plt.tight_layout()
    plt.savefig("./CS2_SimulationResults/EMT22N/PCA_PC1Loadings.png", dpi=400)
    plt.savefig("./CS2_SimulationResults/EMT22N/PCA_PC1Loadings.svg", dpi=400)
    plt.clf()


results_dir = "CS2_SimulationResults"

# Listing the results folders
topo_result_dirs = sorted(glob.glob(results_dir + "/*/"))
print(topo_result_dirs)

# Loop through the results
for topo_result in topo_result_dirs:
    topo_name = os.path.basename(topo_result.rstrip("/")).replace("_BAK", "")
    print(topo_name)
    replicate_dirs = [
        path
        for path in glob.glob(topo_result + "/*/")
        if os.path.basename(path.rstrip("/")).isnumeric()
    ]
    # List of dataframe to concatenate at the end
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
        # Getting the gk columns
        gk_cols = [col for col in sol_df.columns if col.startswith("gk_")]
        gene_cols = [col.replace("gk_", "") for col in gk_cols]
        # Keep only rows where 'State' does NOT contain the substring 'nan'
        sol_df = sol_df[~sol_df["State"].str.contains("nan", na=False)]
        # Keep only the rows which reach steady state
        sol_df = sol_df[sol_df["SteadyStateFlag"] == 1]
        # Resetting the index to avoid merginf problem witht eh new redorder state column
        sol_df = sol_df.reset_index(drop=True)
        # Creating a state df
        nstate_df = sol_df["State"].str.strip("'").apply(list)
        nstate_df = pd.DataFrame(nstate_df.to_list(), columns=gene_cols)
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
        nstate_df = nstate_df[new_order].astype(str).apply("".join, axis=1)
        sol_df["State"] = nstate_df
        print(sol_df)
        #     param_df = pd.read_parquet(
        #         os.path.join(rep_dir, f"{topo_name}_params_{rep_num}.parquet")
        #     )
        #     # sol_df = gk_normalise_solutions(sol_df, param_df, keep_gk=True)
        # Adding Rep number
        sol_df["RepNum"] = int(rep_num)
        soldf.append(sol_df)
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
        steady_type_df.append(multi_stable_df)
    soldf = pd.concat(soldf, axis=0)
    # gk_cols = [cl for cl in soldf.columns if "gk_" in cl]
    # print(soldf[gk_cols])
    # sns.histplot(
    #     sol_df[gk_cols].values.flatten(),
    #     bins=100,
    # )
    # plt.show()
    # plt.close()
    steady_counts_df = pd.concat(steady_counts_df, axis=0)
    steady_counts_df["Motif"] = topo_name
    # C print(steady_counts_df)
    steady_counts_df.to_csv(
        os.path.join(topo_result, f"{topo_name}_StateCounts.csv"), index=False
    )
    steady_type_df = pd.concat(steady_type_df, axis=0)
    steady_type_df["Motif"] = topo_name
    # print(steady_type_df)
    steady_type_df.to_csv(
        os.path.join(topo_result, f"{topo_name}_StateTypes.csv"), index=False
    )
    # Create the barplot
    plt.figure(figsize=(5, 6))
    steady_counts_df["State"] = steady_counts_df["State"].str.strip("'")
    steady_counts_df = steady_counts_df[steady_counts_df["Fraction"] >= 0.05]
    # Sort by Fraction
    steady_counts_df = steady_counts_df.sort_values("Fraction", ascending=False)
    ax = sns.barplot(
        data=steady_counts_df,
        x="State",
        y="Fraction",
        estimator="mean",  # Use mean across replicates
        errorbar="sd",  # For Seaborn >=0.12; use ci="sd" for older versions
        capsize=0.2,
        edgecolor="black",
        color=sns.color_palette("deep")[4],
    )
    plt.xticks(rotation=90)
    # Annotate each bar individually
    for bar in ax.patches:
        height = bar.get_height()
        x = bar.get_x() + bar.get_width() / 2
        label = f"{height:.2f}"

        # Adjust label position slightly above the bar (e.g., 5% of height)
        offset = height * 0.23
        ax.annotate(
            label,
            (x, height + offset),
            ha="center",
            va="bottom",
            fontsize=14,
        )
    ax.set_ylim(0, ax.get_ylim()[1] * 1.1)  # Increase top y-limit by 10%
    # plt.title(f"Steady State Distribution of {topo_name}")
    plt.ylabel("Fraction of States")
    plt.xlabel("State")
    plt.tight_layout()
    plt.savefig(os.path.join(topo_result, f"{topo_name}_StateCounts.png"), dpi=300)
    plt.savefig(os.path.join(topo_result, f"{topo_name}_StateCounts.svg"), dpi=300)
    plt.close()
    # plt.show()
    plt.figure(figsize=(5.5, 5))
    ax = sns.barplot(
        data=steady_type_df,
        x="Stability Type",
        y="Fraction",
        estimator="mean",  # Use mean across replicates
        errorbar="sd",  # For Seaborn >=0.12; use ci="sd" for older versions
        capsize=0.2,
        edgecolor="black",
        color=sns.color_palette("deep")[2],
    )
    for container in ax.containers:
        ax.bar_label(
            container,
            fmt="%.2f",  # <-- format to 2 decimal points
            padding=52,  # <-- adjust distance from top of bar
            fontsize=14,  # <-- optional: control font size
        )
    ax.set_ylim(0, ax.get_ylim()[1] * 1.1)  # Increase top y-limit by 10%
    # plt.title(f"Multi-Stability Type Distribution of {topo_name}")
    plt.ylabel("Fraction of States")
    plt.xlabel("Number of Stable States")
    plt.tight_layout()
    # plt.savefig(os.path.join(topo_result, f"{topo_name}_StateType.svg"))
    plt.savefig(os.path.join(topo_result, f"{topo_name}_StateType.png"), dpi=300)
    plt.savefig(os.path.join(topo_result, f"{topo_name}_StateType.svg"), dpi=300)
    # plt.show()
    plt.close()
    # Running PCA
    run_and_plot_pca(soldf, new_order, n_components=10)

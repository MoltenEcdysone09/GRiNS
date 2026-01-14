import pandas as pd
import numpy as np
from glob import glob
from scipy.stats import qmc, truncnorm, loguniform, lognorm
import itertools as it
import warnings
from grins.reg_funcs import psH, nsH
from typing import Tuple, List, Union

# Suppress the specific warning
warnings.filterwarnings(
    "ignore",
    message="The balance properties of Sobol' points require n to be a power of 2",
    category=UserWarning,
)


def parse_topos(topofile: str, save_cleaned: bool = False) -> pd.DataFrame:
    """
    Parse and cleans the given topofile and return the dataframe. It is expeced that the topo files is a tab-separated file with three columns: Source, Target, and Type. White spaces can also be used to separate the columns.

    For nodes that are not alphanumeric, the function will replace the non-alphanumeric characters with an underscore and prepend "Node_" if the node name does not start with an alphabet. The cleaned topology file will be saved if the save_cleaned flag is set to True.

    Parameters
    ----------
    topofile : str
        The path to the topofile.
    save_cleaned : bool, optional
        If True, save the cleaned topology file. Defaults to False.

    Returns
    -------
    topo_df : pd.DataFrame
        The parsed dataframe.
    """
    topo_df = pd.read_csv(topofile, sep=r"\s+")
    if topo_df.shape[1] != 3:
        raise ValueError(
            "The topology file should have three columns: Source, Target, and Type."
        )
    if not all(col in topo_df.columns for col in ["Source", "Target", "Type"]):
        raise ValueError(
            "The topology file should have the columns: Source, Target, and Type."
        )
    topo_df = topo_df[["Source", "Target", "Type"]]
    # Clean up node names: replace non-alphanumerics and prepend "Node_" if needed.
    topo_df["Source"] = (
        topo_df["Source"]
        .str.replace(r"\W", "_", regex=True)
        .apply(lambda x: f"Node_{x}" if not x[0].isalpha() else x)
    )
    topo_df["Target"] = (
        topo_df["Target"]
        .str.replace(r"\W", "_", regex=True)
        .apply(lambda x: f"Node_{x}" if not x[0].isalpha() else x)
    )
    if topo_df["Type"].nunique() > 2:
        raise ValueError(f"Check the topo file: {topofile}")
    if save_cleaned:
        topo_df.to_csv(
            topofile.replace(".topo", "_cleaned.topo"), sep="\t", index=False
        )
    return topo_df


def _get_regtype(sn: str, tn: str, topo_df: pd.DataFrame) -> str:
    """
    Get the type of regulation for a given source and target node.

    Parameters
    ----------
    sn : str
        The source node.
    tn : str
        The target node.
    topo_df : pd.DataFrame
        The DataFrame containing the topology information.

    Returns
    -------
    str
        The type of regulation, either "ActFld_{sn}_{tn}" for activation or "InhFld_{sn}_{tn}" for inhibition.
    """
    reg_type = topo_df[(topo_df["Source"] == sn) & (topo_df["Target"] == tn)][
        "Type"
    ].iloc[0]
    return f"ActFld_{sn}_{tn}" if reg_type == 1 else f"InhFld_{sn}_{tn}"


def gen_param_names(topo_df: pd.DataFrame) -> Tuple[List[str], List[str], List[str]]:
    """
    Generate parameter names based on the given topology dataframe.

    Parameters
    ----------
    topo_df : pd.DataFrame
        The topology dataframe containing the information about the nodes and edges.

    Returns
    -------
    tuple
        A tuple containing the parameter names, unique target node names, and unique source node names.
    """
    target_nodes = list(topo_df["Target"].unique())
    source_nodes = list(topo_df["Source"].unique())
    unique_nodes = sorted(set(source_nodes + target_nodes))
    param_names = [f"Prod_{n}" for n in unique_nodes] + [
        f"Deg_{n}" for n in unique_nodes
    ]
    for tn in unique_nodes:
        sources = topo_df[topo_df["Target"] == tn]["Source"].sort_values()
        for sn, p in it.product(sources, ["Fld", "Thr", "Hill"]):
            param_names.append(
                f"{p}_{sn}_{tn}" if p != "Fld" else _get_regtype(sn, tn, topo_df)
            )
    return param_names, target_nodes, source_nodes


def scale_array(
    arr: np.ndarray, min_val: float, max_val: float, round_int: bool = False
) -> np.ndarray:
    """
    Scale the array values to the range [min_val, max_val]. If round_int is True, np.ceil the results
    and np.clip to the range.

    Parameters
    ----------
    arr : np.ndarray
        The input array to be scaled.
    min_val : float
        The minimum value of the scaled array.
    max_val : float
        The maximum value of the scaled array.
    round_int : bool, optional
        If True, ceil the results and clip to the range. Defaults to False.

    Returns
    -------
    np.ndarray
        The scaled array with values in the range [min_val, max_val].
    """
    scaled = min_val + (max_val - min_val) * (arr - arr.min()) / (arr.max() - arr.min())
    if round_int:
        scaled = np.round(scaled)
        scaled = np.clip(scaled, min_val, max_val)
    return scaled


def _get_rng(rng: int | np.random.Generator | None) -> np.random.Generator:
    """
    Ensure the input is a valid NumPy Generator.

    This helper standardizes the way random states are handled, supporting
    seed integers, existing Generator instances, or None.

    Parameters
    ----------
    rng : int, np.random.Generator, or None
        The seed or generator to use.
        - If `int`, it's used as a seed to create a new `Generator`.
        - If `np.random.Generator`, it is returned as-is.
        - If `None`, a new `Generator` is created using entropy from the OS.

    Returns
    -------
    np.random.Generator
        A validated or newly initialized NumPy random number generator.
    """
    if isinstance(rng, np.random.Generator):
        return rng
    return np.random.default_rng(rng)


def _gen_sobol_seq(
    dimensions: int,
    num_points: int,
    optimise: bool = False,
    rng: int | np.random.Generator | None = None,
) -> np.ndarray:
    """
    Generate a Sobol sequence.

    Parameters:
        dimensions (int): The number of dimensions for the Sobol sequence.
        num_points (int): The number of points to generate in the sequence.
        optimise (bool, optional): If True, apply Lloyd's optimization to the sequence. Default is False.
        rng (int, Generator, optional): Seed or RNG for scrambling. Default is None.

    Returns:
        np.ndarray: An array containing the generated Sobol sequence.
    """
    sampler = qmc.Sobol(
        d=dimensions,
        scramble=True,
        seed=rng,
        optimization="lloyd" if optimise else None,
    )
    return sampler.random(num_points)


def _gen_uniform_seq(
    dimension: int, num_points: int, rng: int | np.random.Generator | None = None
) -> np.ndarray:
    """
    Generate a sequence of uniformly distributed random points.

    Parameters:
        dimension (int): The number of dimensions for each point.
        num_points (int): The number of points to generate.
        rng (int, Generator, optional): Seed or RNG for scrambling. Default is None.

    Returns:
        np.ndarray: A 2D array of shape (num_points, dimension) containing the generated points.
    """
    # rng = _get_rng(rng)
    return rng.uniform(low=0, high=1, size=(num_points, dimension))


def _gen_loguniform_seq(
    dimension: int, num_points: int, rng: int | np.random.Generator | None = None
) -> np.ndarray:
    """
    Generate a sequence of points sampled from a log-uniform distribution.

    Parameters:
        dimension (int): The dimensionality of the points.
        num_points (int): The number of points to generate.
        rng (int, Generator, optional): Seed or RNG for scrambling. Default is None.

    Returns:
        np.ndarray: An array of shape (num_points, dimension) containing the generated points.
    """
    # rng = _get_rng(rng)
    return _gen_uniform_seq(dimension=dimension, num_points=num_points, rng=rng)


def _gen_latin_hypercube(
    dimension: int,
    num_points: int,
    optimise: bool = False,
    rng: int | np.random.Generator | None = None,
) -> np.ndarray:
    """
    Generate a Latin Hypercube sample.

    Parameters:
        dimension (int): The number of dimensions for the sample.
        num_points (int): The number of points to generate.
        optimise (bool, optional): Whether to use Lloyd's algorithm for optimization. Defaults to False.
        rng (int, Generator, optional): Seed or RNG for scrambling. Default is None.

    Returns:
        np.ndarray: A numpy array containing the generated Latin Hypercube sample.
    """
    sampler = qmc.LatinHypercube(
        d=dimension, scramble=True, seed=rng, optimization="lloyd" if optimise else None
    )
    return sampler.random(num_points)


def _gen_normal(
    dimension: int,
    num_points: int,
    std_dev: float = 1,
    rng: int | np.random.Generator | None = None,
) -> np.ndarray:
    """
    Generate a set of points from a truncated normal distribution.

    This function generates `num_points` points in `dimension`-dimensional space,
    where each point is drawn from a truncated normal distribution with a mean of 0
    and a standard deviation of `std_dev`. The distribution is truncated to avoid
    extreme values, limiting the range to ±3 standard deviations.

    Parameters:
        dimension (int): The number of dimensions for each point.
        num_points (int): The number of points to generate.
        std_dev (float, optional): The standard deviation of the normal distribution. Default is 1.
        rng (int, Generator, optional): Seed or RNG for reproducibility. Default is None.

    Returns:
        np.ndarray: A NumPy array of shape (num_points, dimension) containing the generated points.
    """
    # rng = _get_rng(rng)
    lower, upper = -3, 3
    # Use truncated normal to avoid extreme values (±3 std deviations)
    return truncnorm.rvs(
        a=lower,
        b=upper,
        loc=0,
        scale=std_dev,
        size=(num_points, dimension),
        random_state=rng,
    )


def _gen_lognormal(
    dimension: int,
    num_points: int,
    std_dev: float = 1,
    rng: int | np.random.Generator | None = None,
) -> np.ndarray:
    """
    Generate a log-normal distribution.

    Parameters:
        dimension (int): The number of dimensions for each point.
        num_points (int): The number of points to generate.
        std_dev (float, optional): The standard deviation of the underlying normal distribution. Default is 1.
        rng (int, Generator, optional): Seed or RNG for reproducibility. Default is None.

    Returns:
        np.ndarray: An array of shape (num_points, dimension) containing the generated log-normal points.
    """
    # The log scalign will be handled afterwards
    # return np.exp(_gen_normal(dimension, num_points, std_dev, rng=rng))
    return _gen_uniform_seq(dimension=dimension, num_points=num_points, rng=rng)


def sample_distribution(
    method: str,
    dimension: int,
    num_points: int,
    std_dev: float | None = None,
    optimise: bool = False,
    rng: int | np.random.Generator | None = None,
) -> np.ndarray:
    """
    Generates a sample distribution based on the specified method.

    Parameters
    ----------
    method : str
        The sampling method to use. Options are "Sobol", "LHS", "Uniform", "LogUniform", "Normal", "LogNormal".
    dimension : int
        The number of dimensions for the sample points.
    num_points : int
        The number of sample points to generate.
    std_dev : Optional[float], optional
        The standard deviation for the "Normal" and "LogNormal" distributions. Defaults to 1.0 if not provided.
    optimise : bool, optional
        Whether to optimise the sampling process. Applicable for "Sobol" and "LHS" methods.
    rng : int, np.random.Generator, or None, optional
        The seed or NumPy Generator to use for scrambling and shuffling.
        If None, the global random state is not used; a new Generator
        is initialized.

    Returns
    -------
    np.ndarray
        An array of sample points generated according to the specified method.

    Raises
    ------
    ValueError
        If an unknown sampling method is specified.
    """
    # Initialise the random number generator
    rng_obj = _get_rng(rng)
    if method == "Sobol":
        dist_arr = _gen_sobol_seq(dimension, num_points, optimise, rng_obj)
        # return _gen_sobol_seq(dimension, num_points, optimise)
    elif method == "LHS":
        dist_arr = _gen_latin_hypercube(dimension, num_points, optimise, rng_obj)
        # return _gen_latin_hypercube(dimension, num_points, optimise)
    elif method == "Uniform":
        dist_arr = _gen_uniform_seq(dimension, num_points, rng_obj)
        # return _gen_uniform_seq(dimension, num_points)
    elif method == "LogUniform":
        dist_arr = _gen_loguniform_seq(dimension, num_points, rng_obj)
        # return _gen_loguniform_seq(dimension, num_points)
    elif method == "Normal":
        dist_arr = _gen_normal(
            dimension, num_points, std_dev if std_dev is not None else 1.0, rng_obj
        )
        # return _gen_normal(
        #     dimension, num_points, std_dev if std_dev is not None else 1.0
        # )
    elif method == "LogNormal":
        dist_arr = _gen_lognormal(
            dimension, num_points, std_dev if std_dev is not None else 1.0, rng_obj
        )
        # return _gen_lognormal(
        #     dimension, num_points, std_dev if std_dev is not None else 1.0
        # )
    else:
        raise ValueError(f"Unknown sampling method: {method}")
    # Shuffle the array to avoid the correlation between the parameters
    # Use the generator for the final shuffle to maintain isolation
    flat_arr = dist_arr.flatten()
    rng_obj.shuffle(flat_arr)
    return flat_arr.reshape(dist_arr.shape)
    # return np.random.permutation(dist_arr.flatten()).reshape(dist_arr.shape)


def sample_param_df(
    prange_df: pd.DataFrame,
    num_params: int = 2**10,
    rng: int | np.random.Generator | None = None,
) -> pd.DataFrame:
    """
    Samples parameter values based on the provided parameter range DataFrame.

    This function takes a DataFrame containing parameter ranges and sampling
    methods, and generates a new DataFrame with sampled parameter values. The
    sampling is performed according to the specified methods and standard
    deviations (if provided).

    Parameters
    ----------
    prange_df : pd.DataFrame
        A DataFrame containing parameter ranges and sampling methods. It must
        include at least the columns "Parameter" and "Sampling". Optionally, it
        can include "StdDev", "Minimum", and "Maximum" columns.
    num_params : int, optional
        The number of parameter samples to generate for each parameter. Default
        is 1024 (2**10).
    rng: int, Generator, optional
        Seed or RNG for reproducibility. Default is None.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the sampled parameter values. The columns are
        ordered according to the original order of parameters in `prange_df`.
    """
    # Save the original order of parameters
    original_order = prange_df["Parameter"].tolist()
    # Sort the DataFrame by "Sampling" and "StdDev" columns
    sort_cols = ["Sampling"] + (["StdDev"] if "StdDev" in prange_df.columns else [])
    prange_df = prange_df.sort_values(by=sort_cols)
    # Get the parameter names in the order they appear in the sampling method sorted DataFrame
    param_order = prange_df["Parameter"].tolist()
    # Dictionary to store sampled values for each parameter
    sampled_values = {}

    # Group by sampling settings (including StdDev if available)
    group_cols = ["Sampling"] + (["StdDev"] if "StdDev" in prange_df.columns else [])
    for _, group in prange_df.groupby(group_cols, sort=False, dropna=False):
        method = group["Sampling"].iloc[0]
        std_val = group["StdDev"].iloc[0] if "StdDev" in group.columns else None
        dims = group.shape[0]
        samples = sample_distribution(
            method, dims, num_params, std_dev=std_val, rng=rng
        )
        # Assign each sampled column to its corresponding parameter name.
        for i, param in enumerate(group["Parameter"]):
            sampled_values[param] = samples[:, i]
    # Create a DataFrame from the sampled values
    sampled_df = pd.DataFrame({p: sampled_values[p] for p in param_order})

    # Scale the sampled values to the specified range
    def scale_col(col: pd.Series, param_name: str) -> pd.Series:
        row = prange_df[prange_df["Parameter"] == param_name].iloc[0]
        min_val, max_val = row["Minimum"], row["Maximum"]
        std_val = row.get("StdDev", 1.0)
        # Determine the target range for sampling.
        # For InhFld sample in reciprocal space [1/max, 1/min]
        if "InhFld" in param_name:
            target_min, target_max = 1 / max_val, 1 / min_val
        else:
            target_min, target_max = min_val, max_val
        # Handling LogUniform scaling
        if row["Sampling"] == "LogUniform":
            sample_vals = pd.Series(loguniform.ppf(col.values, target_min, target_max))
        # Handling LogNormal scaling
        elif row["Sampling"] == "LogNormal":
            # For LogNormal, mu is the mean of log-space
            log_m, log_M = np.log(target_min), np.log(target_max)
            mu = (log_m + log_M) / 2
            sample_vals = pd.Series(
                lognorm.ppf(col.values, s=std_val, scale=np.exp(mu))
            )
        # Handling all other non-log scaling
        else:
            sample_vals = pd.Series(
                scale_array(
                    col.values,
                    target_min,
                    target_max,
                )
            )
        # Rounding the hill values
        if "Hill" in param_name:
            sample_vals = np.round(sample_vals)
            # Ensure it stays within [min+1, max] like scale_array does
            sample_vals = np.clip(sample_vals, target_min, target_max)
        if "InhFld" in param_name:
            # Invert the scaled values and return the series
            sample_vals = 1 / sample_vals
        return sample_vals

    # Scale the sampled values for each parameter
    for col in sampled_df.columns:
        sampled_df[col] = scale_col(sampled_df[col], col)
    # Return the sampled DataFrame with the original parameter order
    return sampled_df[original_order]


# Functionto calculate the median g/k values based on the incoming edges. This function will be used by the get_thr_ranges function to get the proper threshold ranges which agrees with the half functional rule.
def _get_updated_gkn_hills(
    gk_n: np.ndarray,
    in_edge_params: pd.DataFrame,
    in_edge_topo: pd.DataFrame,
    num_params: int = 2**10,
    rng: int | np.random.Generator | None = None,
) -> float:
    """
    Update the g/k values for a node using incoming edges and their Hills equations. Used by the get_thr_ranges function to generate threshold ranges which follow half functional rule.

    Parameters:
        gk_n (np.ndarray): Array of g/k values for the node.
        in_edge_params (pd.DataFrame): DataFrame containing parameters for incoming edges.
        in_edge_topo (pd.DataFrame): DataFrame containing the topology of incoming edges.
        num_params (int, optional): Number of parameter samples to generate. Default is 1024.
        rng (int, Generator, optional): Seed or RNG for reproducibility. Default is None.

    Returns:
        float: The median of the product of updated g/k values across all incoming edges.
    """
    # print(f"Incoming edge parameters:\n{in_edge_params}")
    # Create a DataFrame to store the updated g/k values
    gk_hills = pd.DataFrame({"src_gk_n": gk_n})
    # Sample parameters for incoming edges
    inedg_param_samples = sample_param_df(in_edge_params, num_params, rng=rng)
    # Calculate the updated g/k values based on the Hills equation for each incoming edge
    for _, row in in_edge_topo.iterrows():
        # Get the source and target nodes, and the regulation type
        insrc_nd, tgt_nd, type_reg = row["Source"], row["Target"], row["Type"]
        # print(
        #     f"Processing incoming edge from {insrc_nd} to {tgt_nd} of type {type_reg}"
        # )
        # Calculate the updated g/k value based on the regulation type
        if type_reg == 1:
            fld = inedg_param_samples[f"ActFld_{insrc_nd}_{tgt_nd}"].values
            g = inedg_param_samples[f"Prod_{insrc_nd}"].values
            k = inedg_param_samples[f"Deg_{insrc_nd}"].values
            n = inedg_param_samples[f"Hill_{insrc_nd}_{tgt_nd}"].values
            thr = inedg_param_samples[f"Thr_{insrc_nd}_{tgt_nd}"].values
            gk_hills[f"AcH_{insrc_nd}_{tgt_nd}"] = psH(g / k, fld, thr, n)
        else:
            fld = inedg_param_samples[f"InhFld_{insrc_nd}_{tgt_nd}"].values
            g = inedg_param_samples[f"Prod_{insrc_nd}"].values
            k = inedg_param_samples[f"Deg_{insrc_nd}"].values
            n = inedg_param_samples[f"Hill_{insrc_nd}_{tgt_nd}"].values
            thr = inedg_param_samples[f"Thr_{insrc_nd}_{tgt_nd}"].values
            gk_hills[f"InH_{insrc_nd}_{tgt_nd}"] = nsH(g / k, fld, thr, n)
    # Calculate the median of the product of updated g/k values across all incoming edges and return the median g/k value
    return gk_hills.prod(axis=1).median()


# Function to get the threshold ranges for a given source node. This function will be used by the add_thr_rows function to add threshold-related rows to the parameter range DataFrame.
def get_thr_ranges(
    source_node: str,
    topo_df: pd.DataFrame,
    prange_df: pd.DataFrame,
    num_params: int = 2**10,
    # optimise: bool = False,
    rng: int | np.random.Generator | None = None,
) -> float:
    """
    Calculate the median threshold range for a given source node.

    Parameters
    ----------
    source_node : str
        The source node for which the threshold range is calculated.
    topo_df : pd.DataFrame
        DataFrame containing the topology information.
    prange_df : pd.DataFrame
        DataFrame containing the parameter ranges.
    num_params : int, optional
        Number of parameters to sample. Defaults to 1024.
    rng: int, Generator, optional
        Seed or RNG for reproducibility. Default is None.

    Returns
    -------
    float
        The median threshold range for the given source node.
    """
    # Get the production and degradation rates for the source node
    sn_params = prange_df[
        prange_df["Parameter"].str.contains(f"Prod_{source_node}|Deg_{source_node}")
    ]
    # Sample the production and degradation rates for the source node
    sn_gk = sample_param_df(sn_params, num_params, rng=rng)
    # Calculate the g/k values for the source node
    sn_gk_n = (sn_gk[f"Prod_{source_node}"] / sn_gk[f"Deg_{source_node}"]).to_numpy()
    # Get the incoming edges for the source node
    in_edge_topo = topo_df[topo_df["Target"] == source_node]
    # If there are incoming edges, calculate the updated g/k values based on the Hills equation
    if not in_edge_topo.empty:
        sources = in_edge_topo["Source"].unique()
        isn = "|".join(sources)
        # Get the parameters for the incoming edges
        # Regex pattern to match the parameters
        pattern = (
            rf"^((ActFld|InhFld|Thr|Hill)_({isn})_{source_node}|(Prod|Deg)_({isn}))$"
        )
        in_edge_params = prange_df[prange_df["Parameter"].str.contains(pattern)].copy()
        # in_edge_params = prange_df[
        #     prange_df["Parameter"].str.contains(
        #         f"Fld_{isn}_{source_node}|Thr_{isn}_{source_node}|Hill_{isn}_{source_node}"
        #     )
        #     | prange_df["Parameter"].str.contains(f"Prod_{isn}|Deg_{isn}")
        # ]
        # Sample the parameters for the incoming edges
        isn_gk = sample_param_df(
            in_edge_params[in_edge_params["Parameter"].str.contains("Prod_|Deg_")],
            num_params,
            rng=rng,
        )
        # Calculate the updated g/k values based on the Hills equation for the incoming edges
        for in_node in in_edge_topo["Source"].unique():
            # print(f"Processing incoming edge from {in_node} to {source_node}")
            in_gk = isn_gk[f"Prod_{in_node}"] / isn_gk[f"Deg_{in_node}"]
            in_gk_median = np.median(in_gk)
            # print(f"Median g/k value for {in_node}: {in_gk_median}")
            in_edge_params.loc[
                in_edge_params["Parameter"].str.contains(
                    f"Thr_{in_node}_{source_node}"
                ),
                ["Minimum", "Maximum"],
            ] = [0.02 * in_gk_median, 1.98 * in_gk_median]
        # Update the g/k values for the source node based on the incoming edges and return the median g/k value
        return _get_updated_gkn_hills(
            sn_gk_n, in_edge_params, in_edge_topo, num_params, rng=rng
        )
    else:
        # If there are no incoming edges, return the median g/k value for the source node
        return np.median(sn_gk[f"Prod_{source_node}"] / sn_gk[f"Deg_{source_node}"])


# Function to add threshold-related rows to the parameter range DataFrame. This function will be used by the get_param_range_df function to generate a parameter range DataFrame from the topology DataFrame.
def add_thr_rows(
    prange_df: pd.DataFrame,
    topo_df: pd.DataFrame,
    num_params: int = 2**10,
    rng: int | np.random.Generator | None = None,
) -> pd.DataFrame:
    """
    This function modifies the given parameter range DataFrame (`prange_df`) by adding
    threshold-related rows based on the topology DataFrame (`topo_df`). It calculates
    the median threshold values for source nodes and adjusts the minimum and maximum
    values for the corresponding parameters. Additionally, it scales the production
    parameters for the source nodes.

    Parameters
    ----------
    prange_df : pd.DataFrame
        The DataFrame containing parameter ranges.
    topo_df : pd.DataFrame
        The DataFrame containing the network topology.
    num_params : int, optional
        The number of parameters to consider. Default is 1024.
    rng: int, Generator, optional
        Seed or RNG for reproducibility. Default is None.

    Returns
    -------
    pd.DataFrame
        The modified parameter range DataFrame with added threshold-related rows.
    """
    # Get the parameter names, target nodes, and source nodes
    param_names, target_nodes, source_nodes = gen_param_names(topo_df)
    # Iterate over the unique source and target nodes
    for sn in set(target_nodes + source_nodes):
        # print(f"Processing node: {sn}")
        # If the node is a source node, its aplification factor is calculated if the median threshold value is less than 0.01
        if sn in source_nodes:
            median_thr_val = get_thr_ranges(sn, topo_df, prange_df, num_params, rng=rng)
            # print(f"Median threshold value for {sn}: {median_thr_val}")
            amplify_val = 1.0
            # Calculate the amplification factor if the median threshold value is less than 0.01
            if (median_thr_val * 0.02) < 0.010:
                exp_val = np.floor(np.log10(np.abs(median_thr_val * 0.02))).astype(int)
                amplify_val = 10 ** (-exp_val - 2)
            # Update the production and threshold values for the source node
            prange_df.loc[
                prange_df["Parameter"].str.contains(f"Thr_{sn}"), "Minimum"
            ] = median_thr_val * 0.02 * amplify_val
            prange_df.loc[
                prange_df["Parameter"].str.contains(f"Thr_{sn}"), "Maximum"
            ] = median_thr_val * 1.98 * amplify_val
            prange_df.loc[
                prange_df["Parameter"] == f"Prod_{sn}", ["Minimum", "Maximum"]
            ] *= amplify_val
    return prange_df


# Function to generate the parameter range dataframe from the topology of the network.
def gen_param_range_df(
    topo_df: pd.DataFrame,
    num_params: int = 2**10,
    sampling_method: Union[str, dict] = "Sobol",
    thr_rows: bool = True,
    rng: int | np.random.Generator | None = None,
) -> pd.DataFrame:
    """
    Generate a parameter range DataFrame from the topology DataFrame.

    Parameters
    ----------
    topo_df : pd.DataFrame
        The topology DataFrame containing the network structure.
    num_params : int, optional
        The number of parameters to generate. Default is 1024.
    sampling_method : Union[str, dict], optional
        The sampling method to use. Can be a string specifying a single method for all parameters or a dictionary specifying methods for individual parameters. Default is "Sobol".
    thr_rows : bool, optional
        Whether to add threshold rows to the DataFrame. Default is True.
    rng: int, Generator, optional
        Seed or RNG for reproducibility. Default is None.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the parameter ranges and sampling methods.
    """
    # Generate the parameter names, target nodes, and source nodes
    param_names, target_nodes, source_nodes = gen_param_names(topo_df)
    # Create a DataFrame to store the parameter ranges
    prange_df = pd.DataFrame({"Parameter": param_names})
    # Set the default minimum and maximum values for the parameters
    prange_df.loc[
        prange_df["Parameter"].str.contains("Prod_"), ["Minimum", "Maximum"]
    ] = [1.0, 100.0]
    prange_df.loc[
        prange_df["Parameter"].str.contains("Deg_"), ["Minimum", "Maximum"]
    ] = [0.1, 1.0]
    prange_df.loc[
        prange_df["Parameter"].str.contains("ActFld_"), ["Minimum", "Maximum"]
    ] = [1.0, 100.0]
    prange_df.loc[
        prange_df["Parameter"].str.contains("InhFld_"), ["Minimum", "Maximum"]
    ] = [0.01, 1.0]
    prange_df.loc[
        prange_df["Parameter"].str.contains("Hill"), ["Minimum", "Maximum"]
    ] = [1.0, 6.0]
    # Set the sampling method for each parameter, if the sampling method is a dictionary set the specific method for the specified parameters
    if isinstance(sampling_method, str):
        prange_df["Sampling"] = sampling_method
        if sampling_method in ["Normal", "LogNormal"]:
            prange_df["StdDev"] = 1.0
    else:
        for param, method in sampling_method.items():
            prange_df.loc[prange_df["Parameter"].str.contains(param), "Sampling"] = (
                method
            )
        # If the sampling method is not specified for a parameter, set it to "Sobol"
        prange_df["Sampling"] = prange_df["Sampling"].fillna("Sobol")
        if "StdDev" in prange_df.columns:
            prange_df["StdDev"] = 1.0
    # Fill the threshold rows of the parameter range DataFrame
    if thr_rows:
        prange_df = add_thr_rows(prange_df, topo_df, num_params, rng=rng)
    return prange_df


def gen_param_df(
    prange_df: pd.DataFrame = None,
    num_params: int = 2**10,
    topo_df: pd.DataFrame = None,
    sampling_method: Union[str, dict] = "Sobol",
    thr_rows: bool = True,
    rng: int | np.random.Generator | None = None,
) -> pd.DataFrame:
    """
    Generate the final parameter DataFrame by sampling parameters.
    Parameters are grouped by their 'Sampling' (and 'StdDev' if present) to ensure
    that parameters in the same group follow the same distribution in the higher dimensions.
    The final DataFrame columns are arranged in the same order as in prange_df.
    The sampling methods can be: 'Sobol', 'LHS', 'Uniform', 'LogUniform', 'Normal', 'LogNormal'.

    Parameters
    ----------
    prange_df : pd.DataFrame, optional
        DataFrame with columns ["Parameter", "Minimum", "Maximum", "Sampling", ...].
    num_params : int, optional
        Number of samples to generate per parameter. Default is 1024.
    topo_df : pd.DataFrame, optional
        DataFrame containing the network topology information.
    sampling_method : Union[str, dict], optional
        The sampling method to use. Can be a string specifying a single method for all parameters or a dictionary specifying methods for individual parameters. Default is "Sobol". The methods can be: 'Sobol', 'LHS', 'Uniform', 'LogUniform', 'Normal', 'LogNormal'.
    thr_rows : bool, optional
        Whether to add threshold rows to the DataFrame. Default is True.
    rng: int, Generator, optional
        Seed or RNG for reproducibility. Default is None.

    Returns
    -------
    pd.DataFrame
        DataFrame of sampled and scaled parameters.
    """
    # If the parameter range dataframe is not given
    if prange_df is None:
        # Generating the parameter range dataframe from the topology of the network
        prange_df = gen_param_range_df(
            topo_df=topo_df,
            num_params=num_params,
            sampling_method=sampling_method,
            thr_rows=thr_rows,
            rng=rng,
        )
    # Use the sample_param_df function to sample the parameters
    param_df = sample_param_df(prange_df, num_params, rng=rng)
    # Add the ParamNum column to the DataFrame
    param_df = param_df.assign(ParamNum=param_df.index + 1)
    return param_df


def gen_init_cond(
    topo_df: pd.DataFrame,
    num_init_conds: int = 2**10,
    rng: int | np.random.Generator | None = None,
) -> pd.DataFrame:
    """
    Generate initial conditions for each node based on the topology.

    Parameters
    ----------
    topo_df : pd.DataFrame
        DataFrame containing the topology information.
    num_init_conds : int, optional
        Number of initial conditions to generate. Default is 2**10.
    rng: int, Generator, optional
        Seed or RNG for reproducibility. Default is None.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the generated initial conditions for each node.
    """
    _, target_nodes, source_nodes = gen_param_names(topo_df)
    unique_nodes = sorted(set(source_nodes + target_nodes))
    init_conds = _gen_sobol_seq(len(unique_nodes), num_init_conds, rng=rng)
    # Scale initial conditions between 1 and 100
    init_conds = 1 + init_conds * (100 - 1)
    initcond_df = pd.DataFrame(init_conds, columns=unique_nodes)
    # A new columns for the intial condition numbers
    # initcond_df["InitCondNum"] = initcond_df.index + 1
    initcond_df = initcond_df.assign(InitCondNum=initcond_df.index + 1)
    return initcond_df

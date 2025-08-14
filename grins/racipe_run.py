from glob import glob
import sys
import os
from grins.gen_diffrax_ode import gen_diffrax_odesys
from grins.gen_params import (
    gen_param_df,
    gen_init_cond,
    gen_param_range_df,
    parse_topos,
)
from importlib import import_module
from diffrax import (
    diffeqsolve,
    ODETerm,
    SaveAt,
    Tsit5,
    PIDController,
    Event,
    steady_state_event,
)
import subprocess
from multiprocessing import Pool  # noqa: F401
import pandas as pd  # noqa: F401
import jax.numpy as jnp
from jax import jit, vmap
import numpy as np
import time
from typing import Union
from scipy.signal import find_peaks
import warnings


# Function to generate the required directory structure
# def gen_sim_dirstruct(topo_file, save_dir=".", num_replicates=3):
def gen_sim_dirstruct(
    topo_file: str, save_dir: str = ".", num_replicates: int = 3
) -> None:
    """
    Generate directory structure for simulation run.

    Parameters
    -----------
    topo_file : str
        Path to the topo file.
    save_dir : str, optional
        Directory to save the generated structure. Defaults to ".".
    num_replicates : int, optional
        Number of replicates to generate. Defaults to 3.

    Returns
    --------
    None
        Directory structure is created with the topo file name and three folders for the replicates.
    """
    # Get the topo file name
    topo_name = topo_file.split("/")[-1].split(".")[0]
    # Check if the folder with the name of topo file exists
    os.makedirs(f"{save_dir}/{topo_name}", exist_ok=True)
    # Move the topo file to the created folder
    subprocess.run(
        [
            "cp",
            topo_file,
            f"{save_dir}/{topo_name}/{topo_file.split('/')[-1]}",
        ]
    )
    # Make the replicate directories
    for rep in range(1, num_replicates + 1):
        os.makedirs(f"{save_dir}/{topo_name}/{rep:03}", exist_ok=True)
    return None


# Functiont to generate all the parameters related files with replicates
def gen_topo_param_files(
    topo_file: str,
    save_dir: str = ".",
    num_replicates: int = 3,
    num_params: int = 2**10,
    num_init_conds: int = 2**7,
    sampling_method: Union[str, dict] = "Sobol",
):
    """
    Generate parameter files for simulation.

    Parameters
    ----------
    topo_file : str
        The path to the topo file.
    save_dir : str, optional
        The directory where the parameter files will be saved. Defaults to ".".
    num_params : int, optional
        The number of parameter files to generate. Defaults to 2**10.
    num_init_conds : int, optional
        The number of initial condition files to generate. Defaults to 2**7.
    sampling_method : Union[str, dict], optional
        The method to use for sampling the parameter space. Defaults to 'Sobol'. For a finer control over the parameter generation look at the documentation of the gen_param_range_df function and gen_param_df function.

    Returns
    -------
    None
        The parameter files and initial conditions are generated and saved in the specified replicate directories.
    """
    # Get the name of the topo file
    topo_name = topo_file.split("/")[-1].split(".")[0]
    # Parse the topo file
    topo_df = parse_topos(topo_file)
    # # Generate the parameter names
    # param_names = gen_param_names(topo_df)
    # Get the unique nodes in the topo file
    # unique_nodes = sorted(set(param_names[1] + param_names[2]))
    # Generate the required directory structure
    gen_sim_dirstruct(topo_file, save_dir, num_replicates)
    # Specify directory where all the generated ode system file will be saved
    sim_dir = f"{save_dir}/{topo_file.split('/')[-1].split('.')[0]}"
    # Generate the ODE system for diffrax
    gen_diffrax_odesys(topo_df, topo_name, sim_dir)
    # Generate the parameter dataframe and save in each of the replicate folders
    for rep in range(1, num_replicates + 1):
        # Generate the parameter range dataframe
        param_range_df = gen_param_range_df(
            topo_df, num_params, sampling_method=sampling_method
        )
        # Save the parameter range dataframe
        param_range_df.to_csv(
            f"{sim_dir}/{rep:03}/{topo_name}_param_range_{rep:03}.csv",
            index=False,
            sep="\t",
        )
        # # Generate the parameter dataframe with the default values
        param_df = gen_param_df(param_range_df, num_params)
        # print(param_df)
        param_df.to_parquet(
            f"{sim_dir}/{rep:03}/{topo_name}_params_{rep:03}.parquet", index=False
        )
        # Generate the initial conditions dataframe
        initcond_df = gen_init_cond(topo_df=topo_df, num_init_conds=num_init_conds)
        # print(initcond_df)
        initcond_df.to_parquet(
            f"{sim_dir}/{rep:03}/{topo_name}_init_conds_{rep:03}.parquet",
            index=False,
        )
    print(f"Parameter and Intial Condition files generated for {topo_name}")
    return None


# Function to load the ODE system from a specified topology folder as a module so that the ODETerm can be initialized from the ODE system
def load_odeterm(topo_name, simdir):
    """
    Loads an ODE system from a specified topology module and returns an ODETerm object.

    Parameters
    ----------
    topo_name : str
        The name of the topology module to import.
    simdir : str
        The directory path where the topology module is located.

    Returns
    -------
    ODETerm
        An object representing the ODE system.

    Raises
    ------
    ImportError
        If the specified module cannot be imported.
    AttributeError
        If the module does not contain an attribute named 'odesys'.
    """
    sys.path.append(f"{simdir}")
    mod = import_module(f"{topo_name}")
    term = ODETerm(getattr(mod, "odesys"))
    return term


# Function to generate the combinations of initial conditions and parameters
def _gen_combinations(num_init_conds, num_params):
    """
    Generate combinations of initial conditions and parameters.

    This function generates all possible combinations of initial conditions and parameters
    using efficient operations from the JAX library. The combinations are generated as pairs of indices from the initial conditions and parameters dataframes.

    Parameters
    ----------
    num_init_conds : int
        The number of initial conditions.
    num_params : int
        The number of parameters.

    -------
    jnp.ndarray
        A 2D array where each row represents a combination of an initial condition and a parameter. The first column contains indices of initial conditions, and the second column contains indices of parameters.
    """
    # Generate the combinations more efficiently using numpy
    i = jnp.repeat(jnp.arange(num_init_conds), num_params)
    p = jnp.tile(jnp.arange(num_params), num_init_conds)
    icprm_comb = jnp.stack([i, p], axis=1)
    return icprm_comb


# Function to parameterise the ODE system and return the right functions
def parameterise_solveode(
    ode_term,
    solver,
    t0,
    t1,
    dt0,
    saveat,
    stepsize_controller,
    max_steps,
    initial_conditions,
    parameters,
):
    """
    Parameterise the ODE system and return the right functions.

    Parameters
    ----------
    term : ODETerm
        The ODE system to solve.
    solver : object
        The solver to be used for solving the ODE.
    t0 : float
        The initial time.
    t1 : float
        The final time.
    dt0 : float
        The initial time step.
    saveat : SaveAt
        The time steps to save the solution at.
    stepsize_controller : PIDController
        The controller for the step sizes.
    max_steps : int
        The maximum number of steps to take.
    initial_conditions : jnp.ndarray
        The initial conditions.
    Returns
    parameters : jnp.ndarray
        The parameters.

    Returns
    -------
    function
        The functions to solve the ODE system.
    """
    # Check if number of time steps to save is None
    if saveat.subs.ts is None:
        # Function to solve the ODEs
        def solve_steadystate_ode(pi_row):
            sol = diffeqsolve(
                ode_term,  # ODETerm
                solver,  # Solver
                t0,  # Start time
                t1,  # End time
                dt0,  # Time step
                tuple(initial_conditions[pi_row[0]][:-1]),  # Initial conditions
                tuple(parameters[pi_row[1]][:-1]),  # Parameters
                max_steps=max_steps,
                saveat=saveat,
                stepsize_controller=stepsize_controller,
                event=Event(steady_state_event()),
                throw=False,
            )
            # Concatenate the steady state values, time point and event mask
            sol = jnp.concatenate(
                [
                    jnp.array(sol.ys).T,
                    jnp.array(sol.ts[-1])[None, None],
                    jnp.array(sol.event_mask, dtype=jnp.float32)[None, None],
                ],
                axis=1,
            )
            # Add the initial condition and parameter numbers
            sol = jnp.column_stack(
                [
                    sol,
                    initial_conditions[pi_row[0], -1],
                    parameters[pi_row[1], -1],
                ]
            )
            return sol

        # Return the functions to solve the ODEs and format the steady state solutions
        return solve_steadystate_ode

    # If the number of time steps to save is not None
    else:
        # Function to solve the ODEs
        def solve_timeseries_ode(pi_row):
            sol = diffeqsolve(
                ode_term,  # ODETerm
                solver,  # Solver
                t0,  # Start time
                t1,  # End time
                dt0,  # Time step
                tuple(initial_conditions[pi_row[0]][:-1]),  # Initial conditions
                tuple(parameters[pi_row[1]][:-1]),  # Parameters
                max_steps=max_steps,
                saveat=saveat,
                stepsize_controller=stepsize_controller,
                throw=False,
            )
            # Concatenate the time series values and time points
            sol = jnp.concatenate(
                [jnp.array(sol.ys).T, jnp.array(sol.ts)[:, None]], axis=1
            )
            # Add the initial condition and parameter numbers
            sol = jnp.column_stack(
                [
                    sol,
                    jnp.repeat(initial_conditions[pi_row[0], -1], len(sol)),
                    jnp.repeat(parameters[pi_row[1], -1], len(sol)),
                ]
            )
            return sol

        # Return the functions to solve the ODEs and format the time series solutions
        return solve_timeseries_ode


# Functiont to run the simulation for a given topo file
def topo_simulate(
    topo_file,
    replicate_dir,
    initial_conditions,
    parameters,
    t0=0.0,
    tmax=200.0,
    dt0=0.01,
    tsteps=None,
    rel_tol=1e-5,
    abs_tol=1e-6,
    max_steps=2048,
    batch_size=10000,
    ode_term_dir=None,
):
    """
    Simulates the ODE system defined by the topology file and saves the results in the replicate directory. The ode system is loaded as a diffrax ode term and the initial conditions and parameters are passed as jax arrays. The simulation is run for the specified time range and time steps and the results are saved in parquet format in the replicate directory.

    Parameters
    ----------
    topo_file : str
        Path to the topology file.
    replicate_dir : str
        Directory where the replicate results will be saved. The ODETerm file will be checked in the parent directory of this directory.
    initial_conditions : pd.DataFrame
        DataFrame containing the initial conditions.
    parameters : pd.DataFrame
        DataFrame containing the parameters.
    t0 : float, optional
        Initial time for the simulation. Default is 0.0.
    tmax : float, optional
        Maximum time for the simulation. Default is 100.0.
    dt0 : float, optional
        Initial time step size. Default is 0.1.
    tsteps : list, optional
        List of time steps at which to save the results. Default is None.
    rel_tol : float, optional
        Relative tolerance for the ODE solver. Default is 1e-5.
    abs_tol : float, optional
        Absolute tolerance for the ODE solver. Default is 1e-6.
    max_steps : int, optional
        Maximum number of steps for the ODE solver. Default is 2048.
    batch_size : int, optional
        Batch size for processing combinations of initial conditions and parameters. Default is 10000.
    ode_term_dir : str, optional
        Directory where the ODE system file is located. Default is None. If None, the parent directory of the replicate directory is assumed to contain the ODE system file. The ODE system file should be named as the topo file with the .py extension.


    Returns:
    -------
    pd.DataFrame
        DataFrame containing the solutions of the ODE system.
    """
    # Get the name of the topo file
    topo_name = topo_file.split("/")[-1].split(".")[0]
    # Making sure to remove the trailing slash from the replicate directory path if it exists
    replicate_dir = replicate_dir.rstrip("/")
    # Check if the ode term directory is None
    if ode_term_dir is None:
        # Getting the parent directory of the replicate directory to get the ODE system file
        simul_dir = os.path.dirname(replicate_dir)
        print(f"Loading ODE system from: {simul_dir}")
    else:
        # Making sure to remove the trailing slash from the ode term directory path if it exists
        simul_dir = ode_term_dir.rstrip("/")
    # Load the ODE system as a diffrax ode term
    ode_term = load_odeterm(topo_name, simul_dir)
    # Getting the intial conditions dataframe column names
    ic_columns = initial_conditions.columns
    # Converting the initial conditions and parameters to jax arrays
    initial_conditions = jnp.array(initial_conditions.to_numpy())
    parameters = jnp.array(parameters.to_numpy())
    # Get the combinations of initial conditions and parameters
    icprm_comb = _gen_combinations(len(initial_conditions), len(parameters))
    print(f"Number of combinations to simulate: {len(icprm_comb)}")
    # Processing the time steps
    if tsteps is None:
        saveat = SaveAt(t1=True)
        print(
            f"Running steady state simulations for replicate: {replicate_dir.split('/')[-1]}"
        )
    else:
        # Checking if the time steps are in the correct format
        tsteps = sorted(tsteps)
        # If the final step more than tmax, make tmax equal to the final step
        if tsteps[-1] != tmax:
            tmax = tsteps[-1]
        # Convert the time steps to a jax array
        tsteps = jnp.array(tsteps)
        # Make the saveat object be the time steps
        saveat = SaveAt(ts=tsteps)
        print(
            f"Running time series simulations for replicate: {replicate_dir.split('/')[-1]}"
        )
    # Specifying the PID controller for the step sizes
    stepsize_controller = PIDController(rtol=rel_tol, atol=abs_tol)
    # Get the functions to solve the ODE system
    solveode_fn = parameterise_solveode(
        ode_term,
        Tsit5(),
        t0,
        tmax,
        dt0,
        saveat,
        stepsize_controller,
        max_steps,
        initial_conditions,
        parameters,
    )
    # Jit compile the solveode function
    solveode_fn = jit(solveode_fn)
    # # # Solve for one combination of initial conditions and parameters
    # sol = solveode_fn(icprm_comb[0])
    # print(sol)
    # Create an empty array to store the solutions
    if saveat.subs.ts is None:
        solution_matrix = np.zeros(
            (len(icprm_comb), initial_conditions.shape[1] + 3), dtype=np.float32
        )
    else:
        solution_matrix = np.zeros(
            (len(icprm_comb) * len(saveat.subs.ts), initial_conditions.shape[1] + 2),
            dtype=np.float32,
        )
    # Defining the length of the time steps to properly index the solution matrix
    len_tsteps = len(saveat.subs.ts) if saveat.subs.ts is not None else 1
    # Iterate over the combinations array in batches
    for ip in range(0, len(icprm_comb), batch_size):
        # print(ip)
        # Get the chunk of the combinations array
        icprm_chunk = icprm_comb[ip : ip + batch_size]
        # vmap the solveode function over the chunk of the combinations array
        sols = vmap(solveode_fn)(icprm_chunk)
        # Vertically stack the solutions
        sols = jnp.vstack(sols)
        # Round the solutions to 4 decimal places
        sols = jnp.round(sols, 4)
        # print(ip * len(saveat.subs.ts), (ip * len(saveat.subs.ts)) + len(sols))
        # Add the solutions to the solution matrix at the correct index
        solution_matrix[ip * len_tsteps : (ip * len_tsteps) + len(sols)] = np.array(
            sols
        )
    # Convert the solution matrix to a dataframe
    # REmoving the InitCondNum column from the initial conditions
    if saveat.subs.ts is None:
        solution_matrix = pd.DataFrame(
            solution_matrix,
            columns=ic_columns[:-1].tolist()
            + ["Time", "SteadyStateFlag", "InitCondNum", "ParamNum"],
        )
        # Find number of steady state solutions
        # print(solution_matrix["SteadyStateFlag"].value_counts())
        # # Make the steady state flag, init cond and param num columns into integers
        # solution_matrix[["SteadyStateFlag", "InitCondNum", "ParamNum"]] = (
        #     solution_matrix[["SteadyStateFlag", "InitCondNum", "ParamNum"]].astype(int)
        # )
        # # Save the solution matrix as a parquet file
        # solution_matrix.to_parquet(
        #     f"{replicate_dir}/{topo_name}_steadystate_solutions.parquet", index=False
        # )
    else:
        solution_matrix = pd.DataFrame(
            solution_matrix,
            columns=ic_columns[:-1].tolist() + ["Time", "InitCondNum", "ParamNum"],
        )
        # # Make the init cond and param num columns into integers
        # solution_matrix[["InitCondNum", "ParamNum"]] = solution_matrix[
        #     ["InitCondNum", "ParamNum"]
        # ].astype(int)
        # Find number of last time points
        # print(solution_matrix["Time"].value_counts())
        # Save the solution matrix as a parquet file
        # solution_matrix.to_parquet(
        #     f"{replicate_dir}/{topo_name}_timeseries_solutions.parquet", index=False
        # )
    return solution_matrix


# Function to run all the replicate simulations for a given topo file
def run_all_replicates(
    topo_file,
    save_dir=".",
    t0=0.0,
    tmax=200.0,
    dt0=0.01,
    tsteps=None,
    rel_tol=1e-5,
    abs_tol=1e-6,
    max_steps=2048,
    batch_size=10000,
    normalize=True,
    discretize=True,
    gk_threshold=1.01,
):
    """
    Run simulations for all replicates of the specified topo file. The initial conditions and parameters are loaded from the replicate folders. The directory structure is assumed to be the same as that generated by the gen_topo_param_files function, with the main directory with the topo file name which has the parameter range file the ODE system file and the replicate folders with the initial conditions and parameters dataframes.

    Parameters
    ----------
    topo_file : str
        Path to the topology file.
    save_dir : str, optional
        Directory where the replicate folders are saved. Defaults to "." i.e current working directory.
    t0 : float, optional
        Initial time for the simulation. Defaults to 0.0.
    tmax : float, optional
        Maximum time for the simulation. Defaults to 100.0.
    dt0 : float, optional
        Initial time step for the simulation. Defaults to 0.1.
    tsteps : int, optional
        Number of time steps for the simulation. Defaults to None.
    rel_tol : float, optional
        Relative tolerance for the simulation. Defaults to 1e-5.
    abs_tol : float, optional
        Absolute tolerance for the simulation. Defaults to 1e-6.
    max_steps : int, optional
        Maximum number of steps for the simulation. Defaults to 2048.
    batch_size : int, optional
        Batch size for the simulation. Defaults to 1000.
    normalize : bool, optional
        Whether to g/k normalise the solutions. Defaults to True.
    discretize : bool, optional
        Whether to discretize the solutions. Defaults to True.
    gk_threshold : float, optional
        A hard threshold value below which the g/k normalised values, if found, will be cliped to 1. Raises an error during discretization if any  g/k normalised values exceed this value. Default is 1.01.

    Returns
    -------
    None
        The results of the simulation are saved in the replicate folders in the specified directory.

    Note
    ----
    The results of the simulation are saved in the replicate folders in the specified directory. If the simulation is time series, the results are saved as `timeseries_solutions.parquet` and if the simulation is steady state, the results are saved as `steadystate_solutions.parquet`.

    Normalization and discretization of the solutions are optional features, but note that discretization requires normalization to be enabled.

    Behavior based on the `discretize` and `normalize` flags:

    - If `discretize=True`, the normalized solutions are discretized, and state counts are saved to `{topo_name}_steadystate_state_counts_{replicate_base}.csv`. This applies only to steady-state simulations.

    - If `discretize=False`, the solutions are normalized but not discretized.

    Effect on the final solution DataFrame:

    - If only `discretize=True`, a `State` column is added. The order of the levels in the state string will be the same as the order of the node columns.
    - If only `normalize=True`, additional columns are added for each node containing the g/k-normalized values. The column names corresponding to the g/k normalised values will have the format "gk_{node_name}".
    - If both flags are set to `True`, both the `state` column and the normalized value columns are included.

    Example
    --------
    Run the simulation for the specified topo file

        >>> run_all_replicates(topo_file, save_dir, t0, tmax, dt0, tsteps, rel_tol, abs_tol, max_steps, batch_size)
    """
    # # Cheking if discretize is True and normalise is False, if so turn normalise to True
    # if discretize and not normalize:
    #     normalize = True
    # Get the name of the topo file
    topo_name = topo_file.split("/")[-1].split(".")[0]
    # Get the list of replicate folders
    replicate_folders = sorted(
        [
            folder
            for folder in glob(f"{save_dir}/{topo_name}/*/")
            if os.path.basename(folder.rstrip("/")).isdigit()
        ],
    )
    # Loop through the replicate folders and run the simulation for each replicate
    for replicate_dir in replicate_folders:
        # Getting the base name of the replicate directory
        replicate_base = os.path.basename(replicate_dir.rstrip("/"))
        # Load the initial conditions and parameters dataframes
        init_cond_path = (
            f"{replicate_dir}/{topo_name}_init_conds_{replicate_base}.parquet"
        )
        params_path = f"{replicate_dir}/{topo_name}_params_{replicate_base}.parquet"
        # Read the initial conditions and parameters dataframes
        init_conds = pd.read_parquet(init_cond_path)
        params = pd.read_parquet(params_path)
        # Starting the timer
        start_time = time.time()
        # Run the simulation for the specified topo file and given initial conditions and parameters
        sol_df = topo_simulate(
            topo_file=topo_file,
            replicate_dir=replicate_dir,
            initial_conditions=init_conds,
            parameters=params,
            t0=t0,
            tmax=tmax,
            dt0=dt0,
            tsteps=tsteps,
            rel_tol=rel_tol,
            abs_tol=abs_tol,
            max_steps=max_steps,
            batch_size=batch_size,
        )
        # Ending the timer
        print(f"Time taken for replicate {replicate_base}: {time.time() - start_time}")
        if discretize:
            print("Normalising and Discretising the solutions")
            sol_df, state_counts = gk_normalise_solutions(
                sol_df,
                params,
                discretize=discretize,
                threshold=gk_threshold,
            )
            # Remove g/k columns if normalise is false
            if not normalize:
                sol_df = sol_df.drop(
                    columns=[col for col in sol_df.columns if col.startswith("gk_")]
                )
        elif normalize:
            print("Normalising the solutions")
            # G/k normalise and/or discretize the solution dataframe
            sol_df = gk_normalise_solutions(
                sol_df, params, discretize=discretize, threshold=gk_threshold
            )
        else:
            pass
        # Check if the time seires is given or not to name the solution file
        if tsteps is None:
            # Save the solution dataframe
            sol_df.to_parquet(
                f"{replicate_dir}/{topo_name}_steadystate_solutions_{replicate_base}.parquet"
            )
            if discretize:
                state_counts.to_csv(
                    f"{replicate_dir}/{topo_name}_steadystate_state_counts_{replicate_base}.csv",
                    index=False,
                    sep="\t",
                )
        else:
            # Read the solution dataframe
            sol_df.to_parquet(
                f"{replicate_dir}/{topo_name}_timeseries_solutions_{replicate_base}.parquet"
            )
        print(f"Simulation completed for replicate: {replicate_base}\n")
        # # break  ##################################
    return None


# Function to g/k normalise the solution dataframe
def gk_normalise_solutions(sol_df, param_df, threshold=1.01, discretize=False):
    """
    Normalises the solutions in the solution dataframe using the maximum production rate (G) and degradation rate (k) parameters of the individual nodes in the parameter sets.

    Parameters
    ----------
    sol_df : pd.DataFrame
        DataFrame containing the solutions with a `ParamNum` column to join with param_df.
    param_df : pd.DataFrame
        DataFrame containing the parameters with `Prod_` and `Deg_` columns for each node.
    threshold : float, optional
        A threshold value below which the g/k normalised values, if found, will be cliped to 1. Raises an error during discretization if any  g/k normalised values exceed this value. Default is 1.01.
    discretize : bool, optional
        Whether to discretise the solutions. Default is True.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the normalised and discretised solutions.
    pd.DataFrame
        DataFrame containing the counts of each state in the normalised solutions.

    Example
    -------
    First, run the simulation then, normalise and discretise the solutions

        >>> sol_df, state_counts = gk_normalise_solutions(sol_df, params)

    Here, the `sol_df` is the solution dataframe and `params` is the parameter dataframe.
    The `state_counts` dataframe contains the counts of each state in the normalised solutions. The returned `sol_df` is the normalised and discretised solution dataframe.
    """
    # Get the node columns from the solution dataframe
    node_cols = [col.replace("Prod_", "") for col in param_df.columns if "Prod_" in col]
    # Get the production and degradation columns
    prod_cols = [f"Prod_{node}" for node in node_cols]
    deg_cols = [f"Deg_{node}" for node in node_cols]
    # Get the gk columns
    gk_cols = [f"gk_{node}" for node in node_cols]
    # Compute the gk values
    param_df[gk_cols] = param_df[prod_cols].values / param_df[deg_cols].values
    # Join the solution dataframe with the parameter dataframe on the 'ParamNum' column
    norm_df = sol_df.merge(param_df[gk_cols + ["ParamNum"]], on="ParamNum")
    # Divide the node columns by the gk columns
    norm_df[gk_cols] = norm_df[node_cols].values / norm_df[gk_cols].values
    if discretize:
        # Select the node columns and discretise the solutions
        norm_df = pd.concat(
            [norm_df, discretise_solutions(norm_df[gk_cols], threshold)], axis=1
        )
        # Get the State columns and return the state counts
        state_counts = norm_df["State"].value_counts()
        # Make the State column as a column
        state_counts = state_counts.reset_index()
        # Rename the columns
        state_counts.columns = ["State", "Count"]
        # Return the normalised and discretised solution dataframe
        return norm_df, state_counts
    else:
        return norm_df


# Function to process the output to discretise the data
def discretise_solutions(norm_df, threshold=1.01):
    """
    Discretises the solutions in a g/k normalized DataFrame based on histogram peaks and minima.

    Parameters
    ----------
    norm_df : pd.DataFrame
        DataFrame containing normalized values to be discretized. It should include only the g/k normalized columns, as the presence of other columns may lead to spurious results.
    threshold : float, optional
        A threshold value below which the g/k normalised values, if found, will be cliped to 1. Raises an error during discretization if any  g/k normalised values exceed this value. If the parameter sets are in such a way that the maximum possible expression of the node is not production/degradation, then the threshold value needs to be adjusted accordingly. Default is 1.01.

    Returns
    -------
    pd.Series
        A Series with the name "State" containing the discrete state labels for each row in the input DataFrame. The order of the labels in the state string is the same as the one input column order.

    Raises
    ------
    UserWarning
        If any value in the DataFrame exceeds the specified threshold value.

    Example
    -------
    Given a normalized DataFrame ```norm_df```, discretise the values

        >>> lvl_df = discretise_solutions(norm_df)

    The normalized solution DataFrame contains values of the nodes between 0 (lowest) and 1 (highest).
    The returned `lvl_df` will have discrete state labels for each row in the input DataFrame.

    Raises a warning if any node values exceed the threshold value. This can occur when a node starts with a value higher than its g/k ratio and the simulation is stopped before reaching steady state, even though the value is approaching the correct limit.

    In time-series simulations with discretization, similar warnings may occur if initial conditions or intermediate values temporarily exceed the g/k threshold. Additionally, it is important to ensure that the time points and solver tolerance settings are appropriately configured, as improper settings can lead to NaN values in the time series.

    For steady-state simulations, increasing the solver's relative and absolute tolerances can improve convergence and reduce such warnings by allowing the simulation to more accurately reach the true steady state.

    """
    # Flatten the numeric part of the dataframe (columns 4 onwards)
    flat = norm_df.values.flatten()

    # Add dummy values to ensure boundary peaks are detected
    dummy_low = np.full(10, -0.1)
    dummy_high = np.full(10, 1.1)
    data_with_dummy = np.concatenate([dummy_low, flat, dummy_high])

    # Compute histogram over 120 bins
    flat_hist, bin_edges = np.histogram(data_with_dummy, bins=120)

    # Define threshold for peak and minima detection (1% of data length)
    peak_detection_threshold = int(len(flat) * 0.01)

    # Detect peaks in the histogram
    peaks, _ = find_peaks(
        flat_hist,
        height=peak_detection_threshold,
        threshold=peak_detection_threshold,
        prominence=peak_detection_threshold,
        distance=int(len(bin_edges) * 0.1),
    )
    maxima_bins = bin_edges[peaks]

    # Detect minima by inverting the histogram (skip first and last bin)
    minima_indices, _ = find_peaks(
        np.max(flat_hist) - flat_hist[1:-1], height=np.max(flat_hist) * 0.99
    )
    # Adjust indices (because we skipped the first bin)
    minima_indices = [idx + 1 for idx in minima_indices]

    # Initialize minima bins with the boundaries 0.0 and 1.0
    minima_bins = [0.0, 1.0]

    # For multiple peaks, find a representative (median) minimum between each adjacent pair
    if len(maxima_bins) > 1:
        for i in range(1, len(maxima_bins)):
            # Candidate minima between two successive peaks
            candidates = [
                bin_edges[idx]
                for idx in minima_indices
                if maxima_bins[i - 1] <= bin_edges[idx] <= maxima_bins[i]
            ]
            if candidates:
                median_candidate = candidates[len(candidates) // 2]
                minima_bins.append(median_candidate)
    # If only one peak exists, take the median minimum if available
    elif minima_indices:
        median_candidate = bin_edges[minima_indices[len(minima_indices) // 2]]
        minima_bins.append(median_candidate)

    # Ensure the bin edges are sorted
    minima_bins = np.sort(minima_bins)

    # Clip the values to 1.0 as due to small numerical errors, some values may be slightly above 1.0
    norm_df = norm_df.mask((norm_df > 1) & (norm_df < threshold), 1.0)

    # If any value is found to be higher than the hard threshold, raise an error
    if (norm_df > threshold).any().any():
        warnings.warn(
            f"Some g/k values exceed the hard threshold of {threshold}, check your input data. Potential reason if running with default values would be solver accuracy - lower rel_tol and abs_tol values.",
            category=UserWarning,
        )

    # Use vectorized binning to assign each value to a discrete level
    # The number of levels is len(minima_bins)-1 (each interval defines one level)
    lvl_df = norm_df.apply(
        lambda col: pd.cut(
            col,
            bins=minima_bins,
            labels=range(len(minima_bins) - 1),
            include_lowest=True,
        )
    )
    lvl_df = lvl_df.add_prefix("Lvl_")

    # Create a 'State' column by concatenating the discrete levels as strings
    lvl_df["State"] = "'" + lvl_df.astype(str).apply("".join, axis=1) + "'"

    return lvl_df["State"]


if __name__ == "__main__":
    pass
    # # Specify the number of cores to use
    # numCores = 15
    # print(f"Number of cores: {numCores}")
    # # Topo file directory
    # topo_dir = "TOPOS"
    # # Specify the root folder where the generated parameter files and then the simulation files will be saved
    # sim_save_dir = "SimulResults"
    # # Make the directories to store the results
    # os.makedirs(sim_save_dir, exist_ok=True)
    # # Get the list of all the topo files
    # topo_files = sorted(glob(f"{topo_dir}/*.topo"))
    # # Remove topo files with 50N
    # topo_files = [topo_file for topo_file in topo_files if "50N" not in topo_file]
    # print(topo_files)
    # print(f"Number of topo files: {len(topo_files)}")
    # # Specify the number of replicates required
    # num_replicates = 3
    # # Specify the number of parameters required
    # num_params = 10000
    # # Specify the number of initial conditions required
    # num_init_conds = 100
    # # Print the number of replicates, parameters and initial conditions
    # print(f"Number of replicates: {num_replicates}")
    # print(f"Number of parameters: {num_params}")
    # print(f"Number of initial conditions: {num_init_conds}\n")
    # # # Start the pool of worker processes
    # # pool = Pool(int(numCores))
    # # # # Parllelise the generation of the parameter and inital condition files
    # # # print("Generating Parameter and Initial Condition files...")
    # # # pool.starmap(
    # # #     gen_topo_param_files,
    # # #     [
    # # #         (
    # # #             topo_file,
    # # #             sim_save_dir,
    # # #             # sim_ode_dir,
    # # #             num_replicates,
    # # #             num_params,
    # # #             num_init_conds,
    # # #         )
    # # #         for topo_file in topo_files
    # # #     ],
    # # # )
    # # # print("Parameter and Initial Condition files generated.\n")
    # # # # Close the pool of workers
    # # # pool.close()
    # # Loop through the topo files and load the ODE system
    # for topo_file in topo_files:
    #     # Generate the parameter files, intial condition files and the directory structure
    #     gen_topo_param_files(
    #         topo_file,
    #         sim_save_dir,
    #         num_replicates,
    #         num_params,
    #         num_init_conds,
    #         sampling_method="Sobol",
    #     )
    #     # Call the function to run the simulation for the specified topo file - Time series
    #     run_all_replicates(
    #         topo_file,
    #         sim_save_dir,
    #         tsteps=jnp.array([25.0, 75.0, 100.0]),
    #         max_steps=2048,
    #     )
    #     # # Call the function to run the simulation for the specified topo file - Steady state
    #     run_all_replicates(
    #         topo_file,
    #         sim_save_dir,
    #         normalize=False,
    #         discretize=True
    #     )
    #     break

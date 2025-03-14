# Gene Regulatory Interaction Network Simulator (GRiNS)

A Python library for simulating gene regulatory networks (GRNs) using parameter-agnostic frameworks like RACIPE and Ising formalism, with GPU acceleration and efficient ODE solving.

Modeling gene regulatory networks (GRNs) is essential for understanding cellular processes, but parameterizing these networks becomes increasingly difficult as they scale. This Python library provides a simulation framework that unifies **parameter-agnostic** approaches, including **RACIPE** and **Ising formalism**, into a single, flexible tool.  

## Key Features  

- **Simulation Frameworks**: Supports both **ODE-based** (RACIPE) and **coarse-grained** (Ising formalism) methods for studying GRN dynamics.  
- **Parameter-Agnostic Modeling**: Translates network topology into mathematical models without requiring detailed parameter tuning.  
- **Scalable Computation**: Uses the [Jax](https://github.com/jax-ml/jax) ecosystem for GPU acceleration and [Diffrax](https://github.com/patrick-kidger/diffrax) for efficient ODE solving.  
- **Data Processing Tools**: Provides **normalization and discretization** functions to standardize simulation outputs for downstream analysis.  

<!-- You can access the full documentation, including installation instructions, usage examples, and detailed explanations of the simulation frameworks, at [MoltenEcdysone09.github.io/GRiNS](https://MoltenEcdysone09.github.io/GRiNS) -->

## Installation  

### GPU Version Installation (Recommended)  

For optimal performance, it is recommended to install the GPU-accelerated version of the library. This version leverages CUDA for faster computations, making it well-suited for large-scale simulations. If you have a compatible NVIDIA GPU, install the library with:  

```bash
pip install grins[cuda12]
```

### CPU Version Installation

If you do not have a compatible GPU, you can install the CPU version instead:

```bash
pip install grins
```

Compared to the GPU version, the CPU version will be slower, especially for large simulations.

## Citation

<!-- 
## Commands

* `mkdocs new [dir-name]` - Create a new project.
* `mkdocs serve` - Start the live-reloading docs server.
* `mkdocs build` - Build the documentation site.
* `mkdocs -h` - Print help message and exit.

## Project layout

    mkdocs.yml    # The configuration file.
    docs/
        index.md  # The documentation homepage.
        ...       # Other markdown pages, images and other files. -->

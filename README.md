# Gene Regulatory Interaction Network Simulator (GRiNS)

A Python library for simulating gene regulatory networks (GRNs) using parameter-agnostic frameworks like RACIPE and Ising formalism, with GPU acceleration and efficient ODE solving.

Modeling gene regulatory networks (GRNs) is essential for understanding cellular processes, but parameterizing these networks becomes increasingly difficult as they scale. This Python library provides a simulation framework that unifies **parameter-agnostic** approaches, including **RACIPE** and **Ising formalism**, into a single, flexible tool.  

## Key Features  

- **Simulation Frameworks**: Supports both **ODE-based** (RACIPE) and **coarse-grained** (Ising formalism) methods for studying GRN dynamics.  
- **Parameter-Agnostic Modeling**: Translates network topology into mathematical models without requiring detailed parameter tuning.  
- **Scalable Computation**: Uses the [Jax](https://github.com/jax-ml/jax) ecosystem for GPU acceleration and [Diffrax](https://github.com/patrick-kidger/diffrax) for efficient ODE solving.  
- **Data Processing Tools**: Provides **normalization and discretization** functions to standardize simulation outputs for downstream analysis.  

<p align="center">
  <img src="https://github.com/MoltenEcdysone09/GRiNS/blob/main/GRINS_Website_New.png?raw=true" alt="Overview of the simulation frameworks in GRiNS. GRiNS includes implementations of Random Circuit Perturbation (RACIPE) for continuous ODE-based modeling and Ising Boolean formalism for discrete-state simulations."/>
</p>

## Documentation

You can access the full documentation, including installation instructions, usage examples, and detailed explanations of the simulation frameworks, at [MoltenEcdysone09.github.io/GRiNS](https://MoltenEcdysone09.github.io/GRiNS)

## Case Studies

The code to reproduce the analysis done in the manuscript can been accessed in the [Case Studies]("/CaseStudies/") folder.

## Installation  

### GPU Version Installation (Recommended)  

For optimal performance, it is recommended to install the GPU-accelerated version of the library. This version leverages CUDA for faster computations, making it well-suited for large-scale simulations. If you have a compatible NVIDIA GPU (refer to [Jax Installation](https://docs.jax.dev/en/latest/installation.html)), install the library with:  

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

Please cite this package if you have used it.
```
@misc{harlapur2025grinspythonlibrarysimulating,
      title={GRiNS: A Python Library for Simulating Gene Regulatory Network Dynamics}, 
      author={Pradyumna Harlapur and Harshavardhan B V and Mohit Kumar Jolly},
      year={2025},
      eprint={2503.18356},
      archivePrefix={arXiv},
      primaryClass={q-bio.QM},
      url={https://arxiv.org/abs/2503.18356}, 
}
```

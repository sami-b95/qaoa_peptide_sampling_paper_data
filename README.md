This is the data repository for the paper [arXiv:2204.01821](https://arxiv.org/abs/2204.01821).

We provide 3 notebooks reproducing the figures of the paper. These rely on heavier data available [there](.https://zenodo.org/record/6563433/files/data.tar.xz?download=1). To use the notebooks, download this archive and uncompress to the root of the repository.

Execution of the notebooks requires libraries `Pandas` and `RDKit` (the latter, only for [`small_peptide_sampling_probability_distribution_analysis.ipynb`](./`small_peptide_sampling_probability_distribution_analysis.ipynb`)).

- The [`self_avoiding_walk_analysis.ipynb`](./self_avoiding_walk_analysis.ipynb) reproduces all figures for the self-avoiding walk section of the paper.
- The [`small_peptide_sampling_energy_analysis.ipynb`](./small_peptide_sampling_energy_analysis.ipynb) notebook reproduces the figures analyzing the energy achieved by QAOA when sampling small peptides.
- The [`small_peptide_sampling_probability_distribution_analysis.ipynb`](./`small_peptide_sampling_probability_distribution_analysis.ipynb`) notebook reproduces the figure analyzing the distribution of these conformations more finely than through their energies.
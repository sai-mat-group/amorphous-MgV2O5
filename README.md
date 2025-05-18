
# Exploration of Amorphous V₂O₅ as a Cathode for Magnesium Batteries

This repository contains **all data and analysis files** used in the research article, **“Exploration of Amorphous V₂O₅ as Cathode for Magnesium Batteries”**  authored by Vijay Choyal, Debsundar Dey, and Gopalakrishnan Sai Gautam. The work is currently under review, but a preprint is available on [arXiv](https://arxiv.org). 

---

## 🧪 About the Study

The primary goal of this work is to identify promising **amorphous cathode materials for magnesium batteries**, using advanced atomistic simulations.

This study combines:

- Ab Initio Molecular Dynamics **Melt-Quench Simulations**: To generate amorphous V<sub>2</sub>O<sub>5</sub> structures.
- **Moment Tensor Potential (MTP)** construction: A machine-learned interatomic potential trained on AIMD data.
- **Voltages**: Effect of **long-range (dis)order** on average magneisum intercalation voltage.
- **Classical Molecular Dynamics** (CMD) simulations: To investigate **dynamical properties** and quantify magnesium transport in amorphous V<sub>2</sub>O<sub>5</sub>.
- **PyLRO**: **Quantify long-range order** (LRO) in the amorphous structures.
---

## 📁 Repository Structure

- `2_4_6_msd_4ns/`: Contains mean square displacement (MSD) analysis scripts, results and dump files from `LAMMPS' (i.e., CMD simulations).
- `Active_learning_files/`: Data and scripts for active learning during MTP construction.
- `MTP_construction/`: Training and testing Data for MTP and the final optimised MTP potential. 
- `pylro_calc/`: Files related to PyLRO calculations for LRO quantification.
- `videos_and_Gifs/`: Visualizations of the structures and magneisum hopping behavior.

---

## ⚙️ Installation and Usage

For a detailed theoretical explanation of PyLRO and its application, see:  
📄 [DOI: 10.1063/5.0244012](https://doi.org/10.1063/5.0244012)

For installing and using MTP and its interface with LAMMPS, please refer to the following manuals:

- 📘 [MTP Manual (MLIP 2)](https://gitlab.com/ashapeev/mlip-2-paper-supp-info/-/blob/master/manual.pdf)
- ⚙️ [LAMMPS–MLIP Interface Setup](https://gitlab.com/ashapeev/interface-lammps-mlip-2/-/blob/master/README.md)
---

## 📬 Citation

If you use any part of this dataset or constructed potential or scripts, please drop a citation to our work at [arXiv](https://arxiv.org). We will also update the DoI of the manuscript here once it is published.

---


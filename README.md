# amorphous-MgV2O5
Exploration of amorphous MgV2O5 as a cathode for Mg batteries
# Exploration of Amorphous V₂O₅ as a Cathode for Magnesium Batteries

This repository contains **all data and analysis files** used in the research:

**_“Exploration of Amorphous V₂O₅ as Cathode for Magnesium Batteries”_**  
**Authors**: Vijay Choyal, Debsundar Dey, and Gopalakrishnan Sai Gautam  
_(Manuscript currently under review. Link will be updated soon.)_

---

## 🧪 About the Study

The primary goal of this work is to identify promising **amorphous cathode materials for magnesium-ion batteries**, using advanced atomistic simulations.

This study combines:

- **AIMD Melt-Quench Simulations**: To generate amorphous structures.
- **Moment Tensor Potential (MTP)** construction: A machine-learning potential trained on AIMD data.
- **Classical Molecular Dynamics (MD)** simulations: To investigate **dynamical properties**.
- **PyLRO**: Used to **quantify long-range order (LRO)** in the amorphous systems.

For a detailed theoretical explanation of PyLRO and its application, see:  
📄 [DOI: 10.1063/5.0244012](https://doi.org/10.1063/5.0244012)

---

## 📁 Repository Structure

- `1c00134/`: (Details TBD)
- `2_4_6_msd_4ns/`: Contains MSD analysis scripts and results and dumb files.
- `Active_learning_files/`: Data and scripts for active learning during MTP construction.
- `MTP_construction/`: Training data, Testing Data for MTP and constructed MTP potential 
- `pylro_calc/`: Files related to PyLRO calculations for order quantification.
- `videos_and_Gifs/`: Visualizations of the structures and dynamics.

---

## ⚙️ Installation and Usage

For installing and using the **Moment Tensor Potential (MTP)** and its interface with LAMMPS, please refer to the following detailed manuals:

- 📘 [MTP Manual (MLIP 2)](https://gitlab.com/ashapeev/mlip-2-paper-supp-info/-/blob/master/manual.pdf)
- ⚙️ [LAMMPS–MLIP Interface Setup](https://gitlab.com/ashapeev/interface-lammps-mlip-2/-/blob/master/README.md)

---

## 📬 Citation

If you use any part of this dataset or code, please cite the associated manuscript (link will be added once available).

---


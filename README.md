# Antigen-Antibody Residue Hotspot Analysis


![image](https://github.com/user-attachments/assets/4c9cf5d1-f588-4c2d-975e-7c9b0cafd2b7)





[![MIT License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

This repository provides a Python-based analysis tool for identifying antigen-antibody residue hotspots using Ridge regression.

## **Required Inputs**
To use this analysis, you need to provide the following CSV file:

1. **`interaction_energy.csv`**  
   - This file should contain the **non-bonded interaction energies** from a **non-equilibrium pulling simulation (Can be easily obtained using gmx energy using indexing (antigen-antibody)**.
   - The energy is calculated as the **sum of Van der Waals and short-range Coulomb interactions**.
   - **Example format**:
     ```
     Frame,AverageInteractionEnergy
     1,-45.2
     2,-48.7
     3,-50.1
     ```
   - The first column (`Frame`) represents simulation frames, and the second column (`AverageInteractionEnergy`) represents the computed interaction energy.


## **Installation and Running the Analysis**

1. **Clone the repository**:
   ```bash
   git clone https://github.com/<your-username>/antigen-antibody-residue-hotspot.git


2. ```bash
   cd antigen-antibody-residue-hotspot

3. ```bash
   python analysis_script.py

5. Sample outputs are shown here, (The analysis is more robus with multiple replicate results)

![image](https://github.com/user-attachments/assets/506f8cb2-7ce4-44fb-90fb-eb2821aee53b)

   
![image](https://github.com/user-attachments/assets/8470cd80-f8a6-4950-bbb7-9fdc32c8bdd1)

![image](https://github.com/user-attachments/assets/05b17a10-b887-4d72-824f-c057a436221c)

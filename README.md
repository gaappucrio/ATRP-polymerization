# ATRP-polymerization
Parameter Estimation and Kinetic Monte Carlo Simulation of Styrene and n‐Butyl Acrylate Copolymerization through ATRP.

Access the paper at: https://doi.org/10.1021/acs.iecr.1c00943

# Kinetic Monte Carlo Simulation: Styrene and n-Butyl Acrylate Copolymerization via ATRP

This repository contains the Python source code for a Kinetic Monte Carlo (MC) simulation focused on describing the copolymerization of styrene and n-butyl acrylate through Atom Transfer Radical Polymerization (ATRP).

## 🧪 About the Project

The dynamic model simulates the step-by-step kinetic mechanism of the reaction, tracking the growth of polymer chains. It tracks average polymer weight properties, including the number-average molecular weight (Mn) and the polydispersity index (PDI), which are essential for predicting the final product's properties.

## 🗂 File Structure

This repository consists of two main files:

* **`MC_Estireno_ButilAcrilato_V10.py`**: The main script. It contains the initial conditions (such as control volume), kinetic constants from the literature for initiation, propagation, chain transfer, and termination processes, and executes the main stochastic simulation loop.
* **`lib_MC.py`**: An auxiliary library imported by the main file. It contains functions for calculating the molecular weight (`Mn`) and a simple Graphical User Interface (GUI) built with `tkinter` to alert the user when the simulation is complete.

## ⚙️ How the Algorithm Works

The mathematical logic bypasses continuous differential equations and relies on Gillespie's "direct method" for stochastic modeling. In practice, the program works as follows:

1. **Rate Calculation:** It evaluates the reaction rates for all possible chemical events at any given moment, using the number of molecules and constants converted for the Monte Carlo model.
2. **Time Advancement:** It uses random numbers to determine the exact time until the next reaction occurs (time step tau).
3. **Reaction Selection:** It randomly selects which specific chemical reaction will happen, based on probabilities proportional to the rates calculated in step 1.
4. **Update:** It updates the system matrices (e.g., adding a monomer to a chain, activating a radical) and repeats the cycle until the established final time (`tf = 4.e4`) is reached.

## 🚀 How to Run

**Prerequisites:**
Ensure you have Python 3 installed in your environment, along with the following libraries:
```bash
pip install numpy matplotlib pillow
(Note: tkinter, used for the final pop-up, is usually included in standard Python installations).
```

Step-by-step:

Clone this repository to your local machine.
```bash
git clone https://github.com/gaappucrio/ATRP-polymerization.git
```
Navigate to the project's root folder.

Run the simulator in your terminal:

```bash
python MC_Estireno_ButilAcrilato_V10.py
```
The computational time may vary depending on your machine. Upon completion, the terminal will display the total execution time, and a pop-up window will confirm the end of the simulation.

📊 Understanding the Outputs
After the routine finishes, the script automatically exports four .txt files to the same directory. These files contain raw data ready for plotting or statistical analysis:

conversao.txt: Main data table. It lists the time progression, fractional conversion of monomers, Mn, weight-average molecular weight (Mw), and PDI.

distribuicao.txt: A compiled distribution of the formed polymer's molecular weight. It includes the chain length, observed quantity, and calculated masses and fractions.

distribuicao_total.txt: A matrix containing the exact length record of every final inactive chain.

log.txt: A simple log file recording the total execution time and the simulated control volume.

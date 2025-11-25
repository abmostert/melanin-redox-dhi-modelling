# melanin-redox-dhi-modelling
MATLAB code used to model redox and radical equilibria in DHI/melanin systems. Includes analytical modelling, nonlinear least-squares fitting, and species simulations across pH. Accompanies a published research study.

# DHI Redox Modelling in MATLAB

This repository contains the MATLAB code used to model the redox and free-radical behaviour of 5,6-dihydroxyindole (DHI), a melanin precursor. The work combines:

- Analytical expressions for quinone/catechol redox chemistry  
- Literature-derived and estimated thermodynamic parameters (pKa, redox potentials, tautomerisation constants)  
- Non-linear least-squares fitting in MATLAB to match experimental pH-dependent data (e.g. Chio et al., melanin EPR data)

The goal is to understand how quinone–methide tautomerisation and protonation equilibria affect radical concentrations and redox potentials, and how a “simple” DHI model can explain trends in more complex systems such as melanin and polydopamine.

---

## Repository Structure

### 1. Standard redox potential estimation (DHI / DHI-melanin)

These scripts estimate standard and apparent (pH-dependent) one- and two-electron reduction potentials, as well as comproportionation constants, under different modelling assumptions.

- **`StandardRedoxDHImelanin_NoQuinoneMethide.m`**  
  Estimates redox potentials for a DHI/melanin-like system **without** including the quinone–methide tautomer.  
  - Uses literature pKa values and an estimated comproportionation constant.  
  - Computes:
    - K\_comp (comproportionation equilibrium)
    - Deprotonated standard potentials \(E_{SQQ}, E_{Q2SQ}\)
    - Fully protonated potentials for comparison to literature.

- **`StandardRedoxDHImelanin_NoQuinoneImine.m`**  
  Similar analysis, but under an alternative assumption set (no quinone imine).  
  - Scans across different initial fractions of catechol (H\_2Q).  
  - Useful for sensitivity analysis with respect to starting composition.

- **`StandardRedoxDHImelanin.m`**  
  Full model including **quinone–methide tautomerisation**.  
  - Introduces:
    - `Kt` (tautomerisation equilibrium constant)
    - `Kaqi` (quinone methide acid dissociation constant)  
  - Computes:
    - K\_comp including tautomerisation
    - Deprotonated and fully protonated redox potentials
    - Gibbs free energy of the comproportionation reaction (in eV) at different conditions.

These files correspond to the theoretical framework presented in the paper’s methods section and underpin the interpretation of DHI/melanin redox behaviour.

---

### 2. Modelling pH-dependent radical concentration (Chio et al.-type data)

These scripts set up and solve the full equilibrium problem (charge balance, mass balance, electron balance) using **non-linear root finding** (`fsolve`, Levenberg–Marquardt).

Each model solves for:

- System potential \(E\)  
- pH  
- Quinone concentration \([Q]\)

…and then derives the concentrations of all relevant species (H\_2Q, HQ, Q\_2, Q, SQ, HSQ, and, in the extended model, quinone–methide species).

#### Without quinone–methide

- **`Model_Chio_1_Final.m`**  
  Top-level script for the **standard ortho-quinone model**, without any quinone–methide tautomerisation.

  - Defines global parameters:
    - Initial total concentration  
    - Initial catechol fraction  
    - Titrant concentration and volumes  
    - pKa values for DHI  
  - Loops over titrant volume and calls `fsolve` on the model function.
  - Computes:
    - pH vs. total radical concentration (SQ + HSQ) in spins per gram  
    - log concentrations of all species vs. pH  
    - System potential vs. pH  
    - An approximate g-value vs. pH (weighted between HSQ and SQ)

- **`Model_Function_Chio_1_Final.m`**  
  The corresponding function passed to `fsolve`.  
  - Defines:
    - Charge balance  
    - Mass balance (DHI-like moieties)  
    - Electron balance  
  - Expresses species concentrations in terms of the unknowns \(E\), pH, [Q] and equilibrium constants.
  - Returns the vector of residuals \(F(E, \text{pH}, [Q])\) to be minimised.

#### With quinone–methide

- **`Model_Chio_2_Final.m`**  
  Top-level script for the **extended model including quinone–methide** formation and deprotonation.

  - Adds:
    - Tautomerisation constant `Kt`  
    - Quinone–methide pKa `Kaqi`  
    - Additional species: HQI, QI (protonated and deprotonated quinone–methide)
  - As in Model 1, loops over volume, solves with `fsolve`, and stores:
    - pH vs. total radical concentration  
    - Species concentration vs. pH (including quinone–methide forms)  
    - Potential vs. pH  
    - Derived g-value vs. pH

- **`Model_Function_Chio_2_Final.m`**  
  Extended equilibrium model including:
  - Tautomerisation equilibrium  
  - Quinone–methide protonation/deprotonation  
  - Modified mass and electron balances accounting for HQI and QI.

---

## Requirements

- MATLAB (tested originally with a version supporting `fsolve` and `optimset`)  
- Optimization Toolbox (for `fsolve` with Levenberg–Marquardt algorithm)

No external data files are required: all constants are encoded in the scripts, based on literature values and estimated parameters discussed in the paper.

---

## How to Run

1. Clone or download this repository.
2. Open MATLAB and set the working directory to the repo folder.
3. To reproduce standard redox potential calculations:
   - Run `StandardRedoxDHImelanin_NoQuinoneMethide.m`
   - Run `StandardRedoxDHImelanin_NoQuinoneImine.m`
   - Run `StandardRedoxDHImelanin.m`
4. To run the pH-dependent radical concentration models:
   - Run `Model_Chio_1_Final.m` for the standard quinone model.  
   - Run `Model_Chio_2_Final.m` for the model including quinone–methide.
5. Each script generates figures (potential vs. pH, species concentrations vs. pH, total radical concentration vs. pH, etc.) and a `datamatrix` variable you can export for further analysis.

You can adapt initial conditions (e.g. catechol fraction, pKa values, tautomerisation constants) inside the scripts to explore sensitivity to different parameters.

---

## Relation to the Paper

These scripts implement the modelling framework described in the Methods/Theory section of the associated paper on DHI redox chemistry and its implications for melanin and polydopamine. They were used to:

- Compare models with and without quinone–methide tautomerisation  
- Interpret pH-dependent EPR data (e.g. Chio et al.)  
- Estimate thermodynamic quantities such as apparent reduction potentials and the Gibbs free energy of the comproportionation reaction.

**A. Bernarduc Mostert,** *On the free radical redox chemistry of 5,6-dihydroxyindole*, Chemical Physics (Elsevier), (2021), 546, 111158.  
DOI: 10.1016/j.chemphys.2021.111158

---

## How to cite

If you use this code in scientific or academic work, please cite the associated publication:

**A. Bernarduc Mostert,** *On the free radical redox chemistry of 5,6-dihydroxyindole*, Chemical Physics (Elsevier), (2021), 546, 111158.  
DOI: 10.1016/j.chemphys.2021.111158

---

## License

This project is released under the **MIT License**.  
See the `LICENSE` file for details.

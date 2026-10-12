# Transient Thermal Modeling of a Cylindrical Battery in MATLAB

**Course project | Heat Transfer | Sharif University of Technology (2025)**  
**Author:** Sepehr Javanmojarrad  
**Tools:** MATLAB, finite-difference methods

## Overview

This academic project examines transient, one-dimensional radial heat conduction in a cylindrical battery with time-dependent internal heat generation and convective heat loss at the outer surface. The submitted MATLAB scripts implement an explicit finite-difference method and a Crank–Nicolson method to investigate temperature evolution and numerical behavior.

## Files

| File | Description |
| --- | --- |
| `Q1HTProject.m` | Explores the explicit method's behavior at different time steps and reports a trial-and-error stability threshold. |
| `Q2HTProject.m` | Compares explicit and Crank–Nicolson predictions at radial positions of 4 mm and 8 mm; prints temperatures at selected times. |
| `Q3HTProject.m` | Compares radial temperature profiles and produces two-dimensional visualizations at selected times. |

## Running the scripts

1. Open MATLAB.
2. Place the three `.m` files in the same working directory.
3. Run each script separately, for example:

   ```matlab
   Q1HTProject
   Q2HTProject
   Q3HTProject
   ```

Each script defines its own parameters and produces figures. The second script also prints a table. No external MATLAB toolboxes are explicitly called in the source files.

## Numerical-method notes

- The scripts use a 1D radial finite-difference grid, a time-dependent volumetric heat source, and an outer-surface convection condition.
- `Q1HTProject.m` prints a *trial-and-error* critical time-step estimate; that value is **not a verified analytical or spectral stability boundary** and should not be used as such without further checking.
- These are the original course-project scripts, included **without changes**. They have not been packaged as a validated thermal-design tool or independently tested here in MATLAB.

## Scope

This repository documents undergraduate numerical heat-transfer coursework. It is not a battery-management system or an experimentally validated battery model.

## Contact

[LinkedIn](https://www.linkedin.com/in/sepehr-javan-5a0818298/) | [GitHub profile](https://github.com/sprjnm)
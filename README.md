# Advanced Dynamics of Mechanical Systems & Structural Engineering

This repository presents an engineering portfolio focusing on **structural dynamics**, **3D laser scanning experimental modal analysis (EMA)**, and **2D finite element method (FEM) transient moving load simulations**.

---

## Projects Overview

### 1. 🏗️ Experimental Modal Analysis - Milan Citywave Canopy

* **Description**: Dynamic characterization and finite element model updating of a 1:70 scaled aeroelastic model of the **Milan Citywave Canopy** roof structure.
* **Methodology**:
  * 3D Laser Doppler Vibrometry (PSV) velocity measurement over a 202-point spatial grid at the Politecnico di Milano Wind Tunnel (GVPM).
  * Multi-channel $H_1$ Frequency Response Functions ($G_{jk}$) processing.
  * Modal damping estimation ($\zeta_{HP}$, $\zeta_{PD}$) using Half-Power (-3dB) and Phase-Derivative methods.
  * Mode shape extraction via quadrature and spatial correlation against 3D FE models via **Modal Assurance Criterion (MAC)** and **3D 6x6 Cross-MAC matrix**.
* **Folder & Detailed Technical Readme**: `01_Milano_Canopy_Modal_Analysis/`

---

### 2. ⛷️ 2D FEM & Dynamic Moving Load Analysis - Predazzo Olympic Ski Jump Ramp

* **Description**: Complete 2D Euler-Bernoulli beam element modeling and transient dynamic analysis of the **Predazzo Olympic Ski Jump Ramp** for the **Milano-Cortina 2026 Winter Olympic Games**.
* **Methodology**:
  * Global mass ($M$) and stiffness ($K$) matrix assembly incorporating non-uniform column pillars, lattice bracing, curved arch deck ($R=112\text{m}$), and ground spring boundary condition ($k = 5 \times 10^7\text{ N/m}$).
  * Mesh refinement convergence study optimizing element lengths.
  * Eigenmode extraction and natural frequency characterization ($f_1 = 1.56\text{ Hz}$ to $f_4 = 5.29\text{ Hz}$).
  * **Static & Moving Load Simulation**: Static snow load deflection assessment ($q = 7\text{ kN/m}$) and transient dynamic ODE integration (`ode45`) simulating an accelerating athlete ($90\text{ kg}$) traversing the track deck.
* **Folder & Detailed Technical Readme**: `02_Olympic_Ski_Jump_FEM_Dynamics/`

---

### 3. 📐 Theoretical Continuous Beam Dynamics & Parameter Identification

* **Description**: Analytical transverse bending dynamics of slender Euler-Bernoulli cantilever beams and non-linear parameter identification.
* **Methodology**:
  * Analytical PDE solution, boundary matrix determinant roots ($|\det H(\omega)| = 0$), and exact natural frequencies ($f_1 = 4.50\text{ Hz}$ to $f_4 = 154.86\text{ Hz}$).
  * Non-linear complex residual curve-fitting (`lsqnonlin`) absorbing out-of-band inertial ($R^L$) and stiffness ($R^H$) tails.
  * Sensor spatial sensitivity study evaluating antinodal vs nodal proximity placement.
* **Folder & Detailed Technical Readme**: `03_Beam_Vibration_Dynamics/`

---

## Repository Structure

```
ADMS/
├── 01_Milano_Canopy_Modal_Analysis/     # Milan Citywave Canopy 3D Laser EMA & FEM Validation
│   ├── README.md                        # Technical module documentation & data table
│   ├── main.m                           # Master EMA processing & MAC comparison script
│   ├── data/                            # Experimental FRFs, Laser Mesh & FE model (.mat)
│   └── functions/                       # Peak detection, HP & PD damping routines
│
├── 02_Olympic_Ski_Jump_FEM_Dynamics/    # Predazzo Olympic Ski Jump Ramp 2D FEM & Moving Load
│   ├── README.md                        # Technical module documentation & parameters
│   ├── main.m                           # Master FEM, static snow & moving load dynamic script
│   ├── input.inp                        # Geometry input file generated dynamically
│   └── functions/                       # Matrix assembly, 6x6 transformation & solver functions
│
├── 03_Beam_Vibration_Dynamics/          # Theoretical Euler-Bernoulli Beam Dynamics
│   ├── README.md                        # Technical module documentation & equations
│   ├── main_beam_dynamics.m             # Analytical solver & lsqnonlin fitting script
│   └── images/                          # Sensitivity study output figures
│
├── reports/                             # Technical Reports in LaTeX (PDF)
│   ├── Milano_Canopy_Modal_Analysis_Report.pdf
│   ├── Olympic_Ski_Jump_FEM_Report.pdf
│   └── Beam_Dynamics_Theoretical_Report.pdf
│
└── .gitignore                           # Excludes temporary files & course reference slides (docs/)
```

---

## Requirements 

* **Software**: MATLAB R2022b or newer.
* **Toolboxes**: Optimization Toolbox, Signal Processing Toolbox.

To run any simulation, open MATLAB, navigate to the target module directory, and run the master script:

- **Citywave Canopy Modal Analysis**: `01_Milano_Canopy_Modal_Analysis/main.m`
- **Olympic Ski Jump FEM & Moving Load**: `02_Olympic_Ski_Jump_FEM_Dynamics/main.m`
- **Continuous Beam Vibration**: `03_Beam_Vibration_Dynamics/main_beam_dynamics.m`

---

## Technical Reports

For full mathematical derivations, complete experimental setup descriptions, extended FE formulations, and comprehensive discussion of numerical results, please refer to the compiled technical PDF reports located in the `reports/` folder:

* 📄 **Milan Citywave Canopy Modal Analysis Report**: `reports/Milano_Canopy_Modal_Analysis_Report.pdf`
* 📄 **Predazzo Olympic Ski Jump Ramp FEM & Dynamics Report**: `reports/Olympic_Ski_Jump_FEM_Report.pdf`
* 📄 **Continuous Beam Dynamics & System Identification Report**: `reports/Beam_Dynamics_Theoretical_Report.pdf`

---

## Author 

* **Author**: Francesco Valerio Persio Pennesi
* **Institution**: Politecnico di Milano

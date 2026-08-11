# 2D Finite Element Modeling & Moving Load Transient Dynamics of the Predazzo Olympic Ski Jump Ramp

This module presents the 2D finite element structural modeling, mesh convergence optimization, eigenmode extraction, static environmental load assessment, and **moving load dynamic transient simulation** for the **Predazzo Ski Jump Ramp** designed for the **Milano-Cortina 2026 Winter Olympic Games**.

---

## ⛷️ Structural Topology & Technical Data

The ski jump ramp consists of a steel framework supporting two vertical lattice towers, a straight inclined inrun deck, and a circular curved launch arch ($R = 112 \text{ m}$).

### Geometry & Node Coordinates
* **Pillar 1 (Left Support Tower)**: Base at $x = 0 \text{ m}$ ($y = 0 \text{ m}$ to $60 \text{ m}$), ground width $17 \text{ m}$.
* **Pillar 2 (Center Support Tower)**: Positioned at $x = 50.6 \text{ m}$ to $67.6 \text{ m}$.
* **Ramp Inrun & Curved Track**: Point A at $(0, 60)\text{ m}$, transition at $(67.6, 17.5)\text{ m}$, terminal launch tip Point B at $(128.0, 0.0)\text{ m}$.
* **Boundary Conditions**: Clamped ground supports at tower bases; linear vertical spring $k = 5 \times 10^7 \text{ N/m}$ attached at Point B.

### Structural Beam Section Properties
The lattice framework is composed of three distinct structural beam profiles:

| Beam Type | Mass per Length $m$ [kg/m] | Axial Stiffness $EA$ [N] | Bending Stiffness $EJ$ [N$\cdot$m$^2$] | Structural Location |
|:---:|:---:|:---:|:---:|:---|
| **Red Beams** | 122 | $3.21 \times 10^9$ | $2.0 \times 10^8$ | Main column pillars & primary vertical supports |
| **Green Beams** | 80 | $1.74 \times 10^9$ | $4.8 \times 10^7$ | Internal lattice diagonal & horizontal bracing |
| **Blue Beams** | 400 | $2.50 \times 10^9$ | $4.0 \times 10^9$ | Top sliding track deck & primary arch girder |

---

## 🧮 Finite Element Formulation & Dynamic Analysis

1. **Matrix Assembly & Mesh Refinement**:
   * Customized 2D Euler-Bernoulli beam element formulation ($6 \times 6$ element mass $M^e$ and stiffness $K^e$ matrices transformed via direction cosine matrices $\mathbf{T}$).
   * Admissible element length sizing based on maximum frequency $\omega_{\max} = 6 \times 2\pi \text{ rad/s}$ and safety factor $\eta = 4.5$:
     $$L_{\max} = \sqrt{\frac{\pi^2}{\eta \omega_{\max}} \sqrt{\frac{EJ}{m}}}$$
   * Mesh convergence loop verifying natural frequency stabilization within $< 0.1\%$ variance for the first 4 modes.

2. **Eigenmode Extraction**:
   * Free-vibration eigenvalue problem solved for unconstrained DOFs: $[\mathbf{K}_{FF} - \omega^2 \mathbf{M}_{FF}] \mathbf{\phi} = \mathbf{0}$.
   * **Mode 1 ($f_1 = 1.56 \text{ Hz}$)**: First flexural bending mode of the circular curved launch track.
   * **Mode 2 ($f_2 = 3.04 \text{ Hz}$)**: Lateral bending and sway of the main support pillars.
   * **Mode 3 ($f_3 = 4.88 \text{ Hz}$)**: Vertical bending mode of the straight span between Pillar 1 and Pillar 2.
   * **Mode 4 ($f_4 = 5.29 \text{ Hz}$)**: Combined higher-order track flexure.

3. **Rayleigh Damping & Forced FRF Receptance**:
   * Proportional damping matrix $\mathbf{C} = \alpha \mathbf{M} + \beta \mathbf{K}$ with $\alpha = 0.1 \text{ s}^{-1}$ and $\beta = 2 \times 10^{-4} \text{ s}$.
   * Dynamic receptance matrix $\mathbf{H}(\Omega) = [-\Omega^2 \mathbf{M}_{FF} + i \Omega \mathbf{C}_{FF} + \mathbf{K}_{FF}]^{-1}$.
   * Reduced-order validation via **Modal Superposition** (2 modes vs full model).

4. **Static Snow Load Assessment**:
   * Distributed vertical snow load $q = 7 \text{ kN/m}$ applied along the top blue track elements ($p_0 = -7000 \text{ N/m}$).
   * Equivalent nodal load vector $\mathbf{F}_G$ assembled via local-to-global beam element work equivalence.
   * Maximum static vertical deflection: $\delta_{\max} = -0.1801 \text{ m}$ along the curved track.

5. **Moving Load Dynamic Transient Simulation**:
   * **Athlete Model**: Downward vertical force $P = -90 \text{ kg} \times 9.81 \text{ m/s}^2 = -882.9 \text{ N}$.
   * **Motion Profile**: Initial velocity $v_0 = 2 \text{ m/s}$ at top Point A; constant acceleration $a_1 = 3.5 \text{ m/s}^2$ along linear track; constant acceleration $a_2 = 1.5 \text{ m/s}^2$ along circular curved section.
   * **Transient System Integration**: Modal ODEs $\ddot{q}_i + 2\zeta_i\omega_i \dot{q}_i + \omega_i^2 q_i = Q_i(t)$ solved via explicit Runge-Kutta integration (`ode45`).
   * **Post-Departure Free Vibration Verification**: FFT/peak analysis of free vibration tail matching natural frequencies to 4 decimal places ($f_1 = 1.5654 \text{ Hz}$, $f_2 = 3.0349 \text{ Hz}$, $f_3 = 4.8780 \text{ Hz}$, $f_4 = 5.2966 \text{ Hz}$).

---

## 📁 Module File Structure

```
02_Olympic_Ski_Jump_FEM_Dynamics/
├── main.m                      # Master execution script (Mesh convergence, FEM, FRFs, Moving load)
├── input.inp                   # Struct input file generated dynamically by main.m
└── functions/                  # FEM element transformation & solver routines
    ├── generate_inp_file.m     # Input file generator
    ├── loadstructure.m         # Parsing nodal coordinates & beam element properties
    ├── findcard.m & scom.m     # Sub-parsers for input cards
    ├── el_tra.m                # Element mass & stiffness matrix formulation (6x6)
    ├── assem.m                 # Global M and K matrix assembly
    ├── computeSubMatrices.m    # Partitioning into Free and Constrained DOFs (MFF, KFF, CFF)
    ├── trova_centro.m          # Arch geometry center locator
    ├── obtainConnections.m     # Node discretization & connection generator
    ├── dis_stru.m & diseg2.m   # Undeformed and deformed 2D frame plotting
    ├── displayFRF.m            # Magnitude and phase FRF visualization
    ├── computeForces.m         # Equivalent nodal forces for moving point load
    ├── sortElements.m          # Geometric sorting of track deck elements
    ├── modal_eq.m              # Modal ODE differential system for ode45
    ├── findPeaksAfterForcing.m # Free-vibration frequency & damping extractor
    └── computeDisplaisment.m   # Global displacement peak locator
```

---

## 🚀 How to Run

1. Open MATLAB (R2022b or newer).
2. Set the working directory to `02_Olympic_Ski_Jump_FEM_Dynamics/`.
3. Run the master script:
   ```matlab
   main
   ```
4. The script automatically executes the mesh convergence loop, displays mode shape deformations, plots FRF receptances, evaluates static snow load deflection, simulates the transient response of the accelerating athlete traversing the track, and plots modal coordinate time histories $q(t)$ and tip displacement.

---

## 📄 Detailed Technical Report

For complete 2D FEM beam element matrix derivations, full moving load dynamic transient ODE formulations, static snow load calculations, and extended discussion of results, please consult the complete technical report located in the root `reports/` folder:
* 📄 **Predazzo Olympic Ski Jump Ramp Report**: `reports/Olympic_Ski_Jump_FEM_Report.pdf`


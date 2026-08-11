# Analytical Continuous Beam Vibration Dynamics & Non-Linear Parameter Identification

This module focuses on the continuous analytical modeling of flexural bending vibrations in slender beams based on fourth-order partial differential equations (PDEs), Frequency Response Function (FRF) synthesis, and non-linear least-squares modal parameter estimation (`lsqnonlin`).

---

## 📐 Geometric & Material Parameters

The reference structure is a slender aluminum cantilever beam clamped at $x = 0$ and free at $x = L$.

| Parameter | Symbol | Value | Unit |
|:---|:---:|:---:|:---:|
| **Length** | $L$ | $1.200$ | m |
| **Thickness** | $h$ | $0.008$ ($8\text{ mm}$) | m |
| **Width** | $b$ | $0.040$ ($40\text{ mm}$) | m |
| **Density** | $\rho$ | $2700$ | kg/m$^3$ |
| **Young's Modulus** | $E$ | $68 \times 10^9$ ($68\text{ GPa}$) | Pa |
| **Mass per Unit Length** | $m = \rho b h$ | $0.864$ | kg/m |
| **Second Moment of Area** | $J = \frac{1}{12} b h^3$ | $1.7067 \times 10^{-9}$ | m$^4$ |
| **Bending Stiffness** | $EJ$ | $116.0533$ | N$\cdot$m$^2$ |

---

## 🧮 Mathematical Formulations & Analytical Methodology

1. **Euler-Bernoulli Partial Differential Equation**:
   $$EJ \frac{\partial^4 w(x,t)}{\partial x^4} + m \frac{\partial^2 w(x,t)}{\partial t^2} = 0$$
   Applying a standing-wave ansatz $w(x,t) = \phi(x) g(t)$ yields the spatial eigenvalue equation:
   $$\phi'''' (x) - \gamma^4 \phi(x) = 0, \quad \text{where } \gamma^4 = \frac{\omega^2 m}{EJ}$$
   General spatial solution:
   $$\phi(x) = A \cos(\gamma x) + B \sin(\gamma x) + C \cosh(\gamma x) + D \sinh(\gamma x)$$

2. **Boundary Conditions & Transcendental Frequency Equation**:
   Imposing cantilever boundary conditions ($\phi(0)=0$, $\phi'(0)=0$, $\phi''(L)=0$, $\phi'''(L)=0$) yields a $4 \times 4$ homogeneous linear system $\mathbf{H}(\gamma, \omega) \mathbf{z} = \mathbf{0}$. Non-trivial vibrating solutions exist where the determinant vanishes:
   $$\det [\mathbf{H}(\gamma, \omega)] = 0 \implies \cos(\gamma L) \cosh(\gamma L) + 1 = 0$$

3. **Analytical Natural Frequencies**:
   * **Mode 1**: $f_1 = 4.50 \text{ Hz}$ ($\omega_1 = 28.27 \text{ rad/s}$)
   * **Mode 2**: $f_2 = 28.22 \text{ Hz}$ ($\omega_2 = 177.31 \text{ rad/s}$)
   * **Mode 3**: $f_3 = 79.03 \text{ Hz}$ ($\omega_3 = 496.56 \text{ rad/s}$)
   * **Mode 4**: $f_4 = 154.86 \text{ Hz}$ ($\omega_4 = 973.01 \text{ rad/s}$)

4. **FRF Receptance Synthesis**:
   With structural modal damping ratio $\zeta_n = 1\%$, the receptance between force at $x_k$ and displacement at $x_j$ is synthesized via modal superposition:
   $$G_{jk}(\Omega) = \sum_{n=1}^N \frac{\phi_n(x_j) \phi_n(x_k)}{m_n \left( \omega_n^2 - \Omega^2 + i 2 \zeta_n \omega_n \Omega \right)}$$

5. **Non-Linear Parameter Identification (`lsqnonlin`)**:
   Fitting single-mode complex FRF models in narrow bands around each resonance including inertial ($R_{jk}^L / \Omega^2$) and stiffness ($R_{jk}^H$) residual terms:
   $$H_{jk}^{\text{fit}}(\Omega) \approx \frac{A_{jk}^{(n)}}{\omega_n^2 - \Omega^2 + i 2 \zeta_n \omega_n \Omega} + \frac{R_{jk}^L}{\Omega^2} + R_{jk}^H$$
   The complex residual optimization minimizes:
   $$\varepsilon = \sum_{j} \sum_{i} \left\{ \text{Re}\left[H^{\text{exp}} - H^{\text{fit}}\right]^2 + \text{Im}\left[H^{\text{exp}} - H^{\text{fit}}\right]^2 \right\}$$

6. **Sensor Position Sensitivity Study**:
   * **Optimal Antinodal Setup**: Force at $x_k = 0.2\text{ m}$; outputs at antinodal positions $x_j = \{0.4, 0.7, 1.2\}\text{ m}$. High modal participation, clear resonance peaks.
   * **Nodal Proximity Challenge**: Force at $x_k = 0.2\text{ m}$; outputs near mode shape nodes $x_j = \{0.2, 0.6, 1.1\}\text{ m}$. Output at $0.6\text{ m}$ coincides with a zero-crossing node of Mode 3, attenuating the resonance peak. Optimization robustly identifies modal properties with $\text{MAC} = 1.0000$ across all modes.

---

## 📁 Module File Structure

```
03_Beam_Vibration_Dynamics/
├── main_beam_dynamics.m                # Master analytical & parameter identification script
└── images/                             # Parameter study output figures
    ├── Parametri_1/ ... Parametri_6/   # Sensitivity study plot outputs
```

---

## 🚀 How to Run

1. Open MATLAB (R2022b or newer with Optimization Toolbox).
2. Set the working directory to `03_Beam_Vibration_Dynamics/`.
3. Run the master script:
   ```matlab
   main_beam_dynamics
   ```
4. The script computes exact natural frequencies, plots normalized mode shapes, synthesizes FRFs for both sensor configurations, executes non-linear `lsqnonlin` residual fitting, prints identified frequencies and damping ratios, and evaluates MAC values against analytical mode shapes.

---

## 📄 Detailed Technical Report

For exact PDE standing-wave derivations, transcendental boundary determinant matrix proofs, complex residual formulation details, and extended parameter sensitivity discussion, please consult the complete technical report located in the root `reports/` folder:
* 📄 **Continuous Beam Dynamics Report**: `reports/Beam_Dynamics_Theoretical_Report.pdf`


# Experimental Modal Analysis & FE Model Validation of the Milan Citywave Canopy Structure

This module focuses on the experimental dynamic characterization and finite element model validation of the **Milan Citywave Canopy** roof structure, conducted on a 1:70 scale aeroelastic physical model using 3D Laser Scanning Vibrometry data.

---

## 🏗️ Project Technical Context

The **Citywave Canopy** is an iconic large-span roof structure designed for the CityLife district in Milan. To evaluate its dynamic response under aeroelastic and wind-induced forces prior to full-scale construction, experimental modal testing was performed on a scaled physical model in the Wind Tunnel facility of Politecnico di Milano (*Galleria del Vento GVPM*).

### Scaling & Froude Similitude Criteria
* **Length Scale Factor**: $\lambda_L = 1/70$
* **Frequency Scaling**: Following Froude similitude laws for aeroelastic physical modeling, experimental resonance frequencies ($f_{\text{exp}}$) relate to full-scale prototype frequencies ($f_{\text{full}}$) via:
  $$f_{\text{full}} = f_{\text{exp}} \cdot \sqrt{\lambda_L} \approx 0.11952 \cdot f_{\text{exp}}$$

---

## 🔬 Experimental Methodology & Instrumentation

1. **Excitation & Sensing**:
   * **Excitation**: Electrodynamic shaker applying a swept-sine force signal.
   * **Vibrometry Mesh**: 3D Scanning Laser Doppler Vibrometer (PSV) measuring out-of-plane vertical velocities across a spatial grid of **202 measurement points** distributed over the canopy roof.
   * **Receptance FRFs**: 202 complex Frequency Response Functions ($G_{jk}(\Omega)$) computed using the $H_1$ estimator across the $2 - 8 \text{ Hz}$ frequency band.

2. **Spectral Peak Detection & Windowing**:
   * Averaged absolute magnitude spectrum: $|\bar{G}(\Omega)| = \frac{1}{202} \sum_{j=1}^{202} |G_{jk}(\Omega)|$.
   * Anti-nodal response selection (top 15% amplitude points) to isolate high signal-to-noise ratio channels for each structural mode.

3. **Damping Ratio ($\zeta$) Estimation**:
   * **Half-Power (-3dB) Bandwidth Method**: $\zeta_{HP} = \frac{\Omega_2^2 - \Omega_1^2}{4 \omega_n^2}$.
   * **Phase-Derivative Method**: $\zeta_{PD}$ extracted from the local phase slope at resonance: $\frac{d\theta}{d\omega} = -\frac{2}{\zeta \omega_n}$.
   * **Combined Damping**: Average modal loss factor $\zeta_{\text{used}} = \frac{\zeta_{HP} + \zeta_{PD}}{2}$.

4. **Mode Shape Extraction & FE Spatial Mapping**:
   * Experimental mode shapes extracted from the imaginary part of receptance at resonance: $\mathbf{\Phi}(x_j) \propto \text{Im}[G_{jk}(f_0)]$.
   * Spatial interpolation onto the 3D FE structural node mesh using 2D natural-linear scattered interpolation (`scatteredInterpolant`).
   * Spatial correlation evaluated using the **Modal Assurance Criterion (MAC)**:
     $$\text{MAC}(\mathbf{\phi}_{\text{exp}}, \mathbf{\phi}_{\text{FEM}}) = \frac{\left| \mathbf{\phi}_{\text{exp}}^T \mathbf{\phi}_{\text{FEM}} \right|^2}{\left( \mathbf{\phi}_{\text{exp}}^T \mathbf{\phi}_{\text{exp}} \right) \left( \mathbf{\phi}_{\text{FEM}}^T \mathbf{\phi}_{\text{FEM}} \right)}$$

---

## 📊 Summary of Experimental vs Numerical Results

| Mode | $f_{\text{exp}}$ [Hz] | $f_{\text{full}}$ [Hz] | $f_{\text{FEM}}$ [Hz] | Frequency Error [%] | $\zeta_{HP}$ [%] | $\zeta_{PD}$ [%] | $\zeta_{\text{used}}$ [%] | MAC Value | Deformation Characteristics |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---|
| **1** | 2.9175 | 0.3487 | 0.3390 | **+2.86%** | 1.173% | 1.461% | 1.317% | **0.9360** | Symmetric out-of-plane roof bending |
| **2** | 3.4790 | 0.4158 | 0.4540 | **-8.41%** | 1.077% | 1.734% | 1.406% | **0.6380** | Asymmetric torsion (Modal overlap with Mode 5) |
| **3** | 3.9307 | 0.4698 | 0.4980 | **-5.66%** | 1.009% | 1.755% | 1.382% | **0.6040** | High modal crosstalk & phase distortion |
| **4** | 4.4800 | 0.5355 | 0.5720 | **-6.39%** | 1.913% | 2.282% | 2.097% | **0.6690** | Higher-order roof flexure |

### Key Engineering Insights
* **Fundamental Mode Isolation**: Mode 1 exhibits exceptional correlation ($\text{MAC} = 0.9360$), confirming that the primary dynamic flexibility of the canopy is accurately captured by the numerical FE model.
* **Modal Crosstalk & Mode Pairing**: Modes 2 and 3 occur in a dense spectral region ($3.4 - 4.0 \text{ Hz}$). The 3D 6x6 Cross-MAC matrix identifies strong off-diagonal coupling between Experimental Mode 2 and Numerical Mode 5, explaining the spatial distortion in single-mode quadrature extraction.

---

## 📁 Module File Structure

```
01_Milano_Canopy_Modal_Analysis/
├── main.m                      # Master execution script for EMA processing, MAC & plots
├── data/                       # Experimental & numerical datasets
│   ├── FRF_H1.mat              # 202 experimental complex FRF matrices (G_z)
│   ├── PuntiLaser_FS_FEM.mat   # 3D Laser scanner measurement point coordinates (x, y)
│   ├── connectivity.mat        # FE surface element mesh joint connectivity
│   └── modal_output.mat        # Reference FE model natural frequencies & mode shapes
└── functions/                  # Signal processing & damping algorithms
    ├── selectFRFWindow.m       # Frequency band selection
    ├── selectAntiNodalIdx.m    # Anti-nodal high-SNR point identification
    ├── findPeaks.m             # Spectral peak detection algorithm
    ├── computeHPDamping.m      # Half-Power (-3dB) bandwidth damping estimator
    └── computePDDamping.m      # Phase Derivative slope damping estimator
```

---

## 🚀 How to Run

1. Open MATLAB (R2022b or newer).
2. Set the working directory to `01_Milano_Canopy_Modal_Analysis/`.
3. Run the master script:
   ```matlab
   main
   ```
4. The script automatically loads `data/` and `functions/`, executes peak detection, estimates modal damping, interpolates shapes onto the FE mesh, prints the summary table, and displays side-by-side FEM vs EMA mode shape plots and the 3D Cross-MAC matrix.

---

## 📄 Detailed Technical Report

For full mathematical derivations, complete wind-tunnel experimental setup details, extended MAC matrix discussions, and comprehensive analytical results, please consult the complete technical report located in the root `reports/` folder:
* 📄 **Milan Citywave Canopy Report**: `reports/Milano_Canopy_Modal_Analysis_Report.pdf`


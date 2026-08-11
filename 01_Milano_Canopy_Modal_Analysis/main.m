clear; 
close all; 
clc;

% Add relative functions and data folders to MATLAB path
scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, 'functions'));
addpath(fullfile(scriptDir, 'data'));

%% 0. Load Datasets & System Geometry
load("FRF_H1.mat")            
load("PuntiLaser_FS_FEM.mat") 
load("modal_output.mat")      
load("connectivity.mat")      

f = f(:);          
G = FRF_z;         
nNatF = 4;

N_pts = size(G,2);
fprintf('Measurement points: %d\n', N_pts);
fprintf('Frequency range   : %.3f - %.3f Hz\n', min(f), max(f));

x_laser = x(:);
y_laser = y(:);

[~, Loc1] = ismember(surface_canopy.Joint1, nodes.ID);
[~, Loc2] = ismember(surface_canopy.Joint2, nodes.ID);
[~, Loc3] = ismember(surface_canopy.Joint3, nodes.ID);
nodi123 = [Loc1, Loc2, Loc3]; 
nodi_canopy = unique(nodi123(:));   

% Froude similitude scale factor (1:70 scaled model)
lambda_L = 1/70;             

%% 1. FRF Averaging & Initial Peak Extraction
f_min = 2; f_max = 8;
idx_band = f >= f_min & f <= f_max;
f_band = f(idx_band);
G_band = G(idx_band, :);

G_mean_abs = mean(abs(G_band), 2);
G_mean_abs = G_mean_abs(:);   
f_band = f_band(:);

[nat_f_idx] = findPeaks(f_band, G_mean_abs, 10);

if length(nat_f_idx) > nNatF
    nat_f_idx = nat_f_idx(1:nNatF);
end
nat_f = f_band(nat_f_idx);

for i = 1:nNatF
    fprintf('Mode %d: f_model = %.4f Hz, f_full = %.4f Hz\n', i, nat_f(i), nat_f(i)*sqrt(lambda_L));
end

% Plot experimental FRFs and mean spectrum
figure('Name', 'Experimental FRFs & Mean Spectrum');
subplot(2,1,1)
semilogy(f_band, abs(G_band), 'Color', [0.75 0.75 0.75]); hold on
h1 = semilogy(f_band, G_mean_abs, 'b', 'LineWidth', 2);
h2 = plot(nat_f, G_mean_abs(nat_f_idx), 'ro', 'MarkerSize', 9, 'MarkerFaceColor', 'r');
xlabel('f [Hz]'); 
ylabel('FRFs [m/N]');
grid on; axis tight
legend([h1, h2], 'Mean |FRF|', 'Identified Peaks', 'Location', 'northeast');

subplot(2,1,2)
plot(f_band, angle(G_band), 'Color', [0.75 0.75 0.75]); hold on
xlabel('f [Hz]'); 
ylabel('Phase [rad]');
title('Phase of 202 Measurement Points')
grid on; axis tight

%% 2. Modal Parameter Identification
results = struct();   
antiNodalPosPerc = [0.85 0.95 0.7 0.8];
PDdampingInterPoint = [2 3 3 4];

for ii = 1:nNatF
    fprintf('\n---------------------------------------------\n');
    fprintf('MODE %d IDENTIFICATION\n', ii);
    fprintf('---------------------------------------------\n');
    
    % Select frequency window & anti-nodal response points
    [f_win, G_win] = selectFRFWindow(nat_f(ii), 0.2, f, G);
    antinodal_idx = selectAntiNodalIdx(f_win, G_win, temp_nat_f, antiNodalPosPerc(ii), N_pts);

    % Refine resonance frequency f0 on averaged anti-nodal FRF
    G_avg = mean(G_win(:, antinodal_idx), 2);   
    G_abs_avg = abs(G_avg);

    nat_f_idx = findPeaks(f_win, G_abs_avg, 10);
    nat_f_idx = nat_f_idx(1);

    G_avg_max = G_abs_avg(nat_f_idx);
    f0 = f_win(nat_f_idx); 
    omega0 = 2*pi*f0;
    fprintf('Refined natural frequency: f0 = %.4f Hz (Full-scale: %.4f Hz)\n', f0, f0*sqrt(lambda_L));

    % Half-Power & Phase-Derivative Damping Estimation
    [h_hp, f1_hp, f2_hp, HP_level] = computeHPDamping(nat_f_idx, G_avg_max, f_win, G_abs_avg);
    [h_pd, phase_avg, i_lo, i_hi, p] = computePDDamping(nat_f_idx, f_win, G_avg, PDdampingInterPoint(ii));
    h_used = mean([h_hp, h_pd], 'omitnan');

    % Extract mode shape from quadrature (imaginary part of FRFs at f0)
    [~, idx_f0_global] = min(abs(f - f0));
    phi_exp_raw = imag(G(idx_f0_global, :)).';   
    [~, idx_max_abs] = max(abs(phi_exp_raw));
    phi_exp = phi_exp_raw / phi_exp_raw(idx_max_abs);

    % Interpolate laser measurements onto FEM mesh grid
    F_interp = scatteredInterpolant(x_laser, y_laser, phi_exp, 'natural', 'linear');
    phi_exp_FEM_raw = F_interp(nodes.X, nodes.Y);
    
    valid_mask = ~isnan(phi_exp_FEM_raw);
    phi_exp_FEM = phi_exp_FEM_raw;
    phi_exp_FEM(~valid_mask) = 0;    

    results(ii).mode = ii;
    results(ii).f0_model = f0;
    results(ii).f0_full = f0 * sqrt(lambda_L);       
    results(ii).h_hp = h_hp;
    results(ii).h_pd = h_pd;
    results(ii).h_used = h_used;
    results(ii).phi_exp = phi_exp;
    results(ii).phi_exp_FEM = phi_exp_FEM;
    results(ii).valid_mask = valid_mask;

    % Mode identification plot
    figure('Name', sprintf('Mode %d - Parameter Identification', ii));
    subplot(2,1,1)
    semilogy(f_win, abs(G_win(:, antinodal_idx)), 'Color', [0.8 0.8 0.8]); hold on
    semilogy(f_win, abs(G_avg), 'b', 'LineWidth', 2);
    plot(f0, G_avg_max, 'ro', 'MarkerSize', 9, 'MarkerFaceColor', 'r');
    if ~isnan(f1_hp) && ~isnan(f2_hp)
        yline(HP_level, 'g:', 'HP Level');
        plot([f1_hp f2_hp], [HP_level HP_level], 'go-', 'LineWidth', 1.5, 'MarkerFaceColor', 'g', 'MarkerSize', 7);
    end
    xlabel('f [Hz]'); ylabel('|FRF| [m/N]');
    title(sprintf('Mode %d: f_0 = %.4f Hz | h_{HP} = %.3f%% | h_{PD} = %.3f%%', ii, f0, h_hp*100, h_pd*100))
    grid on; axis tight

    subplot(2,1,2)
    plot(f_win, phase_avg, 'b', 'LineWidth', 2); hold on
    omega_loc = 2*pi*f_win(i_lo:i_hi);
    plot(f_win(i_lo:i_hi), polyval(p, omega_loc), 'r--', 'LineWidth', 1.7);
    xline(f0, 'k:', 'f_0');
    xlabel('f [Hz]'); ylabel('Phase [rad]');
    title('Mean FRF Phase & Linear Tangent (Phase Derivative)')
    grid on; axis tight
    legend('Mean Phase', 'Linear Tangent', 'Location', 'best')
end

%% 3. FEM Comparison, MAC Matrix & Results Summary
scalaFEM = 10;
scalaLASER = 10;

fprintf('\n\n========================================================================================\n');
fprintf('SUMMARY RESULTS TABLE\n');
fprintf('========================================================================================\n');
fprintf('%-5s %-12s %-12s %-12s %-10s %-12s %-12s %-12s %-8s\n', 'Mode', 'f_exp_mod', 'f_exp_full', 'f_FEM_full', 'err [%]', 'h_HP [%]', 'h_PD [%]', 'h_USED [%]', 'MAC');
fprintf('----------------------------------------------------------------------------------------\n');

PHI_EXP_MATRIX = zeros(202, nNatF); 

for ii = 1:nNatF
    % Extract numerical mode shape from FEM dataset
    mode_sel = modeshapes(modeshapes.No == ii, :);
    [~, Locb] = ismember(mode_sel.ID, nodes.ID);
    modedef = mode_sel{Locb, {'uX', 'uY', 'uZ'}};

    [~, idx_max_FEM] = max(abs(modedef(:,3)));
    sgn_FEM = sign(modedef(idx_max_FEM,3));
    modedef = modedef / (sgn_FEM * abs(modedef(idx_max_FEM,3)));
    phi_FEM_norm = modedef(:,3);

    % Map FEM nodes to laser measurement grid
    idx_canopy = ismember((1:numel(nodes.X))', nodi_canopy);
    F_fem = scatteredInterpolant(nodes.X(idx_canopy), nodes.Y(idx_canopy), phi_FEM_norm(idx_canopy), 'natural', 'linear');
    phi_FEM_at_laser = F_fem(x_laser, y_laser);
    phi_exp_laser = results(ii).phi_exp;
    
    PHI_EXP_MATRIX(:, ii) = phi_exp_laser;

    % Phase alignment
    valid_mac = ~isnan(phi_FEM_at_laser) & ~isnan(phi_exp_laser);
    vec_FEM = phi_FEM_at_laser(valid_mac);
    vec_EXP = phi_exp_laser(valid_mac);
    if dot(vec_EXP, vec_FEM) < 0
        vec_EXP = -vec_EXP;
        results(ii).phi_exp = -results(ii).phi_exp;         
        results(ii).phi_exp_FEM = -results(ii).phi_exp_FEM; 
    end

    % Modal Assurance Criterion (MAC) calculation
    MAC = (vec_EXP.' * vec_FEM)^2 / ((vec_EXP.' * vec_EXP) * (vec_FEM.' * vec_FEM));
    results(ii).MAC = MAC;

    f_FEM_full = modpar.freq(ii);
    err_pct = (results(ii).f0_full - f_FEM_full) / f_FEM_full * 100;

    fprintf('%-5d %-12.4f %-12.4f %-12.4f %-10.2f %-12.3f %-12.3f %-12.3f %-8.4f\n', ...
            ii, results(ii).f0_model, results(ii).f0_full, f_FEM_full, err_pct, ...
            results(ii).h_hp*100, results(ii).h_pd*100, results(ii).h_used*100, MAC);

    % Side-by-side mode shape plots (FEM vs EMA)
    figure('Name', sprintf('Mode %d Comparison: FEM vs EMA', ii));
    subplot(1,2,1)
    cc = phi_FEM_norm;
    patch('Faces', nodi123, 'Vertices', [nodes.X, nodes.Y, nodes.Z] + modedef*scalaFEM, ...
          'CData', cc, 'FaceColor', 'interp', 'EdgeColor', 'none')
    axis equal; axis tight; grid on; view(2)
    colormap(gca, jet(20)); 
    clim([-1 1]); 
    colorbar('eastoutside')
    title(sprintf('FEM - Mode %d: f = %.4f Hz', ii, f_FEM_full))
   
    subplot(1,2,2)
    cc_exp = results(ii).phi_exp_FEM;
    patch('Faces', nodi123, 'Vertices', [nodes.X, nodes.Y, nodes.Z + cc_exp*scalaLASER], ...
          'CData', cc_exp, 'FaceColor', 'interp', 'EdgeColor', 'none')
    hold on
    plot3(x_laser, y_laser, ones(size(x_laser))*1.5, 'k.', 'MarkerSize', 4)
    axis equal; axis tight; grid on; view(2)
    colormap(gca, jet(20)); 
    clim([-1 1]); 
    colorbar('eastoutside')
    title(sprintf('EMA - Mode %d: f = %.4f Hz | h = %.3f%% | MAC = %.3f', ...
            ii, results(ii).f0_full, results(ii).h_used*100, MAC))
end
fprintf('========================================================================================\n');

%% 4. Cross-MAC Matrix & Mode Pairing
n_fem_modes = nNatF + 2; 
PHI_FEM_FULL = zeros(202, n_fem_modes);

idx_canopy = ismember((1:numel(nodes.X))', nodi_canopy);
x_canopy = nodes.X(idx_canopy);
y_canopy = nodes.Y(idx_canopy);

for jj = 1:n_fem_modes
    m_sel = modeshapes(modeshapes.No == jj, :);
    [~, Loc] = ismember(m_sel.ID, nodes.ID);
    m_def = m_sel{Loc, {'uX', 'uY', 'uZ'}};
    
    [~, idxmax] = max(abs(m_def(:,3)));
    m_def = m_def / (sign(m_def(idxmax,3)) * abs(m_def(idxmax,3)));
    
    F_f = scatteredInterpolant(x_canopy, y_canopy, m_def(idx_canopy,3), 'natural', 'boundary');
    phi_f_laser = F_f(x_laser, y_laser);
    phi_f_laser(isnan(phi_f_laser)) = 0; 
    
    PHI_FEM_FULL(:, jj) = phi_f_laser;
end

CrossMAC = zeros(nNatF, n_fem_modes);
for i = 1:nNatF
    v_exp = PHI_EXP_MATRIX(:, i);
    v_exp(isnan(v_exp)) = 0; 
    for jj = 1:n_fem_modes
        v_fem = PHI_FEM_FULL(:, jj);
        CrossMAC(i,jj) = (v_exp.' * v_fem)^2 / ((v_exp.' * v_exp) * (v_fem.' * v_fem));
    end
end

figure('Name', 'Cross-MAC Matrix');
b = bar3(CrossMAC);
for k = 1:length(b)
    zdata = b(k).ZData;
    b(k).CData = zdata;
    b(k).FaceColor = 'interp';
end
colormap(jet); colorbar; clim([0 1]);
xlabel('Numerical Modes (FEM)');
ylabel('Experimental Modes (EMA)');
zlabel('MAC');
set(gca, 'XTick', 1:n_fem_modes, 'YTick', 1:nNatF);
view(-45, 45);
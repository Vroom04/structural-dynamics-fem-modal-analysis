clc; 
clear; 
close all;

% Add relative functions folder to MATLAB path
scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, 'functions'));

%% 1. Structural Geometry & Beam Section Properties
props = [
    1, 122, 3.21e9, 2.0e8;   % Red beams
    2, 80,  1.74e9, 4.8e7;   % Green beams
    3, 400, 2.5e9,  4.0e9    % Blue beams
];

% Column ground attachment geometry
a = 15 - (17 * 42.5 / 67.6);
theta = 42.5 / 67.6;

% Curved arch geometry calculation
P1 = [128, 0];      
P2 = [67.6, 17.5];  
R  = 112;           
[C1, C2] = trova_centro(P1, P2, R);
   
% Key FEM nodal coordinates
f_nodes = [
    % Pillar 1 (Left framework)
    1, 1, 1, 1,  0.0,   0.0;    
    2, 0, 0, 0,  0.0,   15.0;   
    3, 0, 0, 0,  0.0,   30.0;
    4, 0, 0, 0,  0.0,   45.0;   
    5, 0, 0, 0,  0.0,   60.0;   
    6, 1, 1, 1,  17.0,  0.0;    
    7, 0, 0, 0,  17.0,  a;   
    8, 0, 0, 0,  17.0,  a+15;
    9, 0, 0, 0,  17.0,  a+30;   
    10, 0, 0, 0, 17.0,  a+45;  
        
    % Pillar 2 (Center framework)
    11, 1, 1, 1, 50.60, 0.0;    
    12, 0, 0, 0, 50.60, a;   
    13, 0, 0, 0, 50.60, 17.5-a; 
    14, 0, 0, 0, 50.60, 17.5-a+15; 
    15, 1, 1, 1, 67.60, 0.0;    
    16, 0, 0, 0, 67.60, 8.75;   
    17, 0, 0, 0, 67.60, 17.5;   
        
    % Ramp Terminal Point B
    18, 0, 0, 0, 128.0, 0.0  
];

iso_nodes = [];

% Maximum element length bounds based on beam properties
sf = 1.5;
fmax = 6;
nmods = 4;
omegamax = fmax * 2 * pi;
Lrmax = sqrt(pi^2 / sf / omegamax * sqrt(props(1,4)/props(1,2)));
Lgmax = sqrt(pi^2 / sf / omegamax * sqrt(props(2,4)/props(2,2)));
Lbmax = sqrt(pi^2 / sf / omegamax * sqrt(props(3,4)/props(3,2)));

%% 2. Mesh Convergence Study
freqm = [];
nnods = [];
startValue = 0.1;
step = 0.1;
stopValue = 5;

for ii = startValue:step:stopValue
    Lr = Lrmax / ii;
    Lg = Lgmax / ii;
    Lb = Lbmax / ii;
    
    conns = obtainConnections(Lr, Lg, Lb, C2);
    bluebeams = generate_inp_file('input.inp', f_nodes, conns, iso_nodes, props, 3);
    
    [file_i, xy, nnod, sizee, idb, ndof, incid, l, gamma, m, EA, EJ, posiz, nbeam, pr] = loadstructure;
    [M, K] = assem(incid, l, m, EA, EJ, gamma, idb);
    
    % Boundary spring element at Point B (Node 18)
    i_dof18 = idb(18,:);
    k = 5e7;
    K(i_dof18(2), i_dof18(2)) = K(i_dof18(2), i_dof18(2)) + k;
    
    [MFF, MCF, MFC, MCC] = computeSubMatrices(ndof, M);
    [KFF, KCF, KFC, KCC] = computeSubMatrices(ndof, K);
    
    [modes, omega2] = eig(MFF\KFF);
    omega = sqrt(diag(omega2));
    [omega, i_omega] = sort(omega);
    freq0 = omega / 2 / pi;
    modes = modes(:, i_omega);
    
    freq0 = freq0(1:nmods);
    modes = modes(:, 1:nmods);        
    freqm(:, end+1) = freq0;
    nnods(end+1) = nnod;
end

% Plot Mesh Convergence
figure('Color', 'w', 'Name', 'Mesh Convergence Study');
hold on; grid on; box on;
colori = ['r', 'b', 'g', 'm', 'k'];
LmaxIdx = (1 - startValue) / step + 1;

for modo = 1:nmods
    plot(nnods, freqm(modo, :), [colori(modo) '-o'], ...
         'LineWidth', 2, 'MarkerFaceColor', colori(modo), ...
         'DisplayName', sprintf('Mode %d', modo));
    plot(nnods(LmaxIdx), freqm(modo, LmaxIdx), 'o', ...
         'MarkerEdgeColor', colori(5), 'MarkerFaceColor', colori(5), ... 
         'MarkerSize', 8, 'HandleVisibility', 'off');      
end

xlabel('Total Number of FEM Nodes');
ylabel('Natural Frequencies [Hz]');
title('Natural Frequencies Convergence Study');
legend('Location', 'best');

%% 3. Final Mesh FEM Assembly & Eigenmode Extraction
nIterStandard = floor((1 - startValue) / step + 1);
nIter = 30;

Lr = Lrmax / ((nIter - 1) * step + startValue);
Lg = Lgmax / ((nIter - 1) * step + startValue);
Lb = Lbmax / ((nIter - 1) * step + startValue);

conns = obtainConnections(Lr, Lg, Lb, C2);
bluebeams = generate_inp_file('input.inp', f_nodes, conns, iso_nodes, props, 3);

[file_i, xy, nnod, sizee, idb, ndof, incid, l, gamma, m, EA, EJ, posiz, nbeam, pr] = loadstructure;
[M, K] = assem(incid, l, m, EA, EJ, gamma, idb);

i_dof18 = idb(18,:);
k = 5e7;
K(i_dof18(2), i_dof18(2)) = K(i_dof18(2), i_dof18(2)) + k;

[MFF, MCF, MFC, MCC] = computeSubMatrices(ndof, M);
[KFF, KCF, KFC, KCC] = computeSubMatrices(ndof, K);

[modes, omega2] = eig(MFF\KFF);
omega = sqrt(diag(omega2));
[omega, i_omega] = sort(omega);
freq0 = omega / 2 / pi;
modes = modes(:, i_omega);
freq0 = freq0(1:nmods);
modes = modes(:, 1:nmods);        

% Mode shape normalization
for ii = 1:nmods
    tempMax = abs(modes(1, ii));
    for jj = 1:size(modes,1)
        if abs(modes(jj, ii)) > abs(tempMax)
            tempMax = modes(jj,ii);
        end
    end
    modes(:, ii) = modes(:, ii) ./ (-tempMax);
end

% Plot undeformed & deformed mode shapes
scale_factor = 3;
dis_stru(posiz, l, gamma, xy, pr, idb, ndof, 'Undeformed Structure');

for ii = 1:nmods
    figure();
    grid on; grid minor;
    title(sprintf('Mode Shape %d | Frequency: %.2f Hz', ii, freq0(ii)));
    xlabel('X [m]'); ylabel('Y [m]');
    axis equal; 
    diseg2(modes(:,ii), scale_factor, incid, l, gamma, posiz, idb, xy);
    box on; hold off;
end

%% 4. Rayleigh Damping & Harmonic FRF Analysis
ab = [0.1; 2.0e-4];
C = ab(1)*M + ab(2)*K;
CFF = C(1:ndof, 1:ndof);

% Harmonic horizontal force at Node A
F0 = zeros(ndof, 1);
F0(idb(5,1)) = 1;
om = 0:0.01*2*pi:omegamax;

X = zeros(ndof, length(om));
for ii = 1:length(om)
    A = -om(ii)^2*MFF + 1i*om(ii)*CFF + KFF;
    X(:,ii) = A \ F0;
end

displayFRF(idb(5,1), om, X, 'Horizontal Receptance at Node A (Force at A)', 1);
displayFRF(idb(18,2), om, X, 'Vertical Receptance at Point B (Force at A)', 1);

%% 5. Reduced-Order Modeling (Modal Superposition)
F0 = zeros(ndof, 1);
F0(idb(18,2)) = 1;

for ii = 1:length(om)
    A = -om(ii)^2*MFF + 1i*om(ii)*CFF + KFF;
    X(:,ii) = A \ F0;
end

% 2-Mode reduced superposition
ii = 1:2;
Phi = modes(:,ii);
Mmod = Phi'*MFF*Phi; 
Kmod = Phi'*KFF*Phi; 
Cmod = Phi'*CFF*Phi; 
Fmod = Phi'*F0;

xx_mod = zeros(2, length(om));
for ii = 1:length(om)
    xx_mod(:,ii) = (-om(ii)^2*Mmod + 1i*om(ii)*Cmod + Kmod) \ Fmod; 
end
X_m = Phi * xx_mod; 

displayFRF(idb(5,1), om, X, 'Full FEM vs 2-Mode Superposition Receptance', 1);
displayFRF(idb(5,1), om, X_m, 'Full FEM vs 2-Mode Superposition Receptance', 0);
legend('Full FE Model', 'Modal Superposition (2 modes)', 'Location', 'best')

%% 6. Static Analysis Under Snow Distributed Load
p0 = -7000; 
p0G = [0 p0]'; 
FG = zeros(3*nnod, 1); 

for ii = 1:length(bluebeams)
    gammaii = gamma(bluebeams(ii)); 
    Lk = l(bluebeams(ii)); 
    p0L = [cos(gammaii) sin(gammaii); -sin(gammaii) cos(gammaii)] * p0G; 
    
    FkL = [Lk/2; 0; 0; Lk/2; 0; 0]*p0L(1) + [0; Lk/2; Lk^2/12; 0; Lk/2; -Lk^2/12]*p0L(2);
    
    FkG = [cos(gammaii) -sin(gammaii) 0 0 0 0;
           sin(gammaii) cos(gammaii) 0 0 0 0; 
           0 0 1 0 0 0; 
           0 0 0 cos(gammaii) -sin(gammaii) 0; 
           0 0 0 sin(gammaii) cos(gammaii) 0; 
           0 0 0 0 0 1] * FkL;

    Ek = zeros(6, 3*nnod); 
    for jj = 1:6
        Ek(jj, incid(bluebeams(ii),jj)) = 1;
    end 
    FG = FG + Ek'*FkG;
end

F0 = FG(1:ndof, 1);
xf = squeeze(KFF \ F0);

figure('Name', 'Deformed Shape Under Snow Load', 'Color', 'w'); 
hold on;
diseg2(xf, scale_factor, incid, l, gamma, posiz, idb, xy);
grid on; grid minor;
title('Deformed Shape Under Snow Load');
xlabel('X [m]'); ylabel('Y [m]');
axis equal; box on; hold off;

[max_pvd, idx_pd, max_nvd, idx_nd] = computeDisplaisment(idb, nnod, ndof, xf);
fprintf('Max positive vertical displacement: %.4f m at Node %d\n', max_pvd, idx_pd);
fprintf('Max negative vertical displacement: %.4f m at Node %d\n', max_nvd, idx_nd);

%% 7. Transient Dynamic Response Under Moving Skier Load
P = -90 * 9.81; % 90 kg athlete weight

vo = 2;
a1 = 3.5;
a2 = 1.5;
pos = 0;

linearSegmentLength = sqrt(67.6^2 + 42.5^2);
firstSegmentLength = sqrt(17^2 + (15-a)^2);
secondSegmentLength = linearSegmentLength - 2*firstSegmentLength;

line_time = (-vo + sqrt(vo^2 + 2*a1*linearSegmentLength)) / a1;
firstSegTime = (-vo + sqrt(vo^2 + 2*a1*firstSegmentLength)) / a1;
secondSegTime = (-vo + sqrt(vo^2 + 2*a1*(firstSegmentLength + secondSegmentLength))) / a1;

total_time = 0;
dt = 0.0001;
t = 0:dt:50;

bluebeams = sortElements(bluebeams, 1, posiz);
Mmod = modes'*MFF*modes; 
Kmod = modes'*KFF*modes; 
Cmod = modes'*CFF*modes; 
Q_t = zeros(size(modes,2), length(t));

for ii = 1:length(t)
    v_end_line = vo + a1 * line_time; 
    pos_end_line = vo * line_time + 0.5 * a1 * line_time^2;

    if t(ii) < line_time
        pos = vo * t(ii) + 0.5 * a1 * t(ii)^2; 
    else
        pos = pos_end_line + v_end_line * (t(ii) - line_time) + 0.5 * a2 * (t(ii) - line_time)^2;
    end

    n_element = -1;
    space = 0;
    for jj = 1:length(bluebeams)
        space = space + l(bluebeams(jj));
        if space > pos
            n_element = bluebeams(jj);
            break;
        end
    end
   
    if n_element == -1
        if total_time == 0
            total_time = t(ii);
        end
        continue;
    end

    l2 = space - pos;
    l1 = l(n_element) - l2;
    Lk = l(n_element);
    gammaii = gamma(n_element); 
  
    F0 = computeForces(P, l1, l2, Lk, incid, n_element, ndof);
    Q_t(:, ii) = modes' * F0;
end

% Solve modal ODE system using ode45
q_time = zeros(length(t), nmods);
for i = 1:nmods
    m_ii = Mmod(i, i);
    k_ii = Kmod(i, i);
    c_ii = Cmod(i, i);    
    [~, Y_out] = ode45(@(t_ode, q) modal_eq(t_ode, q, m_ii, c_ii, k_ii, t, Q_t(i,:)), t, [0; 0]);
    q_time(:, i) = Y_out(:, 1);
end

% Extract free-vibration frequencies post-departure
max_pts = zeros(nmods, 2);
index = zeros(nmods, 2);
for ii = 1:nmods
    [max_pts(ii, :), index(ii, :)] = findPeaksAfterForcing(t, q_time(:,ii), total_time);
end

% Plot modal coordinates time history
figure('Name', 'Time History of Modal Coordinates q(t)'); 
hold on; grid on; grid minor; 
legend_labels = cell(size(modes, 2), 1);
for i = 1:nmods
    plot(t, q_time(:, i), 'LineWidth', 1.5); 
    legend_labels{i} = sprintf('q%d', i);
end

for i = 1:nmods
    plot(t(index(i, :)), q_time(index(i,:), i), 'o', 'Color', 'b');
    fprintf('Lagrangian modal coordinate frequency (Mode %d): %.4f Hz\n', i, 1/(t(index(i,2))-t(index(i,1))));
end

xline(firstSegTime, '--r', 'Label', 'End of Pillar 1', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
xline(secondSegTime, '--r', 'Label', 'Start of Pillar 2', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
xline(line_time, '--r', 'Label', 'End of Linear Ramp', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
xline(total_time, '--r', 'Label', 'End of Curved Track', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
title('Time History of Modal Coordinates q(t)');
xlabel('Time [s]'); ylabel('q');
legend(legend_labels, 'Location', 'best', 'FontSize', 11);
box on; hold off;

% Plot tip (Point B) vertical displacement time history
gdl_18_v = idb(18,2); 
disp_node_18 = zeros(length(t), 1);
for i = 1:size(modes,2)
    disp_node_18 = disp_node_18 + modes(gdl_18_v, i) * q_time(:, i);
end

figure('Name', 'Tip Vertical Displacement (Point B)');
plot(t, disp_node_18, 'LineWidth', 1.5);
xline(firstSegTime, '--r', 'Label', 'End of Pillar 1', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
xline(secondSegTime, '--r', 'Label', 'Start of Pillar 2', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
xline(line_time, '--r', 'Label', 'End of Linear Ramp', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
xline(total_time, '--r', 'Label', 'End of Curved Track', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
title('Vertical Displacement of Ramp Tip (Point B) under Moving Load');
xlabel('Time [s]'); ylabel('Displacement [m]');
grid on;
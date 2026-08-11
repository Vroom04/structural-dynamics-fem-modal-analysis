%%
clear; 
close all; 
clc;

% Add relative images folder to MATLAB path
scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, 'images'));

%% Input geometrical data

l=1.2; %m length
h=0.008; %m tickness     
b=0.040; %m width
p=2700; %kg/m^3 density
E=68e9; %Pa young's modulus
m=b*h*p; %kg/m mass per unit length
I=1/12*b*h^3; %inertia modulus
EI=E*I;

%% Set frequency range and frequency resolution

df = 1.0000e-02;         % Hz
fmax = 200;              %Hz 

%% Set space domain and resolution

dx = 0.001;               % m
x = 0:dx:l;

%% frequency range

F=0:df:fmax;
omega=F.*2.*pi;
g=(omega.^2.*m./EI).^(1/4); % gamma beam

%% BCs in matrix form

H=@(g,omega)  [ +1                    0                   +1                   0               ;
                0                    +g                    0                  +g               ;
                -g^2*cos(g*l)        -g^2*sin(g*l)        +g^2*cosh(g*l)      +g^2*sinh(g*l)   ;
                +EI*g^3*sin(g*l)     -EI*g^3*cos(g*l)     +EI*g^3*sinh(g*l)   +EI*g^3*cosh(g*l);];

%% eigenfrequencies and eigenvectors

%computes determinant at each frequency
dets=[];
for i=1:length(omega)
    dets(end+1)=det(H(g(i),omega(i)));
end

%plots the det absolute value
figure(1), box on
semilogy(F,abs(dets),'-b')
hold on, grid on, xlabel('f [Hz]')

%finds the natural frequencies indeces
i_nat=[];
for i=2:length(dets)-1
    if abs(dets(i)) < abs(dets(i-1)) && abs(dets(i)) < abs(dets(i+1))
        i_nat(end+1)=i;
    end
end

%plot natural frequencies
plot(F(i_nat),abs(dets(i_nat)),'or')

%computes eigenvectors
C_hat=zeros(4,length(i_nat));
for i=1:length(i_nat)  
    omega_i=omega(i_nat(i)); %omega i-th
    g_i=g(i_nat(i)); %gamma i-th

    Hi=H(g_i,omega_i); %compute H in omega i-th
    Hi_hat=Hi(1:3,1:3); %reduced square matrix
    Ei_hat=Hi(1:3,4); %reduced column vector

    Ci_hat=[ -Hi_hat\Ei_hat;1]; %eigenmode i-th
    C_hat(:,i)=Ci_hat; %matrix of eigenvectors
end

%% Modeshape functions

for i=1:length(i_nat)
    omega_i=omega(i_nat(i));
    g_i=g(i_nat(i));
    phi(i,:) = C_hat(1,i)*cos(g_i*x)  + ...
        C_hat(2,i)*sin(g_i*x) + ...
        C_hat(3,i)*cosh(g_i*x) + ...
        C_hat(4,i)*sinh(g_i*x);
end

%Mode shapes normalization
for i=1:length(i_nat)
    maxvalue = max(abs(phi(i,:)));
    phi(i,:) = phi(i,:) / maxvalue;
end

%plot modeshapes
figure; hold on; grid on
n_modes = size(phi,1);
colors = lines(n_modes);
for i = 1:n_modes
    plot(x, phi(i,:), 'LineWidth',3, 'Color', colors(i,:), 'DisplayName', sprintf('mode %d f_0 = %.1f Hz', i, F(i_nat(i))));
end
legend('Location', 'northwest');
title('Normalized mode shapes');
xlabel('Distance x from the end [m]');
ylabel(sprintf('\\phi (norm.)'));

%% Second part

natfs = omega(i_nat);

%% Sensor positions

x_k = 0.2;
x_j = [0.2 0.6 1.1];

%% Modal Mass and Damping

n_modi = size(phi, 1);
M = zeros(n_modi, n_modi);

csi = 0.01;

for i = 1:n_modi
    for j = 1:n_modi
        M(i,j) = trapz(x, m * phi(i,:) .* phi(j,:));
    end
end

%% Compute FRFs

FRFs = zeros(length(x_k), length(x_j), length(omega));

for k = 1:length(x_k)
    for j = 1:length(x_j)
        
        %spatial indeces
        ik = round(x_k(k)/dx) + 1;
        ij = round(x_j(j)/dx) + 1; 
        
        temp_G = zeros(size(omega));
        
        for i = 1:length(i_nat)
            numerator = (phi(i, ij) * phi(i, ik)) / M(i,i);
            denominator = (natfs(i)^2 - omega.^2 + 2*1i*csi*natfs(i)*omega);
            temp_G = temp_G + (-numerator./ denominator);
        end
        FRFs(k, j, :) = temp_G;
    end
end

%% Display FRFs

n_inputs = length(x_k);
n_outputs = length(x_j);

for k = 1:n_inputs
    figure('Name', sprintf('Input at %.1f m', x_k(k)), 'NumberTitle', 'off');
    
    for j = 1:n_outputs
        subplot(2, n_outputs, j);
        semilogy(F, abs(squeeze(FRFs(k, j, :))), 'LineWidth', 1.2);
        grid on;
        xlim([0 fmax]);
        title(sprintf('Output: %.1f m', x_j(j)));
        xlabel('f [Hz]');
        ylabel('Mag [m/N]');
        
        subplot(2, n_outputs, j + n_outputs);
        plot(F, angle(squeeze(FRFs(k, j, :))), 'LineWidth', 1.2);
        grid on;
        xlim([0 fmax]);
        xlabel('f [Hz]');
        ylabel('Phase [rad]');
        ylim([-5 5]);
    end
end

%% why we can use lsnonlin
band_lo = 0.85;    
band_hi = 1.15;

for j = 1:n_inputs
    for i=1:n_outputs
%        figure();
        for ii=1:length(i_nat)
            %calculate the frequency range Ws
            f_lo = band_lo * F(i_nat(ii));
            f_hi = band_hi * F(i_nat(ii));
            idx = (F >= f_lo) & (F <= f_hi);
            cut_f = F(idx);
            cut_G = FRFs(j,i,idx);
            cut_G = cut_G(:);
%            subplot(2, 2, ii);
%            plot(real(cut_G), imag(cut_G));
        end
    end
end


%% Curve fitting

%vector and matrices of parameters that will be identified base on the
%approximation

W = F*2*pi;

band_lo = 0.95;    
band_hi = 1.05;

omega_guess = [];
omega_id=[];
zeta_id = [];
A_id=[];
RL_id=[];
RH_id=[];

%guess on resonance frequencies and damping ratio
G_mean = squeeze(mean(abs(FRFs), [1,2]));        
for i=2:length(W)-1            
    if G_mean(i) > G_mean(i-1) && G_mean(i) > G_mean(i+1)            
        %guess on frequency
        omega_guess(end+1)=W(i);

        %guess on damping ratio

        %bounds the spectrum 
        w_lo = band_lo * omega_guess(end);
        w_hi = band_hi * omega_guess(end);
        idx = (W >= w_lo) & (W <= w_hi);
        cut_W = W(idx);

        %isolates FRFs in an interval around the pick of the current
        %resonance frequency
        cut_FRFs = FRFs(:, :, idx);
            
        %this part finds the best FRF for computing the damping
        [max_pick, index] = max(abs(cut_FRFs), [], 'all');
        [r, c, k] = ind2sub(size(cut_FRFs), index);
        best_FRF_complex = squeeze(cut_FRFs(r, c, :)); 
        best_FRF_abs = abs(best_FRF_complex);
        
        %prepares data for half-power method
        half_pwr = max_pick/sqrt(2);
        i_lo = find(best_FRF_abs(1:k) <= half_pwr, 1, 'last');
        i_hi_rel = find(best_FRF_abs(k:end) <= half_pwr, 1, 'first');
        if ~isempty(i_lo) && ~isempty(i_hi_rel)
            i_hi = k + i_hi_rel - 1; 
            zeta_id(end+1) = (cut_W(i_hi)^2 - cut_W(i_lo)^2)/(4*cut_W(k)^2);
        else
            %if the half power points are not inside the interval it uses
            %the phase 
            fase = squeeze(unwrap(angle(cut_FRFs(r, c, :)))); 
            pendenza_fase = gradient(fase, cut_W); 
            slope_at_peak = pendenza_fase(k); 
            zeta_id(end+1) = abs(1 / (cut_W(k) * slope_at_peak));
        end

        %guess on the amplitude
        m = length(omega_guess); 
        for r = 1:n_inputs
            for c = 1:n_outputs
                [~, k_loc] = max(abs(squeeze(cut_FRFs(r, c, :))));
                Gpk = cut_FRFs(r, c, k_loc); 
                A_id(r, c, m) = -2 * zeta_id(end) * omega_guess(end)^2 * imag(Gpk);
            end
        end

        %residuals set to zero
        RL0 = zeros(n_inputs, n_outputs, 1);   
        RH0 = zeros(n_inputs, n_outputs, 1);
        
        x0 = [ omega_guess(m); ...
        zeta_id(m); ...
        reshape(real(A_id(:,:,m)), [], 1); ...
        reshape(imag(A_id(:,:,m)), [], 1); ...
        reshape(real(RL0), [], 1); ...
        reshape(imag(RL0), [], 1); ...
        reshape(real(RH0), [], 1); ...
        reshape(imag(RH0), [], 1) ];

        x_opt = lsqnonlin(@(xx) errFRF(xx, cut_W, cut_FRFs, n_inputs, n_outputs), x0, [],[],[]);
            
        %parameter estraction
        P = n_inputs * n_outputs;
       
        omega_id(m) = x_opt(1);
        zeta_id(m)  = x_opt(2);    
       
        A_real = reshape(x_opt(3 : 2+P), [n_inputs, n_outputs]);
        A_imag = reshape(x_opt(3+P : 2+2*P), [n_inputs, n_outputs]);
        A_id(:,:,m) = A_real + 1i * A_imag;    
       
        RL_real = reshape(x_opt(3+2*P : 2+3*P), [n_inputs, n_outputs]);
        RL_imag = reshape(x_opt(3+3*P : 2+4*P), [n_inputs, n_outputs]);
        RL_id(:,:,m) = RL_real + 1i * RL_imag;    
       
        RH_real = reshape(x_opt(3+4*P : 2+5*P), [n_inputs, n_outputs]);
        RH_imag = reshape(x_opt(3+5*P : 2+6*P), [n_inputs, n_outputs]);
        RH_id(:,:,m) = RH_real + 1i * RH_imag;
   
        %print obtained values
        fprintf('\n--- MODO IDENTIFICATO #%d ---\n', m);
        fprintf('Frequenza Naturale (wn): %8.2f rad/s - %.4f Hz\n', omega_id(m), omega_id(m)/2/pi);
        fprintf('Smorzamento (zeta):     %8.4f (%.4f%%)\n', zeta_id(m), zeta_id(m)*100);
    end    
end

f_id = omega_id./(2*pi);

%% Plot the local result
step=5;
for k = 1:n_inputs
    figure('Name', sprintf('Input at %.1f m', x_k(k)), 'NumberTitle', 'off');
    
    for j = 1:n_outputs

        for ii=1:length(omega_id)
            %calculate the frequency range Ws
            w_lo = band_lo * omega_guess(ii);
            w_hi = band_hi * omega_guess(ii);
            idx = (W >= w_lo) & (W <= w_hi);
            cut_W = W(idx);
            
            %costruct the FRF around the resonance
            low_freq_term = RL_id(k, j, ii)./cut_W.^2;
            high_freq_term = RH_id(k, j, ii);
            res_term = A_id(k, j, ii)./(-cut_W.^2 + 1i * 2 * zeta_id(ii)*omega_id(ii).*cut_W + omega_id(ii)^2);
            Gnum = low_freq_term + high_freq_term + res_term;

            subplot(2, n_outputs, j);

            if ii == 1
                semilogy(F, abs(squeeze(FRFs(k, j, :))), 'LineWidth', 1.2);
                grid on;
                hold on;
                xlim([0 fmax]);
                title(sprintf('Output: %.1f m', x_j(j)));
                xlabel('f [Hz]');
                ylabel('Mag [m/N]');
            end
            semilogy(cut_W(1:step:end)/(2*pi), abs(Gnum(1:step:end)), 'or');
            
            subplot(2, n_outputs, j + n_outputs);
            
            if ii==1
                plot(F, angle(squeeze(FRFs(k, j, :))), 'LineWidth', 1.2);
                grid on;
                hold on;
                title(sprintf('Output: %.1f m', x_j(j)));
                xlim([0 fmax]);
                xlabel('f [Hz]');
                ylabel('Phase [rad]');
                ylim([-5 5]);
            end 
            plot(cut_W(1:step:end)/(2*pi), angle(Gnum(1:step:end)), 'or');
        end
    end
end

%% Normalize amplitudes

A_id_norm = zeros(n_outputs, length(omega_id));
if n_inputs > 1
    for ii = 1:length(omega_id)
        A_temp = squeeze(A_id(:, :, ii)).'; 
        [U, S, V] = svd(A_temp);
        s = U(:, 1);
        s_real = abs(s) .* sign(real(s));
        s_norm = s_real / max(abs(s_real));
        A_id_norm(:, ii) = s_norm;
    end
else
    for ii = 1:length(omega_id)
        % Estraggo le ampiezze per il modo ii (vettore riga 1 x 3)
        amp_sperimentali = A_id(1, :, ii);
        
        % Trovo il sensore con il segnale massimo per fare da "ancora"
        [~, index] = max(abs(amp_sperimentali));  
        valore_sperimentale = amp_sperimentali(index);
        
        % Trovo il valore analitico nello stesso punto (uso round per sicurezza)
        valore_analitico = phi(ii, round(x_j(index)*1000));
        
        % Calcolo il fattore di scala
        fattore_scala = valore_analitico / valore_sperimentale;
        
        % Scalo i punti e li salvo nella COLONNA ii (trasponendo la riga in colonna)
        A_id_norm(:, ii) = amp_sperimentali.' * fattore_scala;
    end
end

%% Plot identified and simulated mode shapes

for i = 1:length(omega_id)
    figure('Name','Mode shape: simulated vs identified' ,'Color','w');
    plot(x, real(phi(i,:)), 'b-', 'LineWidth', 1.4); hold on; grid on;
    plot(x_j, real(A_id_norm(:,i)), 'ob','MarkerSize',8,'LineWidth',2, 'MarkerFaceColor','w');
    xlabel('Distance x from the end [m]');
    ylabel(sprintf('\\phi_%d (norm.)', i));
    title(sprintf('Mode shape %d  -  f = %.2f Hz', i, f_id(i)));
    xlim([0 l]); ylim([-1.3 1.3]);
end

%% compute MAC

for i=1:length(omega_id)
    calculateMAC(real(phi(i, x_j.*1000)), real(A_id_norm(:, i)))
end


%% MAX function

function mac_value = calculateMAC(phi_id, phi_an)
    % Assicurati che entrambi i vettori siano vettori colonna
    phi_id = phi_id(:);
    phi_an = phi_an(:);
    
    % Calcolo numeratore e denominatore (l'operatore ' fa il trasposto complesso coniugato)
    numeratore = abs(phi_id' * phi_an)^2;
    denominatore = (phi_id' * phi_id) * (phi_an' * phi_an);
    
    % Risultato
    mac_value = numeratore / denominatore;
end

%% Error function used by sqnonlin

function err = errFRF(x, Om, Gexp, Nin, Nout)
    % Global parameter
    om = x(1); % natural frequency
    z  = x(2); % damping
    
    % Total couples of inputs and outputs
    P = Nin * Nout;
    M = numel(Om); % resonance frequences
    
    % residuas extraction
    A  = reshape(x(3       : 2+P),   [Nin, Nout]) + 1i*reshape(x(3+P     : 2+2*P), [Nin, Nout]);
    RL = reshape(x(3+2*P   : 2+3*P), [Nin, Nout]) + 1i*reshape(x(3+3*P   : 2+4*P), [Nin, Nout]);
    RH = reshape(x(3+4*P   : 2+5*P), [Nin, Nout]) + 1i*reshape(x(3+5*P   : 2+6*P), [Nin, Nout]);
    
    % denominator
    den = -Om.^2 + 1i*2*z*om.*Om + om^2;
    
    err = zeros(2 * M * P, 1);
    block = 0; 
    
    for i = 1:Nin
        for j = 1:Nout
            Gnum = A(i,j)./den + RL(i,j)./Om.^2 + RH(i,j);
            Gnum = Gnum(:);
            
            G_exp = squeeze(Gexp(i,j,:));
            G_exp = G_exp(:);
            
            d = G_exp - Gnum;
            
            offset = block * 2 * M;
            err(offset + (1:M))     = real(d); % real part
            err(offset + M + (1:M)) = imag(d); % imag part
            
            block = block + 1;
        end
    end
end



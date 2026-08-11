function [h_hp, f1_hp, f2_hp, HP_level] = computeHPDamping(idx, max, f, G)
    f0     = f(idx); 
    omega0 = 2*pi*f0;
    HP_level = max / sqrt(2);
    f1_hp = NaN; 
    f2_hp = NaN;
    
    % Ricerca a sinistra del picco
    for kL = idx-1:-1:1
        if G(kL) <= HP_level, f1_hp = interp1(G(kL:kL+1), f(kL:kL+1), HP_level); break; end
    end
    % Ricerca a destra del picco
    for kR = idx+1:numel(G)
        if G(kR) <= HP_level, f2_hp = interp1(G(kR-1:kR), f(kR-1:kR), HP_level); break; end
    end
 
    % Calcolo del rapporto di smorzamento (h) tramite formula quadratica delle frequenze
    if ~isnan(f1_hp) && ~isnan(f2_hp)
        h_hp = ((2*pi*f2_hp)^2 - (2*pi*f1_hp)^2)/(4*omega0^2);
        fprintf('Damping (half-power)      : h = %.5f  (%.3f %%)\n', h_hp, h_hp*100);
    else
        h_hp = NaN;
        warning('Half-power: impossibile trovare entrambi i punti a |G|max/sqrt(2)');
    end
end
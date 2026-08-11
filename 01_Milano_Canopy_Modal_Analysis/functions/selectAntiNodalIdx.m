function [antinodal_idx] = selectAntiNodalIdx(f, G, nf, percentage, N_pts)
    [~, idx_peak_global] = min(abs(f - nf));
    amp_peak = abs(G(idx_peak_global, :));     
    sorted_amps = sort(amp_peak);
    threshold   = sorted_amps(max(1, ceil(percentage * numel(amp_peak))) );
    antinodal_idx = find(amp_peak >= threshold);
    fprintf('Antinodali selezionati: %d / %d punti (soglia %.2f%%)\n', numel(antinodal_idx), N_pts, percentage);
end
function [h_pd, phase_avg, i_lo, i_hi, p] = computePDDamping(idx, f, G, nPoints)
    f0     = f(idx); 
    omega0 = 2*pi*f0;
    phase_avg = unwrap(angle(G));
    
    i_lo = max(1, idx - nPoints); 
    i_hi = min(numel(f), idx + nPoints);

    p = polyfit(2*pi*f(i_lo:i_hi), phase_avg(i_lo:i_hi), 1);
    h_pd = -1/(omega0 * p(1)); 
    fprintf('Damping (phase derivative): h = %.5f  (%.3f %%)\n', h_pd, h_pd*100);
end
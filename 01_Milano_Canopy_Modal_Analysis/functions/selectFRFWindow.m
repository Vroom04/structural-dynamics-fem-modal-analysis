function [f_win, G_win] = selectFRFWindow(nf, df, f, G)
    idx_win = f >= (nf - df) & f <= (nf + df);
    f_win   = f(idx_win);
    G_win   = G(idx_win, :); 
end
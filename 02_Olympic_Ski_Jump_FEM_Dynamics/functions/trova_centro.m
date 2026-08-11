function [C1, C2] = trova_centro(P1, P2, R)
    M = (P1 + P2) / 2;
    
    d = norm(P2 - P1);
    
    if d > 2*R
        error('Il raggio è troppo piccolo per unire i due punti!');
    end
    
    h = sqrt(R^2 - (d/2)^2);
    
    dx = P2(1) - P1(1);
    dy = P2(2) - P1(2);
    
    v_perp = [-dy, dx] / d;
    
    C1 = M + h * v_perp;
    C2 = M - h * v_perp;
end
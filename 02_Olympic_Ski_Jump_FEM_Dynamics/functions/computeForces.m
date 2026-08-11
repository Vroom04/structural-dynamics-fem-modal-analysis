function [F0] = computeForces(P, l1, l2, L, incid, n_element, ndof)
    F0 = zeros(ndof,1);
    gammaii = gamma(n_element);
    % reazioni trasversali (locale) – Hermite
    Rt1 = P*cos(gammaii) * l2^2*(3*l1+l2)/L^3;
    Rt2 = P*cos(gammaii) * l1^2*(l1+3*l2)/L^3;
    M1  = P*cos(gammaii) * l1*l2^2/L^2;
    M2  = -P*cos(gammaii) * l1^2*l2/L^2;
    % reazioni assiali (locale) – interpolazione lineare
    Ra1 = P*sin(gammaii) * l2/L;
    Ra2 = P*sin(gammaii) * l1/L;
    % trasformazione locale -> globale
    Fx1 =  Ra1*cos(gammaii) - Rt1*sin(gammaii);
    Fy1 =  Ra1*sin(gammaii) + Rt1*cos(gammaii);
    Fx2 =  Ra2*cos(gammaii) - Rt2*sin(gammaii);
    Fy2 =  Ra2*sin(gammaii) + Rt2*cos(gammaii);
    % poi assegni alle 6 posizioni incid(n_element, 1..6)
  
    if incid(n_element, 1)<= ndof
        F0(incid(n_element, 1)) = Fx1;
    end
    if incid(n_element,2)<=ndof
        F0(incid(n_element,2)) = Fy1;
    end
    if incid(n_element,3)<=ndof
        F0(incid(n_element, 3)) = M1;
    end
    
    if incid(n_element, 4)<= ndof
        F0(incid(n_element, 4)) = Fx2;
    end
    if incid(n_element,5)<=ndof
        F0(incid(n_element,5))   = Fy2;
    end
    if incid(n_element, 6)<= ndof
        F0(incid(n_element, 6)) = M2;
    end
    
end



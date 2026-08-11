function dqdt = modal_eq(t_ode, q, m, c, k, t_vec, Q_vec)
    dqdt = zeros(2,1); 
    F_t = interp1(t_vec, Q_vec, t_ode, 'linear', 0);
    dqdt(1) = q(2);                                       
    dqdt(2) = (-c*q(2)-k*q(1) + F_t) / m;               
end
function [conns] = obtainConnections(Lr, Lg, Lb, C)
    conns = {
        % TRAVI ROSSE 
        1, 2, 'line', [], Lr, 1; 
        2, 3, 'line', [], Lr, 1; 
        3, 4, 'line', [], Lr, 1; 
        4, 5, 'line', [], Lr, 1;
        6, 7, 'line', [], Lr, 1; 
        7, 8, 'line', [], Lr, 1; 
        8, 9, 'line', [], Lr, 1; 
        9, 10, 'line', [], Lr, 1;
        11, 12, 'line', [], Lr, 1; 
        12, 13, 'line', [], Lr, 1;
        13, 14, 'line', [], Lr, 1; 
        15, 16, 'line', [], Lr, 1;
        16, 17, 'line', [], Lr, 1;
        
        % TRAVI VERDI 
        %sinistra
        1, 7, 'line', [], Lg, 2; 
        2, 8, 'line', [], Lg, 2; 
        3, 9, 'line', [], Lg, 2; 
        4, 10,'line', [], Lg, 2;
        2, 7, 'line', [], Lg, 2; 
        3, 8, 'line', [], Lg, 2;  
        4, 9, 'line', [], Lg, 2;
        %destra
        12, 15, 'line', [], Lg, 2; 
        13, 16, 'line', [], Lg, 2; 
        13, 17, 'line', [], Lg, 2; 
        12, 16, 'line', [], Lg, 2;
        
        % TRAVI BLU 
        5, 10, 'line', [], Lb, 3; 
        10, 14, 'line', [], Lb, 3; 
        14, 17, 'line', [], Lb, 3;
        17, 18, 'arc', [C, 1], Lb, 3 
    };
end
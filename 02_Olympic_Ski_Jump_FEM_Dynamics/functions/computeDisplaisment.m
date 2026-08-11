function [max_d1, idx_d1, max_d2, idx_d2] = computeDisplaisment(idb, nnod, ndof, disp)
    max_d1=0;
    idx_d1=0;
    max_d2=0;
    idx_d2=0;
    for ii=1:nnod
        idx = idb(ii,2);
        if idx>ndof
            continue;
        end
        d = disp(idb(ii, 2));
        if d>max_d1
            max_d1 = d;
            idx_d1 = ii;
        end
        if d<max_d2
            max_d2 = d;
            idx_d2 = ii;
        end
    end
end
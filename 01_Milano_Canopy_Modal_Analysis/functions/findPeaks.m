function [nat_f_idx] = findPeaks(f, G, nPoints)
    nat_f_idx=[];
    for i = nPoints+1:length(f)-1-nPoints
        temp = true;
        for j=1:nPoints
            if G(i) < G(i-j) || G(i) < G(i+j)
                temp = false;
                break;
            end
        end
        if temp == true
            nat_f_idx(end+1) = i;
        end
    end
end
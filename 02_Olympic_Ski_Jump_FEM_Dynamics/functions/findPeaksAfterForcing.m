function [max, index] = findPeaksAfterForcing(t, q, totalTime)
    max = zeros(1,2);
    index = zeros(1,2);
    counter=1;
    for ii=2:length(t)-1
        if q(ii)>q(ii-1) && q(ii)>q(ii+1) && t(ii)>totalTime
            max(counter) = q(ii);
            index(counter) = ii;
            counter = counter+1;
            if counter > 2
                break; 
            end
        end
   end
end
function [beams] = sortElements(beams, coord, xy)
    for ii=1:length(beams)
        for jj=1:length(beams)-1
            if xy(beams(jj,coord)) > xy(beams(jj+1,coord))
                temp = beams(jj,:);
                beams(jj,:) = beams(jj+1,:);
                beams(jj+1,:) = temp;
            end
        end
    end
end
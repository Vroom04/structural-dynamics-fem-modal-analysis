function [MFF, MCF, MFC, MCC] = computeSubMatrices(ndof, M)
    MFF = M(1:ndof,1:ndof);
    MCF = M(ndof+1:end,1:ndof);
    MFC = M(1:ndof,ndof+1:end);
    MCC = M(ndof+1:end,ndof+1:end);
end
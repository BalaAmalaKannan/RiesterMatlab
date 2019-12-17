% --------------------------------------------------------------------
%% %% VARIED from graycomatrix - no [r,c] included in input
function oneGLRLM = computeGLRLM(si,offset,nl)
% For given direction, compute the run length matrix
switch offset
    case 1
        % 0 degree
        oneGLRLM = rle_0(si,nl);
    case 2
        % 45 degree
        seq = zigzag(si);
        oneGLRLM  = rle_45(seq,nl);
    case 3
        % 90 degree
        oneGLRLM = rle_0(si',nl);
    case 4
        % 135 degree
        seq = zigzag(fliplr(si));
        oneGLRLM = rle_45(seq,nl);
    otherwise
        error('Only 4 directions supported')
end

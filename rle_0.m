function oneglrlm = rle_0(si,NL)
% RLE   image gray level Run Length matrix for 0degree
% Assure row number is exactly the gray level

[m,n]=size(si);
oneglrlm=zeros(NL,n);
for i=1:m
    x=si(i,:);
    % run length Encode of each vector
    index = [ find(x(1:end-1) ~= x(2:end)), length(x) ];
    len = diff([ 0 index ]); % run lengths
    val = x(index);          % run values
    temp =accumarray([val;len]',1,[NL n]);% compute current numbers (or contribution) for each bin in GLRLM
    oneglrlm = temp + oneglrlm; % accumulate each contribution
end
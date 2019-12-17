function stats = grayCoprops(varargin)

allStats = {'Contrast','Correlation','Energy','Homogeneity','Variance','Entropy'};

[GLCM, requestedStats] = ParseInputs(allStats, varargin{:});

% Initialize output stats structure.
numStats = length(requestedStats);
numGLCM = size(GLCM,3);
empties = repmat({zeros(1,numGLCM)},[numStats 1]);
stats = cell2struct(empties,requestedStats,1);

for p = 1 : numGLCM
    
    if numGLCM ~= 1 %N-D indexing not allowed for sparse.
        tGLCM = normalizeGLCM(GLCM(:,:,p));
    else
        tGLCM = normalizeGLCM(GLCM);
    end
    
    % Get row and column subscripts of GLCM.  These subscripts correspond to the
    % pixel values in the GLCM.
    s = size(tGLCM);
    [c,r] = meshgrid(1:s(1),1:s(2));
    r = r(:);
    c = c(:);
    
    % Calculate fields of output stats structure.
    for k = 1:numStats
        name = requestedStats{k};
        switch name
            
            case 'Energy'
                stats.(name)(p) = calculateEnergy(tGLCM);
                
            case 'Contrast'
                stats.(name)(p) = calculateContrast(tGLCM,r,c);
                
            case 'Correlation'
                stats.(name)(p) = calculateCorrelation(tGLCM,r,c);
                
            case 'Homogeneity'
                stats.(name)(p) = calculateHomogeneity(tGLCM,r,c);
                
            case 'Variance'
                stats.(name)(p) = calculateVariance(tGLCM,r,c);
                
            case 'Entropy'
                stats.(name)(p) = calculateEntropy(tGLCM);
                
        end
    end
    
end


%-----------------------------------------------------------------------------
function glcm = normalizeGLCM(glcm)

% Normalize glcm so that sum(glcm(:)) is one.
if any(glcm(:))
    glcm = glcm ./ sum(glcm(:));
end

%-----------------------------------------------------------------------------
function C = calculateContrast(glcm,r,c)

term_1 = abs(r - c).^2;
term_2 = glcm;
contrast = term_1 .* term_2(:);

C = sum(contrast);

%-----------------------------------------------------------------------------
function Correlation = calculateCorrelation(glcm,r,c)

% means of px and py
mean_r = mean_Index(r,glcm);
mean_c = mean_Index(c,glcm);

% standard deviation of px and py
Std_r = stdIndex(r,glcm,mean_r);
Std_c = stdIndex(c,glcm,mean_c);
%--- mean and standard deviation of pixel value in the column direction, e.g.,for glcm = [0 0;1 0] mean_c is 1 and Sc is 0.

numerator_term = (r - mean_r) .* (c - mean_c) .* glcm(:);
numerator = sum(numerator_term);

Correlation = numerator / (Std_r * Std_c);

%-----------------------------------------------------------------------------
function SD = stdIndex(index,glcm,mean_index)
% Varx = Varx+(glcm_normalized(i,j)*(i-Meanx).^2);
% Vary = Vary+(glcm_normalized(i,j)*(j-Meany).^2);
standard_deviation = (index - mean_index).^2 .* glcm(:);
SD = sqrt(sum(standard_deviation)); % sqrt - as .^2 is not required

%-----------------------------------------------------------------------------
function Mean = mean_Index(index,glcm)
% Meanx = Meanx + (glcm_normalized(i,j)*i);
% Meany = Meany + (glcm_normalized(i,j)*j);
Mean = index .* glcm(:);
Mean = sum(Mean);

%-----------------------------------------------------------------------------
function E = calculateEnergy(glcm)

energy = glcm.^2;
E = sum(energy(:));

%-----------------------------------------------------------------------------
function H = calculateHomogeneity(glcm,r,c)

denominator = (1 + (r - c).^2);
homogeneity = glcm(:) ./ denominator;
H = sum(homogeneity);

%-----------------------------------------------------------------------------
function V = calculateVariance(glcm,r,c)

% means of px and py
mean_r = mean_Index(r,glcm);
mean_c = mean_Index(c,glcm);
mean_xy = (mean_r + mean_c);
% mean_xy = mean(glcm(:)); % using built-in function on whole matrix - verified with mean_r + mean_c (same result)
variance = ((r - mean_xy).^2).*(glcm(:));
V = sum(variance(:));

%-----------------------------------------------------------------------------
function Entropy = calculateEntropy(GLCM)

%--- calculate histogram counts
e = imhist(GLCM(:));
% remove zero entries in p
e(e==0) = [];
% normalize p so that sum(p) is one.
e = e ./ numel(GLCM);
entropy_1 = e(:).*log2(e(:));

Entropy = - sum(entropy_1(:));


%-------------------- To create sturcture with field name and their values -------------------------------------------
function [glcm,reqStats] = ParseInputs(allstats,varargin)

numstats = length(allstats);
narginchk(1,numstats+1);

reqStats = '';
glcm = varargin{1};

% The 'nonnan' and 'finite' attributes are not added to validateattributes because the
% 'integer' attribute takes care of these requirements.
validateattributes(glcm,{'logical','numeric'},{'real','nonnegative','integer'}, ...
    mfilename,'GLCM',1);

if ndims(glcm) > 3
    error(message('images:graycoprops:invalidSizeForGLCM'))
end

% Cast GLCM to double to avoid truncation by data type. Note that GLCM is not an
% image.
if ~isa(glcm,'double')
    glcm = double(glcm);
end

list = varargin(2:end);

if isempty(list)
    % GRAYCOPROPS(GLCM) or GRAYCOPROPS(GLCM,PROPERTIES) where PROPERTIES is empty.
    reqStats = allstats;
else
    if iscell(list{1}) || numel(list) == 1
        % GRAYCOPROPS(GLCM,{...})
        list = list{1};
    end
    
    if ischar(list)
        %GRAYCOPROPS(GLCM,SPACE-SEPARATED STRING)
        C = textscan(list, '%s');
        list = C{1};
    end
    
    anyprop = allstats;
    anyprop{end+1} = 'all';
    
    for k = 1 : length(list)
        match = validatestring(list{k}, anyprop, mfilename, 'PROPERTIES', k+1);
        if strcmp(match,'all')
            reqStats = allstats;
            break;
        end
        reqStats{k} = match;
    end
    
end

% Make sure that reqStats are in alphabetical order.
reqStats = sort(reqStats);

if isempty(reqStats)
    error(message('images:graycoprops:internalError'))
end

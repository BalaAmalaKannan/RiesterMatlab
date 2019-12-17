function stats = grayCopropsSumDiff(varargin)

allStats = {'SumAverage','SumEntropy','SumVariance','DifferenceVariance','DifferenceEntropy','IMC1','IMC2'};

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
            
            case 'SumAverage'
                stats.(name)(p) = calculateSumAverage(tGLCM,GLCM);
                
            case 'SumEntropy'
                stats.(name)(p) = calculateSumEntropy(tGLCM,GLCM);
                
            case 'SumVariance'
                stats.(name)(p) = calculateSumVariance(tGLCM,GLCM);

            case 'DifferenceVariance'
                stats.(name)(p) = calculateDifferenceVariance(tGLCM,GLCM);
                
            case 'DifferenceEntropy'
                stats.(name)(p) = calculateDifferenceEntropy(tGLCM,GLCM);
                
            case 'IMC1'
                stats.(name)(p) = calculateInformationMeasureOfCorrelation1(tGLCM,GLCM);
                
            case 'IMC2'
                stats.(name)(p) = calculateInformationMeasureOfCorrelation2(tGLCM,GLCM);
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
function SA = calculateSumAverage(glcm,GLCM)

S = size(GLCM,1);
SumAverage = zeros(1,2*S);
pxy = Pxplusy(glcm,GLCM);

for i = 2:2*S
    SumAverage(i) = i.*pxy(i);
end

% Sum Average
%SumAverage(1) = [];       % not necessary f_6(1) is zero anyway
SA = sum(SumAverage);

%-----------------------------------------------------------------------------
function SE = calculateSumEntropy(glcm,GLCM)

S = size(GLCM,1);
SumEntropy = zeros(1,2*S);
pxy = Pxplusy(glcm,GLCM);

for i = 2:2*S
    SumEntropy(i) = pxy(i).*log(pxy(i)+eps);
end

% Sum Average
SE = -sum(SumEntropy);

%-----------------------------------------------------------------------------
function SV = calculateSumVariance(glcm,GLCM)

SA = calculateSumAverage(glcm,GLCM);
S = size(GLCM,1);
SumVariance = zeros(1,2*S);
pxy = Pxplusy(glcm,GLCM);

for i=2:2*S
    SumVariance(i) = (i - SA).^2.*pxy(i);
end

% Sum Average
SV = sum(SumVariance);

%-----------------------------------------------------------------------------
function DV = calculateDifferenceVariance(glcm,GLCM)

S = size(GLCM,1);
DifferenceVariance = zeros(1,S);
px_y = Pxminusy(glcm,GLCM);

for k = 1:S
    DifferenceVariance(k) = (k-1).^2*px_y(k);
end

DV = sum(DifferenceVariance);

%-----------------------------------------------------------------------------
function DE = calculateDifferenceEntropy(glcm,GLCM)

S = size(GLCM,1);
DifferenceEntropy = zeros(1,S);
px_y = Pxminusy(glcm,GLCM);

for k = 1:S
    DifferenceEntropy(k) = px_y(k)*log(px_y(k)+eps);
end

DE = -sum(DifferenceEntropy);


%-----------------------------------------------------------------------------
function pxy = Pxplusy(glcm,GLCM)

S=size(GLCM,1);
pxy=zeros(1,2*S);

for i=1:S
    for j=1:S
        pxy(i+j) = pxy(i+j)+glcm(i,j);
    end
end

%-----------------------------------------------------------------------------
function px_y = Pxminusy(glcm,GLCM)

S=size(GLCM,1);
px_y=zeros(1,S);

for i=1:S
    for j=1:S
        px_y(abs(i-j)+1) = px_y(abs(i-j)+1)+glcm(i,j);
    end
end

%-----------------------------------------------------------------------------
function IMC1 = calculateInformationMeasureOfCorrelation1(glcm,GLCM)

% Information Measures of Correlation 1
[Hx,Hy] = HxHy(glcm,GLCM);
[HXY1,HXY2] = HXY12(glcm,GLCM);
% entropy of p(i,j)
HXY = calculateEntropy(GLCM);

numerator = HXY - HXY1;
denominator = max(Hx,Hy);

IMC1 = numerator ./ denominator; 

%-----------------------------------------------------------------------------
function IMC2 = calculateInformationMeasureOfCorrelation2(glcm,GLCM)

% Information Measures of Correlation 2
[HXY1,HXY2] = HXY12(glcm,GLCM);
% entropy of p(i,j)
HXY = calculateEntropy(GLCM);

term = (1 - exp((-2).*(HXY2 - HXY)));

IMC2 = term.^0.5; 

%-----------------------------------------------------------------------------
function [Hx,Hy] = HxHy(glcm,GLCM)
S = size(GLCM,1);
% prealloting
HX_ = zeros(1,S);
HY_ = zeros(1,S);

py = sum(glcm,1);
px = sum(glcm,2);

for i = 1:S
    for j = 1:S
        HX_(i)= px(i)*log(px(i)+eps);
        HY_(j)= py(j)*log(py(j)+eps);
    end
end
% entorpies of px and py
Hx = -sum(HX_);
Hy = -sum(HY_);

%-----------------------------------------------------------------------------
function [HXY1,HXY2] = HXY12(glcm,GLCM)
S = size(GLCM,1);

% prealloting
HXY_1 = zeros(S);
HXY_2 = zeros(S);

py = sum(glcm,1);
px = sum(glcm,2);

for i = 1:S
    for j = 1:S       
        HXY_1(i,j) = glcm(i,j)*log(px(i)*py(j)+eps);
        HXY_2(i,j) = px(i)*py(j)*log(px(i)*py(j)+eps);
    end
end
HXY1 = -sum(HXY_1(:));
HXY2 = -sum(HXY_2(:));

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

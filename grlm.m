clc;clear all;close all;

%% Read input image
dirpath = 'D:\Input Images\Bad_images\';
% dirpath = 'C:\Users\bala-amala.kannan\Dropbox\matlab\Good_Images\';
% dirpath = 'C:\Users\bala-amala.kannan\Dropbox\matlab\RiesterEpipoleEyeTest\crop\';
imagesfiles = dir(strcat (dirpath,'*.JPG')); % Creation of directory for continuous feeding of image to the system
total_images = length(imagesfiles);
tic
for i = 1:total_images
    currentfilename = imagesfiles(i).name; %read the name of the first file
        currentfilename = strcat(num2str(i),'_bad.jpg');
    currentimage = imread(strcat(dirpath,currentfilename)); %read the image
    image{i} = currentimage;
    obj=image{i};
    cropped_img = crop_img(obj);
    gray = rgb2gray(cropped_img);
    
%  Description 
%  -------------------------------------------
%%   Computes the graylevel run length (GLRL) matrix used for textural
%   basic steps:
%       Step 1 determine direction
%       Step 2 zigzag scan
%       Step 3 obtain new sequences
%       Step 4 calculate run-length matrix
%   -----------------------------------------
%     [GLRLMS,SI]= grayrlmatrix(gray);
    gray_limits = [min(gray(:)) max(gray(:))];
    % SI - scaled image. The values in SI are between 1 and NumLevels.
    [GLRLMS,SI] = grayrlmatrix(gray,'NumLevels',8,'GrayLimits',gray_limits);
% [GLRLMS1,SI1] = grayrlmatrix(gray,'NumLevels',182,'GrayLimits',gray_limits);


%%  Current supported statistics include:
%  -------------------------------------------
%   Short Run Emphasis (SRE)
%   Long Run Emphasis (LRE)
%   Gray-Level Nonuniformity (GLN)
%   Run Length Nonuniformity (RLN)
%   Run Percentage (RP)
%   Low Gray-Level Run Emphasis (LGRE)
%   High Gray-Level Run Emphasis (HGRE)
%   Short Run Low Gray-Level Emphasis (SRLGE)
%   Short Run High Gray-Level Emphasis (SRHGE)
%   Long Run Low Gray-Level Emphasis (LRLGE)
%   Long Run High Gray-Level Emphasis (LRHGE)
%  --------------------------------------------
    stats = grayrlprops(GLRLMS,SI);
%     stats1 = grayrlprops(GLRLMS1,SI1);


%% To get the final feature, take average in all directions
feature_name = char('SRE', 'LRE', 'GLN', 'RLN', 'RP', 'LGRE', 'HGRE', 'SRLGE', 'SRHGE', 'LRLGE', 'LRHGE');
final_features{i} = mean(stats);

end

%% selected feature in the base paper - dominant feature in each iteration%%
%------- SRE, LGRE, HGRE, LRHGE, SRLGE -------%

%%--- Decide 'NumLevels' value in the matrix - for 8 bit GL imgs, use 16 graylevels 
%     else max intensity ((8 for numeric, 2 for binary (default) | integer)) ---%%

%%--- verify the feature formula again -------%%
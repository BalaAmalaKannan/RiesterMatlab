clc;clear all;close all;

%% Read input image
% dirpath = 'D:\Input Images\Bad_images\';
dirpath = 'D:\Input Images\Good_images\';
% dirpath = 'C:\Users\bala-amala.kannan\Dropbox\matlab\RiesterEpipoleEyeTest\crop\';
imagesfiles = dir(strcat (dirpath,'*.JPG')); % Creation of directory for continuous feeding of image to the system
total_images = length(imagesfiles);
tic

for i = 1:total_images
    currentfilename = imagesfiles(i).name; %read the name of the first file
        currentfilename = strcat(num2str(i),'_good.jpg');
    currentimage = imread(strcat(dirpath,currentfilename)); %read the image
    image{i} = currentimage;
    obj=image{i};
    cropped_img = crop_img(obj);
    gray_img = rgb2gray(cropped_img);

    % To Create gray-level co-occurrence matrix from image
glcms = graycomatrix(gray_img,'NumLevels',max(gray_img(:)),'GrayLimits',[],'Offset',[0 1; -1 1; -1 0; -1 -1]);

% Properties of gray-level co-occurrence matrix - features extracted for texture analysis
statistical_metric(i) = grayCoprops(glcms);
other_metric(i) = grayCopropsSumDiff(glcms);
end
toc

%% verify IMC1 and IMC2(observed result is 1 always) result
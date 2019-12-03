clc;clear all;close all;

%% Read input image
dirpath = 'C:\Users\bala-amala.kannan\Dropbox\matlab\Good_Images\';
imagesfiles = dir(strcat (dirpath,'*.JPG')); % Creation of directory for continuous feeding of image to the system
n = length(imagesfiles);

for i = 1:n
    currentfilename = imagesfiles(i).name; %read the name of the first file
    currentfilename = strcat(num2str(i),'_good.JPG');
    currentimage = imread(strcat(dirpath,currentfilename)); %read the image
    image{i} = currentimage;
    obj=image{i};
    cropped_img = crop_img(obj);
    
    rgb_img = cropped_img;
a = im2double(rgb_img);
[m,n,k] = size(a);
ar = a(:,:,1);
ag = a(:,:,2);
ab = a(:,:,3);
%--------- when color is adjusted by adding scalar coefficient, the Value index is changed(increased) ----%
br = 1.4*ar;
bg = 1.4*ag;
bb = 1.4*ab;

c = zeros(m,n,k);
c(:,:,1) = br;
c(:,:,2) = bg;
c(:,:,3) = bb;

rgb_img_adj = c; 
    %--------- convert rgb cropped image to HSV color space to extract color information -----------
    % Intensities are in the range [0,1]
    
    hsv_img = (rgb2hsv(rgb_img_adj));  
    
    %% Extracting information from each layer in hsv color space (hsv - represents human perceive color)
    
    %-------- Hue - pure spectrum color. Values btw 0 to 360 ---------------
    Hue = hsv_img(:,:,1);
    [Hue_count,Hue_binLoc] = imhist(Hue);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %           Hue index is not chosen because people with different ethnic
    %           origins have different pigmentation on the retina.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %-------- Saturation - purity of the hue. pure color has saturation of 1 and tints of color have saturation less than 1 ----- ---------------
    Saturation = hsv_img(:,:,2);
    [Saturation_count,S_binLoc] = imhist(Saturation);
    %     figure
    %     imhist(Saturation);
    %     title('Saturation histogram');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %** since saturation data type is double, the bin locations are in the range of [0,1].
    % So there is no need normalization here.**
    % The index of the maximum frequency is estimated here as a feature for saturation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [max_freq_Saturation,Saturation_index] = max(Saturation_count);
    if Saturation_index == 1
        [second_max_freq_Saturation,second_Saturation_index] = max(Saturation_count(Saturation_count < max(Saturation_count)));
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The index (between [0,1]) of the maximum frequency is estimated here as a feature for saturation
    % The index values in Saturation_index([1,256]) is replaced with the values in S_binLoc([0,1])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    norm_Saturation_index(i) = S_binLoc(Saturation_index);
    norm_second_Saturation_index(i) = S_binLoc(second_Saturation_index);
    
    
    
    %-------- Value - lightness/darkness of a color. dark - closer to 0,bright - closer to 1
    Value = hsv_img(:,:,3);
    [Value_count,V_binLoc] = imhist(Value);
    %     figure
    %     imhist(Value);
    %     title('Value histogram')
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %** Explanation for normalization is same as discussed in Saturation **%
    % The index (between [0,1]) of the maximum frequency is estimated here as a feature for Value
    % The index values in Value_index([1,256]) is replaced with the values in V_binLoc([0,1]) for normalization
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [max_freq_Value,Value_index] = max(Value_count);
    if Value_index == 1
        [second_max_freq_Value,second_Value_index] = max(Value_count(Value_count < max(Value_count)));
    end
    
    norm_Value_index(i) = V_binLoc(Value_index);
    norm_second_Value_index(i) = V_binLoc(second_Value_index);
end

%% illustrate Saturation and Value indexes of different images

% plot 1 - norm_Saturation_index vs norm_Value_index
x_axis_pl1 = norm_Saturation_index;
y_axis_pl1 = norm_Value_index;
figure;
plot(x_axis_pl1,y_axis_pl1);
xlabel('saturation')
ylabel('Value')

% plot 2 - norm_second_Saturation_index vs norm_second_Value_index
x_axis_pl2 = norm_second_Saturation_index;
y_axis_pl2 = norm_second_Value_index;
figure;
labels = cellstr(num2str([1:18]'));
plot(x_axis_pl2(:),y_axis_pl2(:),'*')
text(x_axis_pl2(:),y_axis_pl2(:),labels,'VerticalAlignment','bottom','HorizontalAlignment','right')
xlabel('saturation')
ylabel('Value')

%--- Overexposed >> Saturation index - closer to 0;
%--- Underexposed >> Value index - closer to 0;
%--- normal >> Saturation and Value indexes - range [0.5 - 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------- Images here are directly converted to hsv model--------------%
%----------- Images are converted to hsv model after preprocessing --------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% verify once if index require normalization since here data type double is used    

function outputImage1 = sobel_edge_detection(cropped_img)  
    % % % ----------------------------------------------------- Green channel
    G = cropped_img(:,:,2);
    gray_d = G;
    % Apply median filter 32x32 - to reduce noise in certain images
    gray_d = double(medfilt2(gray_d,[16 16]));
    %gaussian filter with sigma 3
    gray_d = imgaussfilt(gray_d,3);
    %% Perform sharpening the image
    % gray_d = imsharpen(gray_d);   
    
    %get image high/length and width  
    [imglength_r, imgwidth_r] = size(gray_d);
    
    %zero padding in the input image
    ImageToAlgo = zeros (imglength_r + 2 , imgwidth_r + 2);
    ImageToAlgo(2:imglength_r + 1,2:imgwidth_r + 1) = gray_d;
    [imglength, imgwidth] = size(ImageToAlgo);
    
    %sobel edge matrix
    SobelFilterMatrixGx = [
        -1 0 1;
        -2 0 2;
        -1 0 1
        ];
    
    SobelFilterMatrixGy = [
        1 2 1;
        0 0 0;
        -1 -2 -1
        ];
    
    
    GxK = zeros (imglength-2,imgwidth-2);
    GyK = zeros (imglength-2,imgwidth-2);
    
    for k = 2:imglength-1
        for j=2:imgwidth-1
            
            %get image
            MatrixImage = ImageToAlgo(k-1:k+1,j-1:j+1);
            
            %Sobel mask for x-direction:
            Gx = (sum(sum(MatrixImage .* SobelFilterMatrixGx)));
            
            %Sobel mask for y-direction:
            Gy = (sum(sum(MatrixImage .* SobelFilterMatrixGy)));
            
            %The gradient of the image
            GxK(k-1,j-1) = Gx;
            GyK(k-1,j-1) = Gy;
        end
    end
    
    maggray = GxK.*GxK + GyK.*GyK;
    Thresh = sqrt(mean(mean(maggray)));
    
    maggraytemp = sqrt(maggray);
    outputImage = maggraytemp >  Thresh;
    
    se = strel('line',4,40);
    outputImage1 =imerode(outputImage,se);

end




%% ---------------- sobel_save ----------------
% % % figure(i) - medfilt2(I,[3 3]); imgaussfilt(I,2);
% % % figure1_(i) - medfilt2(I,[16 16]); imgaussfilt(I,3); better results from obseved figures
% % % figure2_(i) - medfilt2(I,[16 16]); imgaussfilt(I,2);
%%--------------------------------------------
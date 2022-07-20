clear all
close all
clc

%% 1

%% i)
I = imread('bad.bmp');
figure () ; 
imshow(I)
title ('Original Picture')
 
figure();
imshow(I); 
title ('Picture after ImRead')

figure();
imhist(I); 
title ('Original Picture Histogram')

%% ii) 

figure()
J = histeq(I);
imhist(J)
title ('Equalized Histogram')

figure();
imshow(J)
title ('Equalized Image')

%% iii) 

% Create Black and White Picture.
[counts,bins] = imhist(I);

i = 1;
tmp_sum = 0;
total_pixels = sum(counts);
half_sum = total_pixels/2;

while(1)
    tmp_sum = sum(counts(1:i));
    if(tmp_sum >= half_sum)
        break;
    end
    i = i + 1;
end

Threshold = i;
BW = I ; 

BW(find(I < Threshold)) = 0 ; 
BW(find(I >= Threshold)) = 255 ; 

NumOfWhites = length(find(BW == 0)) ;
NumOfBlacks  = length(find(BW == 255)) ;

figure() ; 
imshow(BW)
title ('Black and White Image')

figure()
imhist(BW)
title ('Black and White Histogram')


%% 2

k1 = 3; % k1 = 1 so we make it 3
k2 = 5;

I = imread('brain.gif');
figure();
imshow(I); % Print the original picture
title ('Original Picture')

% Implement median filters
M1 = medfilt2(I,[k1 k1]);
M2 = medfilt2(I,[k2 k2]); 

% Implement Gaussian filters
G1 = fspecial('Gaussian',[k1 k1],1);
G2 = fspecial('Gaussian',[k2 k2],1);
Res_G1 = imfilter(I,G1);
Res_G2 = imfilter(I,G2);

% Print the four versions of the image
figure();

subplot(2,2,1);
imshow(M1);
title ('Median Filter (k = 3)')
subplot(2,2,2);
imshow(M2);
title ('Median Filter (k = 5)')
subplot(2,2,3);
imshow(Res_G1);
title ('Gaussian Filter (k = 3)')
subplot(2,2,4);
imshow(Res_G2);
title ('Gaussian Filter (k = 5)')

%% 3
clearvars M1 M2

% Add Salt & Pepper noise
S_P_I = imnoise(I,'salt & pepper');
figure();
imshow(S_P_I);
title ('Image with Salt & Pepper Noise')

% Implement median filters
M1 = medfilt2(S_P_I,[k1 k1]);
M2 = medfilt2(S_P_I,[k2 k2]);

% Implement Average filters
A1 = fspecial('average',[k1 k1]);
A2 = fspecial('average',[k2 k2]);
Res_A1 = imfilter(S_P_I,A1);
Res_A2 = imfilter(S_P_I,A2);

% Print the four versions of the image
figure();
subplot(2,2,1);
imshow(M1);
title ('Median Filter (k = 3)')
subplot(2,2,2);
imshow(M2);
title ('Median Filter (k = 5)')
subplot(2,2,3);
imshow(Res_A1);
title ('Average Filter (k = 3)')
subplot(2,2,4);
imshow(Res_A2);
title ('Average Filter (k = 5)')

% Calculate MSE
n = size(I,1);
m = size(I,2);

mse_M1 = sum(sum(((S_P_I-M1).^2),1),2)/(m*n);
mse_M2 = sum(sum(((S_P_I-M2).^2),1),2)/(m*n);

mse_A1 = sum(sum(((S_P_I-Res_A1).^2),1),2)/(m*n);
mse_A2 = sum(sum(((S_P_I-Res_A2).^2),1),2)/(m*n);

%% 4 
clearvars M1 M2

h = [1 1 1; 1 -8 1 ; 1 1 1] ;

Frame1 = imread('Frame1.png');
Frame2 = imread('Frame2.png');

M1 = medfilt2(rgb2gray(Frame1),[k1 k1]);
M2 = medfilt2(rgb2gray(Frame2),[k1 k1]);

Res_H1 = imfilter(rgb2gray(Frame1),h);
Res_H2 = imfilter(rgb2gray(Frame2),h);

% Print the four versions of the image
figure();
subplot(2,2,1);
imshow(M2 - M1);
title ('Frame Difference with Median Filter (k = 3)')
subplot(2,2,2);
imshow(Res_H2 - Res_H1);
title ('Frame Difference with Laplacian Filter')
clear all;
close all ;
clc ;

%% 1

[I] = imread('tools.bmp');
[X,map] = imread('tools.bmp');
Fourier_A = abs(fftshift(fft2(I)));
figure()
imshow(Fourier_A);
title('FT before normalization')

tmp = log(1+ Fourier_A) ;
max_value = max(tmp(:));
c = 255/max_value ;
DFourier_A = c*log(1+ Fourier_A);

figure()
imshow((DFourier_A)./255);
title('FT after normalization (gray-scale)')

figure();
colormap(map);
colormap(jet());
image(DFourier_A);
title('FT after normalization (colored)')

%% 2

[I] = imread('horizontal.jpg');

Fourier_A = abs(fftshift(fft2(I)));
figure();
subplot(1,2,1)
imshow(I);
title('Horizontal Image FT before normalization')

tmp = log(1+ Fourier_A) ;
max_value = max(tmp(:));
c = 255/max_value ;
DFourier_A = c*log(1+ Fourier_A);

subplot(1,2,2)
imshow((DFourier_A)./255);
title('Horizontal Image FT after normalization')

[I] = imread('vertical.jpg');

Fourier_A = abs(fftshift(fft2(I)));
figure();
subplot(1,2,1)
imshow(I);
title('Vertical Image FT before normalization')

tmp = log(1+ Fourier_A) ;
max_value = max(tmp(:));
c = 255/max_value ;
DFourier_A = c*log(1+ Fourier_A);

subplot(1,2,2)
imshow((DFourier_A)./255);
title('Vertical Image FT after normalization')

%% 3

black_I = zeros(256);
middle = 256/2;
step_1 = 10;
step_2 = 30;
axes_1 = [(middle - step_1/2):((middle + step_1/2)-1)];
axes_2 = [(middle - step_2/2):((middle + step_2/2)-1)];
black_1 = black_I;
black_2 = black_I;
black_1(axes_1,axes_1) = 1;
black_2(axes_2,axes_2) = 1;
figure()
imshow(black_1);
title('Image (size = 10)')
figure()
imshow(black_2);
title('Image (size = 30)')

Fourier_A = abs(fftshift(fft2(black_1)));
figure()
subplot(1,2,1)
imshow(Fourier_A);
title('Size 10 Image FT before normalization')

tmp = log(1+ Fourier_A) ;
max_value = max(tmp(:));
c = 255/max_value ;
DFourier_A = c*log(1+ Fourier_A);

subplot(1,2,2)
imshow((DFourier_A)./255);
title('Size 10 Image FT after normalization')

Fourier_A = abs(fftshift(fft2(black_2)));
figure()
subplot(1,2,1)
imshow(Fourier_A);
title('Size 30 Image FT before normalization')

tmp = log(1+ Fourier_A) ;
max_value = max(tmp(:));
c = 255/max_value ;
DFourier_A = c*log(1+ Fourier_A);

subplot(1,2,2)
imshow((DFourier_A)./255);
title('Size 30 Image FT after normalization')

%% 4

black_1_rot = imrotate(black_1,45);
black_2_rot = imrotate(black_2,45);

figure()
imshow(black_1_rot);
title('Image (size = 10) after rotation')

figure()
imshow(black_2_rot);
title('Image (size = 30) after rotation')

Fourier_A = abs(fftshift(fft2(black_1_rot)));
figure()
subplot(1,2,1)
imshow(Fourier_A);
title('Image (size = 10) after rotation FT before normalization')

tmp = log(1+ Fourier_A) ;
max_value = max(tmp(:));
c = 255/max_value ;
DFourier_A = c*log(1+ Fourier_A);

subplot(1,2,2)
imshow((DFourier_A)./255);
title('Image (size = 10) after rotation FT after normalization')

Fourier_A = abs(fftshift(fft2(black_2_rot)));
figure()
subplot(1,2,1)
imshow(Fourier_A);
title('Image (size = 30) before rotation FT before normalization')

tmp = log(1+ Fourier_A) ;
max_value = max(tmp(:));
c = 255/max_value ;
DFourier_A = c*log(1+ Fourier_A);

subplot(1,2,2)
imshow((DFourier_A)./255);
title('Image (size = 30) after rotation FT before normalization')

%% 5

black_I = zeros(256);
middle = 256/2;
movement = 50 ;
step_1 = 10;
step_2 = 30;
axes_1 = [(middle - step_1/2) - movement:((middle + step_1/2)-1 - movement)];
axes_2 = [(middle - step_2/2) - movement:((middle + step_2/2)-1 - movement)];
black_1 = black_I;
black_2 = black_I;
black_1(axes_1,axes_1) = 1;
black_2(axes_2,axes_2) = 1;
figure()
imshow(black_1);
title('Image (size = 10) after movement')
figure()
imshow(black_2);
title('Image (size = 30) after movement')

Fourier_A = abs(fftshift(fft2(black_1)));
figure()
subplot(1,2,1)
imshow(Fourier_A);
title('Image (size = 10) after movement FT before normalization')

tmp = log(1+ Fourier_A) ;
max_value = max(tmp(:));
c = 255/max_value ;
DFourier_A = c*log(1+ Fourier_A);

subplot(1,2,2)
imshow((DFourier_A)./255);
title('Image (size = 10) after movement FT after normalization')

Fourier_A = abs(fftshift(fft2(black_2)));
figure()
subplot(1,2,1)
imshow(Fourier_A);
title('Image (size = 30) after movement FT before normalization')

tmp = log(1+ Fourier_A) ;
max_value = max(tmp(:));
c = 255/max_value ;
DFourier_A = c*log(1+ Fourier_A);

subplot(1,2,2)
imshow((DFourier_A)./255);
title('Image (size = 30) after movement FT after normalization')

%% 6

[I] = imread('tools.bmp');

figure();
imhist(I);
title('Image Histrogram before Equalization')

I_eq = histeq(I);
% black_2_eq = histeq(black_2);
figure();
imhist(I_eq);
title('Image Histrogram after Equalization')

Fourier_A = abs(fftshift(fft2(I_eq)));
tmp = log(1+ Fourier_A) ;
max_value = max(tmp(:));
c = 255/max_value ;
DFourier_A = c*log(1+ Fourier_A);

figure();
colormap(map);
colormap(jet());
image(DFourier_A);
title('Image FT after Equalization with normalization (colored)')

[X,map] = imread('tools.bmp');
Fourier_A = abs(fftshift(fft2(I)));
tmp = log(1+ Fourier_A) ;
max_value = max(tmp(:));
c = 255/max_value ;
DFourier_A = c*log(1+ Fourier_A);
figure();
colormap(map);
colormap(jet());
image(DFourier_A);
title('Image FT before Equalization with normalization (colored)')
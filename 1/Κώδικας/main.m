%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%    LAB1:MONTELA XRWMATWN                                    %
%                                                              %
%    OMADA ERGASIAS:LAB31239685                                %     
%                                                              %
%    ZAXARIOUDAKIS XRHSTOS AM:2014030056                       %
%    KOLOMVAKIS XRHSTOS    AM:2013030103                       %
%    ROPI MIKEL            AM:2011030083                       %
%                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

%% 1

[A,map] = imread('image.bmp');
factors = [0.222 0.707 0.071]; % Y = (222*R+707*G+71*B)/1000
n_levels = max(size(map)); % 256

new_map = zeros(n_levels,3);
mask = zeros(n_levels,3);

for i=1:3
    mask(:,i) = factors(i);
end

gray_column = sum(mask .* map,2);  % Calculate Y  (,2 : sum all value in same line)

for i=1:3
    new_map(:,i) = gray_column; % Gray is the same in every column of the pallete 
end

imwrite(A,new_map,'gray_image.bmp');

% Show the gray scale picture
figure(1);
imshow(A,new_map)
title('Gray Scaled Image');

%% 2

H = [0:(360/255):360]' ; 
S = (1-(1/10)).*ones(n_levels,1); % LAB31239685 => 6 + 8 + 5 = 19 => k = 10 = > k = 1 + 0 = 1
I = 0.5.*ones(n_levels,1);
HSI_map = [H S I];

RGB_map = HSI_to_RGB(H,S(1),I(1),n_levels); % HSI pallete to RGB

% Create the pointer array so that we can see the color pallete in a picture
new_image = zeros(n_levels,n_levels); % 256x256
for i=1:n_levels
    new_image(i,:) = i; 
end

imwrite(new_image,RGB_map,'rgb_colomap.bmp');

figure(2);
imshow(new_image,RGB_map)
title('Color Pallete');

%% 3

[B,map_B]=imread('image.bmp');
HSI = rgb2hsv(map_B);
Sfact = 0:0.2:1; % Fade_Factors

S_vect= HSI(:,2,:); % Saturation Vector
S_vect_new = zeros(length(S_vect));% Empty Table for new Saturation

for i=1:length(Sfact)
 
   S_vect_new=S_vect-Sfact(i); % We decrease the saturation by fade factor
   
   for j=1:length(S_vect_new)% Check if we have negative values (are out of bounds)
        if(S_vect_new(j)<0)          
         S_vect_new(j)=0;
        end
     end
                                    
   HSI(:,2,:) = S_vect_new; % We set the new Saturation to the HSI
                            
   mapB_new = hsv2rgb(HSI); % Change HSI to RGB
   
   name = sprintf('imag_Sfact%.1f.bmp',Sfact(i)); 
   imwrite(B,mapB_new,name); % We write the faded images
   
   str=sprintf('New Image When S Reduced By a Factor of  %.1f',Sfact(i));
   figure(2 + i);            % Show the image for each fade factor
   imshow(B,mapB_new)
   title(str);
   
end

%% 4

%[B,map2] = imread('land.bmp');
[B,map2] = imread('Pic_not_wb.bmp');

figure(9);
imshow(B,map2);
title('Picture Before White Balance');

r_avg = sum(map2(:,1))/(n_levels);
g_avg = sum(map2(:,2))/(n_levels);
b_avg = sum(map2(:,3))/(n_levels);

if ((r_avg~=g_avg) | (r_avg~=b_avg)) 
   map2(:,1) = (g_avg/r_avg)*map2(:,1);
   %map2(:,2) = map2(:,2);
   map2(:,3) = (g_avg/b_avg)*map2(:,3);
end

for i=1:256 %% Correct out of bounds values   
 for j=1:3 
     if(map2(i,j)>1)
        map2(i,j)=1;
     end
 end
end 

figure(10);
imwrite(B,map2,'White_balanced_image.bmp'); % We write the white balanced image
imshow(B,map2);
title('Picture After White Balance');

% % Testing White Balance (Image does not change)
% 
% r_avg = sum(map2(:,1))/(n_levels);
% g_avg = sum(map2(:,2))/(n_levels);
% b_avg = sum(map2(:,3))/(n_levels);
% 
% if ((r_avg~=g_avg) | (r_avg~=b_avg)) 
%    map2(:,1) = (g_avg/r_avg)*map2(:,1);
%    %map2(:,2) = map2(:,2);
%    map2(:,3) = (g_avg/b_avg)*map2(:,3);
% end
% 
% figure(11);
imwrite(B,map2,'White_balanced_image.bmp'); % We write the white balanced image
% imshow(B,map2);
% title('Picture After White Balance 2');
clear all;
clc; 
close all;

%% Create the corrupted image %%

f=imread('lena.bmp');
x = im2double(f);

[A,B]=size(f);

figure()
imshow(f)
title('Original image');

awgn_var = 0.01 ;
awgn_mu = 0 ; 

ImageWithNoise = imnoise(f,'gaussian',awgn_mu,awgn_var); %% Add noise

h = fspecial('average'); % Add blur
g=imfilter(ImageWithNoise,h);

figure()
imshow(g)
title('Corrupted Image');

%% Wiener Filter %%

H_uv = fft2(h,A,B);
F_uv =  fft2(f) ;
G_uv= fft2(g) ;

H_conj =conj(H_uv);
H_abs =(abs(H_uv)).^2;

Sf = (abs(F_uv)).^2;
Sf = log(Sf) ;
Sn = awgn_var ;   
P = Sn./Sf ; 

F_hat = (H_conj./(H_abs+ P)).*G_uv;
f_hat1 = ifft2(F_hat);

figure()
imshow((uint8(f_hat1)))
title('Wiener Filter');

%% CLSR filtering %%

P=[0 -1 0; -1 4 -1; 0 -1 0];
M=A+2;
N=B+2;

Pe=zeros(M,N);
Pe(M/2-1:M/2+1,N/2-1:N/2+1) = P;

P_abs =(abs(fft2(Pe))).^2;
var_est= (M-1)*(N-1)*(awgn_var + awgn_mu*awgn_mu) ;

H_uv=fft2(h,M,N);
F_uv =  fft2(f) ;
G_uv=fft2(g,M,N);
H_conj=conj(H_uv);
H_abs=(abs(H_uv)).^2;

gamma = 0.0000001;
iterations = 1;
alpha=0.01*var_est;

while true
    
    cor_factor = (1/3)*gamma ;
    
    F_hat=(H_conj./(H_abs+(gamma.*P_abs))).*G_uv;
    f_hat2=ifft2(F_hat);
    phi_gamma=(norm(ifft(G_uv-(H_uv.*F_hat)))).^2;
    
    if (phi_gamma< var_est-alpha)
        gamma=gamma+cor_factor;
    end
        
    if (phi_gamma >= var_est + alpha) 
        gamma=gamma-cor_factor;
    end
      
    if ((phi_gamma> var_est-alpha) && (phi_gamma< var_est+alpha))
        break;
    end
    
    iterations=iterations+1;
end

figure()
imshow(uint8(f_hat2));
title('Constrained Least Squares Restoration');

MSE_Wiener = sum(sum(((double(f)-f_hat1).^2),1),2)/(A*B)
f_hat2 =f_hat2(1:A,1:B); % Decrease dimensions to AxB
MSE_CLSR = sum(sum(((double(f)-f_hat2).^2),1),2)/(A*B)
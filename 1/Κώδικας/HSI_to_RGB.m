function [RGB_map] = HSI_to_RGB(H,S,I,n_levels)

    RGB_map = zeros(n_levels,3);
    
    if(I<=0.5)
        M = I*(1+S);
    else
        M = I + S - I*S;
    end
    
    m = 2*I - M;
    
    for i=1:n_levels
        % Filling R component (RGB_map(:,1))
        if (H(i) >= 0) && (H(i) < 60) 
            RGB_map(i,1) = m + (M - m)*(H(i)/60);
        elseif (H(i) >= 60) && (H(i) < 180)
            RGB_map(i,1) = M;
        elseif (H(i) >= 180) && (H(i) < 240)
            RGB_map(i,1) = m + (M - m)*((240 - H(i))/60);
        else
            RGB_map(i,1) = m;
        end
        
        % Filling G component (RGB_map(:,2))
        if (H(i) >= 0) && (H(i) < 120)
            RGB_map(i,2) = m;
        elseif (H(i) >= 120) && (H(i) < 180)
            RGB_map(i,2) = m + (M - m)*((H(i) - 120)/60);
        elseif (H(i) >= 180) && (H(i) < 300)
            RGB_map(i,2) = M;
        else
            RGB_map(i,2) = m + (M - m)*((360 - H(i))/60);
        end
    
        % Filling B component (RGB_map(:,3))
        if (H(i) >=0) && (H(i) < 60) 
            RGB_map(i,3) = M;
        elseif (H(i) >= 60) && (H(i) < 120)
            RGB_map(i,3) = m + (M - m)*((120 - H(i))/60);
        elseif (H(i) >= 120) && (H(i) < 240)
            RGB_map(i,3) = m;
        elseif (H(i) >= 240) && (H(i) < 300)
            RGB_map(i,3) = m + (M - m)*((H(i) - 240)/60);
        else
            RGB_map(i,3) = M;
        end
end


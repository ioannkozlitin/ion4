function [xe_part,xi] = concentrations(Z,T,V,x, xe)
global PHI;

el_num = length(Z);

if(el_num == 1)
    if( xe < 0 )
        xi = NaN;
        xe_part = NaN;
    else
        K = floor(xe);
        xi = zeros(Z,1);
        if( K > 0 )
            xi(K) = K + 1 - xe;
        else
            xi(1) = 0;
        end
        xi(K+1) = xe - K;
        xe_part = xe;
    end
else
    
    xe_part = zeros(1,el_num);
    xi = zeros(max(Z),el_num);
    
    a = 2*V*(T/2/pi)^1.5;
    for j = 1:el_num
        Phi_temp = T*log(a/xe);
        
        if( Phi_temp <= PHI(1,j) )
            temp0 = (Phi_temp - PHI(1,j))/T;
            gamma = exp( temp0 );
            xe_part(j) = gamma/(gamma+1);
        end
        if( Phi_temp >= PHI(Z(j),j) )
            temp0 = (Phi_temp - PHI(Z(j),j))/T;
            gamma = exp( temp0 );
            xe_part(j) = Z(j) - 1/(gamma+1);
        end
        
        for k = 1:Z(j)-1
            check1 = PHI(k,j) <= Phi_temp;
            check2 = Phi_temp <= 0.5*PHI(k+1,j) + 0.5*PHI(k,j);
            case1 = check1 && check2;
            if( case1 )
                temp0 = (Phi_temp - PHI(k,j))/T;
                gamma = exp( temp0 );
                temp1 = ( PHI(k+1,j) - PHI(k,j) )/2/T;
                eps = 1/( exp( temp1 ) - 1 );
                xe_part(j) = k - ( 1 - eps*(gamma-1) )/( gamma + 1 );
            end
        end
        for k = 2:Z(j)
            check3 = PHI(k,j) >= Phi_temp;
            check4 = Phi_temp >= 0.5*PHI(k-1,j) + 0.5*PHI(k,j);
            case2 = check3 && check4;
            if( case2 )
                temp0 = (Phi_temp - PHI(k,j))/T;
                gamma = exp( temp0 );
                temp1 = ( PHI(k,j) - PHI(k-1,j) )/2/T;
                eps = 1/( exp( temp1 ) - 1 );
                xe_part(j) = k - ( 1 - eps*(gamma-1) )/( gamma + 1 );
            end
        end
        K = floor( xe_part(j) );
        if( K > 0)
            xi(K,j) = x(j)*( K + 1 - xe_part(j) );
        else
            xi(1,j) = 0;
        end
        xi(K+1,j) = x(j)*( xe_part(j) - K );
    end
    
end
end
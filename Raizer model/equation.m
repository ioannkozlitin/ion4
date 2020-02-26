function f = equation(Z,T,V,x, xe)
global PHI;

improved = 1;

Phi_temp = 0;

el_num = length(Z);

if(el_num == 1)
    if( improved == 0 )
        if( (0 <= xe) && (xe <= 0.5) )
            Phi_temp = 2*PHI(1);
        end
        if( (Z-0.5 <= xe) && (xe <= Z) )
            Phi_temp = ( PHI(Z) - PHI(Z-1) )*(xe - Z+0.5) + PHI(Z);
        end
        for k = 1:Z-1
            if( (k-0.5 <= xe) && (xe <= k) )
                Phi_temp = ( PHI(k+1) - PHI(k) )*(xe - k+0.5) + PHI(k);
            end
        end
    else
        if( (0 <= xe) && (xe <= 0.5) )
            Phi_temp = PHI(1) + T*log(xe/(1-xe));
        end
        if( (Z-0.5 <= xe) && (xe <= Z) )
            Phi_temp = PHI(Z) + T*log((xe+1-Z)/(Z-xe));
        end

        for k = 1:Z-1
            if( (k-0.5 <= xe) && (xe <= k) )
                temp1 = (PHI(k+1) - PHI(k))/(2*T);
                eps = 1/(exp(temp1) - 1);
                temp2 = (xe+1-k+eps)/(k-xe+eps);
                Phi_temp = PHI(k) + T*log(temp2);
            end
        end
        for k = 2:Z
            if( ( k-1 <= xe) && (xe <= k-0.5) )
                temp = (PHI(k) - PHI(k-1))/(2*T);
                eps = 1/(exp(temp) - 1);
                temp2 = (xe+1-k+eps)/(k-xe+eps);
                Phi_temp = PHI(k) + T*log(temp2);
            end
        end
    end

    a = 2*V*(T/2/pi)^1.5;
    f = Phi_temp - T*log( a/xe );
else
    xe_part = zeros(1,el_num);
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
    end
    sum = 0;
    for j = 1:el_num
       sum = sum + x(j)*xe_part(j);
    end
    f = xe - sum;    
end
end
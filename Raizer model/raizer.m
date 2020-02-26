function [xe,xi,R,ZH,ZD] = raizer(Z,x, T,rho)
global M_EL;

rho = rho/11.2;
T = T/27.26;

el_num = length(Z);
data(Z);

mass = 0;
for j = 1:el_num
    mass = mass + x(j)*M_EL(j);% + x(j)*M_EL(j);
end

V = mass/rho;
R = (3*V/4/pi)^(1/3); 

Left = 0;
Right = max(Z);
eps = 1e-6;

[answer, delta, iternum] = dichotomy(Z,T,V,x, Left,Right,eps); %#ok<ASGLU>
xe = answer;
[xe_part,xi] = concentrations(Z,T,V,x, xe); %#ok<ASGLU>

check = isnan(xi);
if( check ~= 1 )
    sum_H = 0;
    sum_D = 0;
    for j = 1:el_num
        for i = 1:Z(j)
            sum_H = sum_H + xi(i,j)*i^1.5;
            sum_D = sum_D + xi(i,j)*i^ 2 ;
        end
    end
    sum_H = sum_H + xe;
    sum_D = sum_D + xe;
    ZH = sum_H^(2/3);
    ZD = sum_D^(1/2);
else
    ZH = NaN;
    ZD = NaN;
end
end
clear all; close all; clc;

Z = [18,36];
x = [0.5,0.5]; 
%Z = 1;
%x = 1;

compare = 0;

% Rho_all = [0.1, 1];
Rho_all = 1;

for r = 1:length(Rho_all)
rho = Rho_all(r);
if( compare == 1 )
    Saha;
else
    T_range = [3 2];%:1:3;
end
N = length(T_range);
xe_ans = zeros(1,N);
ZH_ans = zeros(1,N);
ZD_ans = zeros(1,N);
xi_ans = zeros(N,sum(Z));
el_num = length(Z);

T_max_R = zeros(el_num,max(Z));
xi_max_R  = zeros(el_num,max(Z));
delta_T_R = zeros(el_num,max(Z));

for n = 1:N
    T = 10^T_range(n);
    [xe,xi, R, ZH,ZD] = raizer(Z,x, T,rho);
    check = isnan(xi);
    if( check ~= 1 )
        
        xe_ans(n) = xe;
        ZH_ans(n) = ZH;
        ZD_ans(n) = ZD;
        
        add = 0;
        for num = 1:el_num
            for j = 1:Z(num)
                pos = j + add;
                xi_ans(n, pos ) = ( xi(j,num) );
            end
            
            add = Z(num);
        end
        
    else
        
        xe_ans(n) = NaN;
        ZH_ans(n) = NaN;
        ZD_ans(n) = NaN;
        for num = 1:el_num
            for j = 1:Z(num)
                xi_ans(n, j ) = NaN;
            end
        end
        
    end
end

add = 0;
for num = 1:el_num
    for j = 1:Z(num)
        pos = j + add;
        xi_temp = xi_ans(:,pos);
        [T_max_R(num,j), xi_max_R(num,j), delta_T_R(num,j)] = analyze(xi_temp,T_range);
    end
    add = Z(num);
end
density = Rho_all(r)
Temp_16 = T_max_R(1,16)
Temp_17 = T_max_R(1,17)

if(N >= 2)
    figure; hold on;
    plot(T_range,xe_ans,'-b')
    plot(T_range,ZH_ans,'-k',T_range,ZD_ans,'-r')
    add = 0;
    for num = 1:el_num     
        figure; hold on;
        for j = 1:Z(num)
            pos = j + add;
            plot(T_range,xi_ans(:,pos))
        end
        add = Z(num);
    end
%     plot(T_range,xi_ans(:,16));
%     plot(T_range,xi_ans(:,17))
end

title = num2str(Z(2));
file = strcat(title,'.txt');
fid = fopen(file,'w');

for n = 1:N
    fprintf(fid,'%2.4f ',T_range(n));
    for num = 1:Z(2)
%         fprintf(fid,'%2.4f ',xi_ans(n,num));
        fprintf(fid,'%2.4f ',xi_ans(n,num+Z(1)));
    end
    fprintf(fid,'\r\n');
end
fclose(fid);
end

xe_ans

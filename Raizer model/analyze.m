function [T_max, xi_max, delta_T] = analyze(xi,T_range)
points_num = length(xi);

[xi_max,Num_max] = max(xi);
T_max = T_range(Num_max);

check1 = zeros(1,Num_max);
for m = 1:Num_max
    check1(m) = ( xi(m) - xi_max/2 )^2;
end
[~,Num_right] = min(check1);
T_right = T_range(Num_right);
xi_right = xi(Num_right); %#ok<NASGU>

% figure; hold on;
% plot(T_range,xi)

if( points_num > Num_max )
    check2 = zeros(1,points_num - Num_max);
    for m = 1:points_num - Num_max
        check2(m) = ( xi(Num_max + m) - xi_max/2 )^2;
    end

    [~,Num_left] = min(check2);
    T_left = T_range(Num_left+Num_max);
    xi_left = xi(Num_left+Num_max); %#ok<NASGU>

    delta_T = abs( T_right - T_left );
else
    delta_T = 0;
end
end
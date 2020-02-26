% lgTeV lgRho xe zd zh x0 x1 ... x18

file_title = 'argonSaha01.txt';
fid = fopen(file_title,'r');
S2 = fscanf(fid,'%g');
fclose(fid);

points_num = length(S2)/24;

data = zeros(points_num,24);
for m = 1:points_num
    for n = 1:24
        data(m,n) = S2( n + 24*(m-1) );
    end
end
T_range = data(:,1);
xe = data(:,3);
ZD = data(:,4);
ZH = data(:,5);

%figure; hold on;
xi = zeros(points_num,18);
T_max_S = zeros(1,18);
xi_max_S  = zeros(1,18);
delta_T_S = zeros(1,18);

for n = 1:18
    for m = 1:points_num
        xi(m,n) = data(m,6+n);
    end
    
    xi_temp = xi(:,n);
    [T_max_S(n), xi_max_S(n), delta_T_S(n)] = analyze(xi_temp,T_range);
    
    plot(T_range,xi(:,n))
end

% figure; hold on;
% plot(T_range, xe, '-b')
% plot(T_range, ZH, '-k', T_range, ZD, '-r')
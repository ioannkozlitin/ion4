clear; clc; close all;


Z = [18,36];
X = [0.5,0.5]; 

lg_T = -2.5:0.1:4.6;
lg_rho = -6:0.1:6;

xe = [];
for i = 1:length(lg_T)
  for j = 1:length(lg_rho)
      T = 10^lg_T(i);
      rho  = 10^lg_rho(j);
      xe(i, j) = raizer(Z, X, T, rho);
  end
end

save('xe_res', 'xe');
clear; clc; close all;

res_path = '../debug';

current_dir = pwd;
cd(res_path);
CuNew;
xe1 = xe_Saha;
Raizer;
xe2 = xe_Saha;
cd(current_dir);

figure(); hold on;
plot(lgT, xe1(:, 1), '-b', 'LineWidth', 2);
plot(lgT, xe2(:, 1), '-r', 'LineWidth', 2);
xlabel('lgT'); ylabel('x_e');
title(['lgV = ' num2str(lgV(1))]);
legend('Saha', 'Raizer');
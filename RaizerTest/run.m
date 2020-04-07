clear; clc; close all;
pkg load io;

program_path = '../debug/ion4';
output_file_name = 'res.xlsx';

for i = 1:100
  Z{i} = i;
  X{i} = 1;
end

Z{end + 1} = [18, 36];
X{end + 1} = [0.5, 0.5];
Z{end + 1} = [3, 1];
X{end + 1} = [0.5, 0.5];

lgT_min = -2.5;
lgT_max = 4.6;
lgT_step = 0.1;

lgV = 1:6;

printf(" lgV  |");
for i = 1:length(lgV)
  len = printf(" %.2f", lgV(i));
  printf("%s |", 32 * ones(1, 15 - len));
end
printf("\n%s\n", '-' * ones(1, 17 * length(lgV) + 7));

output{1, 1} = 'lgV';
for j = 1:length(lgV)
  output{1, j + 1} = sprintf('%.2f', lgV(j));
end

for i = 1:length(Z)
  nameStr = get_name(Z{i});
  printf(" %4s |", nameStr);
  output{i + 1, 1} = nameStr;
  for j = 1:length(lgV)
    [rms_error, max_error] = get_error(program_path, Z{i}, X{i}, lgT_min, lgT_max, lgT_step, lgV(j));
    len = printf(" %.2f (%.2f)", rms_error, max_error);
    printf("%s |", ' '* ones(1, 15 - len));
    output{i + 1, j + 1} = sprintf('%.2f (%.2f)', rms_error, max_error);
  end
  printf("\n"); fflush(stdout);
end

xlswrite(output_file_name, output);

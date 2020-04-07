clear; clc; close all;

program_path = '../debug/ion4';

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
printf("\n");

printf("%s\n", '-' * ones(1, 17 * length(lgV) + 7));

for i = 1:length(Z)
  nameStr = get_name(Z{i});
  printf(" %4s |", nameStr);
  for j = 1:length(lgV)
    [rms_error, max_error] = get_error(program_path, Z{i}, X{i}, lgT_min, lgT_max, lgT_step, lgV(j));
    len = printf(" %.2f (%.2f)", rms_error, max_error);
    printf("%s |", ' '* ones(1, 15 - len));
  end
  printf("\n"); fflush(stdout);
end
function [rms_error, max_error] = get_error(program_path, Z, X, lgT_min, lgT_max, lgT_step, lgV)
  command = [program_path];
  command = [command " " num2str(lgT_min) " " num2str(lgT_max) " " num2str(lgT_step)];
  command = [command " " num2str(lgV)];
  
  for i = 1:length(Z)
    command = [command " " num2str(Z(i)) " " num2str(X(i))];
  end

  [~, output] = system(command);
  
  ind = strfind(output, 'output = ');
  output = output(ind:end);
  error = str2num(output);
  rms_error = error(1);
  max_error = error(2);
end
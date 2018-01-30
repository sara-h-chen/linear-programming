function [x, gval] = shipping(values)
  if ~isvector(values) || length(values) ~= 3
		error('Input must be a vector of length 3');
  end
  if any(mod(values,10))
    error('Inputs must be divisible by 10');
  end
  % Change unit of measurement to 100s of kilos
  % As that is the atomic weight
  % So the values should be for 100 kilos
  values = values / 10;
  
  % Make all x integer
  intcon = 1:12;
  % Make all x >=0
  lb = zeros(12,1);
  % Set upper bound of x (required)
  ub(1:12) = Inf;
  
  % Utility
  % Coefficients are [1F, 1C, 1B, 2F, 2C ... 4C, 4B]
  resource_indices = zeros(4, 12);
  resource_indices(1, :) = [1 1 1 0 0 0 0 0 0 0 0 0];
  resource_indices(2, :) = [0 0 0 1 1 1 0 0 0 0 0 0];
  resource_indices(3, :) = [0 0 0 0 0 0 1 1 1 0 0 0];
  resource_indices(4, :) = [0 0 0 0 0 0 0 0 0 1 1 1];
  
  compartment_indices = zeros(3, 12);
  compartment_indices(1, :) = [1 0 0 1 0 0 1 0 0 1 0 0];
  compartment_indices(2, :) = [0 1 0 0 1 0 0 1 0 0 1 0];
  compartment_indices(3, :) = [0 0 1 0 0 1 0 0 1 0 0 1];
  
  % Volumes are divided by 10 because of the change of unit
  volumes = [7 10 8.5 6];
  volume_mask = repelem(volumes, 3);
  
  % (1) No more than available weight
  A = [resource_indices(1, :);
       resource_indices(2, :);
       resource_indices(3, :);
       resource_indices(4, :)];
  b = [200; 160; 250; 130];
  
  % (2) Equal weight distribution
  Aeq = [15 * compartment_indices(1, :) - 10 * compartment_indices(2, :);
       15 * compartment_indices(1, :) - 18 * compartment_indices(3, :)];
  beq = [0; 0];
  
  % (3) No more than allowed weight per compartment
  A = [A;
       compartment_indices(1, :);
       compartment_indices(2, :);
       compartment_indices(3, :)];
  b = [b;
       120;
       180;
       100];
   
  % (4) No more than allowed volume per compartment
  A = [A;
       volume_mask .* compartment_indices(1, :);
       volume_mask .* compartment_indices(2, :);
       volume_mask .* compartment_indices(3, :)];
  b = [b;
       1000;
       1300;
       700];
  
  x = zeros(3, 4, values(1));
  gval = zeros(values(1), 1);
  for i=1:values(1)
    values(4) = i;
    % 3 components per cargo - for front, center and back
    % So coefficients are [1F, 1C, 1B, 2F, 2C ... 4C, 4B]
    f_coeff = repelem(values, 3);
    % Minimizes -dot(f_coeff, x), so maximizes dot(f_coeff, x)
    f = -f_coeff;
    
    result = intlinprog(f, intcon, A, b, Aeq, beq, lb, ub);
    gval(i) = f_coeff * result;
    % Scale coefficients down to get from 100kg -> 1 ton
    x(:, :, i) = reshape(result / 10, [3 4]);
  end
  
  % Plot graph  
  plot(gval, '-o');
  grid on;
  set(gca, 'XTickLabel', 10:10:values(1)*10);
  xlabel('value(4)')
  ylabel('gval')
  title('Optimal shipping value');
  
end

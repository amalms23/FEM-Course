% Initialize global stiffness matrix and load vector
K = zeros(num_elements + 1, num_elements + 1);
F_load = zeros(num_elements + 1, 1);

% Assemble stiffness matrix and load vector
for i = 1:num_elements
    % Element stiffness matrix
    k_elem = (E * I) / dx^3 * [12, 6 * dx, -12, 6 * dx; ...
                                6 * dx, 4 * dx^2, -6 * dx, 2 * dx^2; ...
                                -12, -6 * dx, 12, -6 * dx; ...
                                6 * dx, 2 * dx^2, -6 * dx, 4 * dx^2];
    
    % Assemble into global stiffness matrix
    K(i, i) = K(i, i) + k_elem(1, 1);
    K(i, i+1) = K(i, i+1) + k_elem(1, 2);
    K(i+1, i) = K(i+1, i) + k_elem(2, 1);
    K(i+1, i+1) = K(i+1, i+1) + k_elem(2, 2);

    % Apply point load
    if i == num_elements
        F_load(i+1) = F;
    end
end

% Apply boundary conditions (fixed end)
K(1, :) = 0;
K(:, 1) = 0;
K(1, 1) = 1;

% Solve for displacements
U = K \ F_load;

% Calculate strains and stresses
strains = diff(U) / dx;
stresses = (E * I) * strains;

% Calculate the final node positions for the deformed configuration
x_deformed = x + U;

% Plot the initial and deformed configurations
figure;
subplot(2, 1, 1);
plot(x, zeros(size(x)), 'b-', 'LineWidth', 2); % Initial configuration
xlabel('Position along the beam');
ylabel('Displacement');
title('Initial and Deformed Configurations');
grid on;

subplot(2, 1, 2);
plot(x_deformed, zeros(size(x_deformed)), 'r-', 'LineWidth', 2); % Deformed configuration
xlabel('Position along the beam');
ylabel('Displacement');
legend('Initial', 'Deformed');
grid on;
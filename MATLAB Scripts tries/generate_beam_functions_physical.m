% --- GENERATOR SCRIPT (using PHYSICAL coordinates) ---
% Derives analytical expressions for f_int and K_T using the coordinate x.
% Run this script ONCE to create the 'beam_functions.mat' file.

clear; clc;

fprintf('Generating symbolic beam functions using physical coordinates...\n');

% 1. Declare symbolic variables
syms x L E A I real; % Integration variable is now 'x'
syms u1 w1 th1 u2 w2 th2 real;
d_sym = [u1; w1; th1; u2; w2; th2];

% 2. Get shape functions and derivatives symbolically in 'x'
[Nu, Nw, dNu_dx, dNw_dx, d2Nw_dx2] = shape_functions_beam_physical(x, L);

% 3. Form symbolic B-matrices
Bu_vec = [dNu_dx(1), 0, 0, dNu_dx(2), 0, 0];
Bw1_vec = [0, dNw_dx(1), dNw_dx(2), 0, dNw_dx(3), dNw_dx(4)];
Bw2_vec = [0, d2Nw_dx2(1), d2Nw_dx2(2), 0, d2Nw_dx2(3), d2Nw_dx2(4)];

% 4. Define strains, stresses, and integrands symbolically
du_dx_sym = Bu_vec * d_sym;
dw_dx_sym = Bw1_vec * d_sym;
d2w_dx2_sym = Bw2_vec * d_sym;

Nxx_sym = E*A * (du_dx_sym + 0.5 * dw_dx_sym^2);
Mxx_sym = E*I * d2w_dx2_sym;

% Integrand for the internal force vector
B_NL_u_sym = Bu_vec + dw_dx_sym * Bw1_vec;
f_int_integrand = B_NL_u_sym' * Nxx_sym + Bw2_vec' * Mxx_sym; % No Jacobian needed

% Integrands for the tangent stiffness matrix
g_sym = dw_dx_sym;
KT_material_integrand = Bu_vec'*(E*A)*Bu_vec + Bw2_vec'*(E*I)*Bw2_vec;
KT_geometric_integrand = Bw1_vec' * Nxx_sym * Bw1_vec;
KT_nonlinear_integrand = Bw1_vec'*(E*A*g_sym^2)*Bw1_vec + ...
                         Bu_vec'*(E*A*g_sym)*Bw1_vec + ...
                         Bw1_vec'*(E*A*g_sym)*Bu_vec;
KT_integrand = KT_material_integrand + KT_geometric_integrand + KT_nonlinear_integrand;

% 5. Perform symbolic integration over [0, L]
fprintf('Performing symbolic integration over [0, L]...\n');
f_int_expr = int(f_int_integrand, x, 0, L);
K_T_expr = int(KT_integrand, x, 0, L);

% 6. Convert expressions to optimized MATLAB functions
fprintf('Converting expressions to MATLAB functions...\n');
vars = {L, E, A, I, u1, w1, th1, u2, w2, th2};
f_int_handle = matlabFunction(simplify(f_int_expr), 'Vars', vars);
K_T_handle = matlabFunction(simplify(K_T_expr), 'Vars', vars);

% 7. Save the function handles
save('beam_functions.mat', 'f_int_handle', 'K_T_handle');

fprintf('Symbolic functions saved to beam_functions.mat\n');
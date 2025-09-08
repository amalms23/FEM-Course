function [f_int, K_T] = element_routine_beam_symbolic(L, E, A, I, d_elem)
% Computes f_int and K_T using SYMBOLIC integration.
% NOTE: This is computationally slow and for demonstration only.

    % 1. Declare symbolic variables
    syms xi_sym real; % Symbolic integration variable
    d_sym = sym('d', [6, 1]); % Symbolic displacement vector [d1; d2; ... d6]
    
    % 2. Get shape functions and their derivatives symbolically
    [~, ~, dNu_dxi, dNw_dxi, d2Nw_dxi2] = shape_functions_beam(xi_sym, L);
    
    J = L / 2; % Jacobian is constant for a straight beam
    dNu_dx = (1/J) * dNu_dxi;
    dNw_dx = (1/J) * dNw_dxi;
    d2Nw_dx2 = (1/J)^2 * d2Nw_dxi2;
    
    % 3. Form symbolic B-matrices
    Bu_vec = [dNu_dx(1), 0, 0, dNu_dx(2), 0, 0];
    Bw1_vec = [0, dNw_dx(1), dNw_dx(2), 0, dNw_dx(3), dNw_dx(4)];
    Bw2_vec = [0, d2Nw_dx2(1), d2Nw_dx2(2), 0, d2Nw_dx2(3), d2Nw_dx2(4)];

    % 4. Form symbolic strains and stress resultants
    du_dx_sym = Bu_vec * d_sym;
    dw_dx_sym = Bw1_vec * d_sym;
    d2w_dx2_sym = Bw2_vec * d_sym;

    Nxx_sym = E*A * (du_dx_sym + 0.5 * dw_dx_sym^2);
    Mxx_sym = E*I * d2w_dx2_sym;
    
    % 5. Define the integrands for f_int and K_T
    B_NL_u_sym = Bu_vec + dw_dx_sym * Bw1_vec;
    
    % Integrand for the internal force vector
    f_int_integrand = (B_NL_u_sym' * Nxx_sym + Bw2_vec' * Mxx_sym) * J;
    
    % Integrands for the tangent stiffness matrix
    g_sym = dw_dx_sym;
    KT_material_integrand = (Bu_vec'*(E*A)*Bu_vec + Bw2_vec'*(E*I)*Bw2_vec) * J;
    KT_geometric_integrand = (Bw1_vec' * Nxx_sym * Bw1_vec) * J;
    KT_nonlinear_integrand = (Bw1_vec'*(E*A*g_sym^2)*Bw1_vec + ...
                              Bu_vec'*(E*A*g_sym)*Bw1_vec + ...
                              Bw1_vec'*(E*A*g_sym)*Bu_vec) * J;

    KT_integrand = KT_material_integrand + KT_geometric_integrand + KT_nonlinear_integrand;
    
    % 6. Perform symbolic integration
    f_int_expression = int(f_int_integrand, xi_sym, -1, 1);
    K_T_expression = int(KT_integrand, xi_sym, -1, 1);
    
    % 7. Substitute numerical values for the current state
    f_int = double(subs(f_int_expression, d_sym, d_elem));
    K_T   = double(subs(K_T_expression, d_sym, d_elem));
    
end
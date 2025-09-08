function [f_int, K_T] = element_routine_beam(L, E, A, I, d_elem)
% Computes the internal force vector and tangent stiffness matrix
% for a 2-node nonlinear Euler-Bernoulli beam element.

    ndof = 6;
    f_int = zeros(ndof, 1);
    K_T = zeros(ndof, ndof);
    
    % 2-point Gauss Quadrature
    gauss_points = [-1/sqrt(3), 1/sqrt(3)];
    gauss_weights = [1, 1];
    
    for i = 1:length(gauss_points)
        xi = gauss_points(i);
        w_gp = gauss_weights(i);
        
        % Get shape functions and derivatives in natural coordinate system
        [Nu, Nw, dNu_dxi, dNw_dxi, d2Nw_dxi2] = shape_functions_beam(xi, L);
        
        % Jacobian of transformation: dx = J * dxi
        J = L / 2;
        
        % Derivatives with respect to physical coordinate x
        dNu_dx = (1/J) * dNu_dxi;
        dNw_dx = (1/J) * dNw_dxi;
        d2Nw_dx2 = (1/J)^2 * d2Nw_dxi2;
        
        % B-matrices (mapping nodal DOFs to strains/curvatures)
        % These map the full 6-DOF element vector to a single value at the GP
        Bu_vec = zeros(1, ndof);
        Bw1_vec = zeros(1, ndof);
        Bw2_vec = zeros(1, ndof);
        
        Bu_vec([1, 4]) = dNu_dx; % For du/dx
        Bw1_vec([2, 3, 5, 6]) = dNw_dx; % For dw/dx
        Bw2_vec([2, 3, 5, 6]) = d2Nw_dx2; % For d2w/dx2
        
        % Displacements and gradients at the current Gauss point
        du_dx = Bu_vec * d_elem;
        dw_dx = Bw1_vec * d_elem;
        
        % Stress resultants (Nxx, Mxx) at the Gauss point
        Nxx = E*A * (du_dx + 0.5 * dw_dx^2);
        Mxx = E*I * (Bw2_vec * d_elem); % Note: Mxx = EI*kappa = -EI*w''
                                       % For weak form, we need Mxx and Bw2
                                       % which leads to -Mxx and -Bw2, so signs cancel.
        
        % --- Assemble Internal Force Vector ---
        B_NL_u = Bu_vec + dw_dx * Bw1_vec;
        f_int = f_int + (B_NL_u' * Nxx + Bw2_vec' * Mxx) * J * w_gp;

        % --- Assemble Tangent Stiffness Matrix ---
        g = dw_dx;
        
        % Material part of stiffness
        KT_material = Bu_vec' * (E*A) * Bu_vec + Bw2_vec' * (E*I) * Bw2_vec;
        
        % Geometric part of stiffness
        KT_geometric = Bw1_vec' * Nxx * Bw1_vec;
        
        % Additional nonlinear terms for tangent matrix
        KT_nonlinear = Bw1_vec' * (E*A*g^2) * Bw1_vec + ...
                       Bu_vec' * (E*A*g) * Bw1_vec + ...
                       Bw1_vec' * (E*A*g) * Bu_vec;
                     
        K_T = K_T + (KT_material + KT_geometric + KT_nonlinear) * J * w_gp;
    end
end
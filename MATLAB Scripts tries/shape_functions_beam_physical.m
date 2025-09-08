function [Nu, Nw, dNu_dx, dNw_dx, d2Nw_dx2] = shape_functions_beam_physical(x, L)
% Provides shape functions and their derivatives in the PHYSICAL coordinate x
% for a 2-node beam element of length L.

    % --- Axial Shape Functions (Linear Lagrange) ---
    % For DOFs [u1, u2] at x=0 and x=L
    Nu = [ 1 - x/L, x/L ];
    dNu_dx = [ -1/L, 1/L ];
    
    % --- Transverse Shape Functions (Cubic Hermite) ---
    % For DOFs [w1, theta1, w2, theta2] where theta = dw/dx
    xL = x/L; % Normalized physical coordinate
    
    Nw = [ 1 - 3*xL^2 + 2*xL^3, ...        % Shape function for w1
           x * (1 - xL)^2, ...              % Shape function for theta1
           3*xL^2 - 2*xL^3, ...              % Shape function for w2
           x * (xL^2 - xL) ];                % Shape function for theta2
           
    % First derivatives w.r.t. x
    dNw_dx = [ (-6*x/L^2) + (6*x^2/L^3), ...
               (1 - xL)^2 - (2*x/L)*(1-xL), ...
               (6*x/L^2) - (6*x^2/L^3), ...
               (xL^2 - xL) + x*(2*x/L^2 - 1/L) ];
                
    % Second derivatives w.r.t. x
    d2Nw_dx2 = [ (-6/L^2) + (12*x/L^3), ...
                 (-2/L)*(1-xL) - (2/L)*(1-xL) - (2*x/L)*(-1/L), ...
                 (6/L^2) - (12*x/L^3), ...
                 (2*x/L^2 - 1/L) + (2*x/L^2 - 1/L) + x*(2/L^2) ];
end
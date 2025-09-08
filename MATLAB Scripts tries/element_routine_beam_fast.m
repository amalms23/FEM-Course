function [f_int, K_T] = element_routine_beam_fast(L, E, A, I, d_elem, f_int_handle, K_T_handle)
% Computes f_int and K_T by substituting numerical values into
% the pre-compiled symbolic function handles.

    % Unpack the displacement vector for the function handle
    u1 = d_elem(1); w1 = d_elem(2); th1 = d_elem(3);
    u2 = d_elem(4); w2 = d_elem(5); th2 = d_elem(6);

    % Call the fast, pre-compiled functions
    f_int = f_int_handle(L, E, A, I, u1, w1, th1, u2, w2, th2);
    K_T   = K_T_handle(L, E, A, I, u1, w1, th1, u2, w2, th2);
end
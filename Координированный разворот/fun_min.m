function f = fun_min(x)
% x(1) = k_omegay, x(2) = k_psi

    A_1 = 2800*x(2);
    A_0 = 142.8*x(3);

    f = A_1 * A_0^(-1);

end


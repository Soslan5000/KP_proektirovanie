function f = fun_min(x)
% x(1) = k_omegay, x(2) = k_psi

    A_1 = 1000*x(2) + 200*x(1) + 1216;
    A_0 = 200.0*x(2);

    f = A_1 * A_0^(-1);

end
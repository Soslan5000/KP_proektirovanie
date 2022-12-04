function [c, seq] = nonclon(x)
% x(1) = k_omegay, x(2) = k_psi

    A_5 = 1;
    A_4 = 28.6842;
    A_3 = 414.3537;
    A_2 = 1000*x(1) + 245.98;
    A_1 = 1000*x(2) + 200*x(1) + 1216;
    A_0 = 200.0*x(2);
  
    lambd_star = 2.15; % Желаемое значение лямбды
    lambd1 = lambd_star - (A_0^(-1) * A_1 * A_2 * A_3^(-1));
    lambd2 = lambd_star - (A_1^(-1) * A_2 * A_3 * A_4^(-1));
    lambd3 = lambd_star - (A_2^(-1) * A_3 * A_4 * A_5^(-1));
    
    sigm_star = 2; % Желаемое значение сигмы
    sigm1 = sigm_star - (A_0^(-1) * A_1^2 * A_2^(-1));
    sigm2 = sigm_star - (A_1^(-1) * A_2^2 * A_3^(-1));
    sigm3 = sigm_star - (A_2^(-1) * A_3^2 * A_4^(-1));

    c = [lambd1; lambd2; lambd3; sigm1; sigm2; sigm3];
    seq = [];

end
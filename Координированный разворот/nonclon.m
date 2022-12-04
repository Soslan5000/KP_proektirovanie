function [c, seq] = nonclon(x)
% x(1) = k_omegay, x(2) = k_psi

    A_5 = 1;
    A_4 = 29.28;
    A_3 = 428.28;
    A_2 = 400.0 + 2800.0*x(1);
    A_1 = 2800*x(2);
    A_0 = 142.8*x(3);
  
    lambd_star = 2.15; % Желаемое значение лямбды
    lambd1 = lambd_star - (A_0^(-1) * A_1 * A_2 * A_3^(-1));
    lambd2 = lambd_star - (A_1^(-1) * A_2 * A_3 * A_4^(-1));
    lambd3 = lambd_star - (A_2^(-1) * A_3 * A_4 * A_5^(-1));
    
    sigm_star = 2.1; % Желаемое значение сигмы
    sigm1 = sigm_star - (A_0^(-1) * A_1^2 * A_2^(-1));
    sigm2 = sigm_star - (A_1^(-1) * A_2^2 * A_3^(-1));
    sigm3 = sigm_star - (A_2^(-1) * A_3^2 * A_4^(-1));

    c = [lambd1; lambd2; lambd3; sigm1; sigm2; sigm3];
    seq = [];

end
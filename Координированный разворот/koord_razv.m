clc
clear

syms s k_psi_e k_omx k_sigm

% Зададим параметры в системе
H = 5;
M = 0.6;
V = 192;
Zb = -0.2;
sin = 0.08;
cos = 1;
g_del_V = 0.051;

M_x_b = -5.8;
M_x_omx = -1;
M_x_omy = -0.2;
M_x_sigme = -7;
M_y_b = -3;
M_y_omx = -0.05;
M_y_omy = -0.2;
M_y_sigmn = -2.5;

om_pr = 20;
kzi_pr = 0.707;

% Передаточная функция привода
W_pr = om_pr^2 / (s^2 + 2*om_pr*kzi_pr*s + om_pr^2);
% Передаточная функция ОУ
W1 = M_x_sigme / (s - M_x_omx);
% Передаточная функция ОУ и привода вместе
W2 = collect(W_pr*W1);
% Передаточная функция с обратной связью
W3 = collect(W2 / (1 - W2*k_omx));
% Умножаем на интегратор
W4 = collect(W3 / s);
% Ещё раз оборачиваем обратной связью
W5 = collect(W4 / (1 - W4*k_sigm));
% ПФ разомкнутой цепи
W6 = collect(k_psi_e * W5 * g_del_V / s);
% ПФ всей системы
W7 = collect(W6 / (1 - W6));

% Выделим числитель и знаменатель итоговой ПФ
[num, den] = numden(W7);
num = collect(num);
den = collect(den);
% Выделим коэффициенты числителя и знаменателя
% При этом сами коэффициенты будут записаны в порядке от младшего к
% старшему, то есть перевернутся
koefs_num = coeffs(num, 's');
koefs_den = coeffs(den, 's');
% Выделим старший коэффициент знаменателя. В списке он будет последний
starsh_koef = koefs_den(length(koefs_den));
% Поделим все коэффициенты на это число
koefs_num = vpa(koefs_num / starsh_koef);
koefs_den = vpa(koefs_den / starsh_koef);
% Развернём оба списка с коэффициентами, чтобы они шли от старшего ко
% младшему
koefs_num = rot90(rot90(koefs_num));
koefs_den = rot90(rot90(koefs_den));

% Преобразованные коэффициенты передаточной функции
B0 = koefs_num(1);
A5 = koefs_den(1);
A4 = koefs_den(2);
A3 = koefs_den(3);
A2 = koefs_den(4);
A1 = koefs_den(5);
A0 = koefs_den(6);

% Передаточная функция в нормальном виде (старший коэффициент знаменателя
% равен 1
num = B0;
den = (A5*s^5 + A4*s^4 + A3*s^3 + A2*s^2 + A1*s + A0);
PF = (num / den);
disp('Числитель ПФ')
pretty(num)
disp('Знаменатель ПФ')
pretty(den)

% Зададим параметры для оптимизации в fmicon
k0 = [1;1;1];
A = [];
b = [];
Aeg = [];
beg = [];
lb = 0.1*ones(2,1);
ub = 300*ones(2,1);

% Осуществим параметрическую оптимизацию
[x, fval] = fmincon('fun_min', k0, A, b, Aeg, beg, lb, ub, 'nonclon');

% Полученные коэффициенты
k_omx = x(1)
k_sigm = x(2)
k_psi_e = x(3)

% Преобразовывает коэффициенты характеристического полинома, 
% подставляя в них найденные значения
B_0 = vpa(-142.8*x(3));
A_5 = vpa(1);
A_4 = vpa(29.28);
A_3 = vpa(428.28);
A_2 = vpa(400.0 + 2800.0*x(1));
A_1 = vpa(2800.0*x(2));
A_0 = vpa(142.8*x(3));
num = B_0;
den = A_5*s^5 + A_4*s^4 + A_3*s^3 + A_2*s^2 + A_1*s + A_0;
num = sym2poly(num);
den = sym2poly(den);

% Получаем передаточную функцию и строим переходный процесс
disp('Передаточная функция с полученными коэффициентами')
W_sys = tf(num, den)
figure;
step(-1*W_sys)
grid on
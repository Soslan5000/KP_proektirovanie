clc
clear

syms k_psi k_omegay s

% Определение параметров системы
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
kzi_pr = sqrt(2)/2;

% Передаточная функция привода
W_pr = om_pr^2 / (s^2 + 2*om_pr*kzi_pr*s + om_pr^2);
% Передаточная функция по угловой скорости рысканья с учётом привода
W_sigmn_omy_raz = collect((M_y_sigmn*s + (-M_y_sigmn*Zb)) / (s^2 + (-M_y_omy-Zb)*s + (M_y_omy*Zb-M_y_b)) * W_pr);
% Обратная связь с коэффициентом обратной связи k_omegay
W_sigmn_omy_zamk = collect(W_sigmn_omy_raz / (1 - k_omegay*W_sigmn_omy_raz));
% Разомкнутая передаточная функция по углу рысканья с замкнутым внутренним
% контуром
W_psizad_psi_raz = collect(k_psi * W_sigmn_omy_zamk / s);
% Передаточная функция замкнутой системы по углу рысканья
W_psizad_psi_zamk = collect(W_psizad_psi_raz / (1 - W_psizad_psi_raz));

% Выделим числитель и знаменатель итоговой ПФ
[num, den] = numden(W_psizad_psi_zamk);
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
B1 = koefs_num(1);
B0 = koefs_num(2);
A5 = koefs_den(1);
A4 = koefs_den(2);
A3 = koefs_den(3);
A2 = koefs_den(4);
A1 = koefs_den(5);
A0 = koefs_den(6);

% Передаточная функция в нормальном виде (старший коэффициент знаменателя
% равен 1
num = (B1*s + B0);
den = (A5*s^5 + A4*s^4 + A3*s^3 + A2*s^2 + A1*s + A0);
PF = (num / den);
disp('Числитель ПФ')
pretty(num)
disp('Знаменатель ПФ')
pretty(den)

% Зададим параметры для оптимизации в fmicon
k0 = [10, 10];
A = [];
b = [];
Aeg = [];
beg = [];
lb = 0.1*ones(2,1);
ub = 300*ones(2,1);

% Осуществим параметрическую оптимизацию
[x, fval] = fmincon('fun_min', k0, A, b, Aeg, beg, lb, ub, 'nonclon');

% Полученные коэффициенты
k_omegay = x(1)
k_psi = x(2)

% Преобразовывает коэффициенты характеристического полинома, 
% подставляя в них найденные значения
B_1 = vpa(-1000.0*x(2));
B_0 = vpa(-200.0*x(2));
A_5 = vpa(1);
A_4 = vpa(28.6842);
A_3 = vpa(414.3537);
A_2 = vpa(1000*x(1) + 245.98);
A_1 = vpa(1000*x(2) + 200*x(1) + 1216);
A_0 = vpa(200.0.*x(2));
num = B_1*s + B_0;
den = A_5*s^5 + A_4*s^4 + A_3*s^3 + A_2*s^2 + A_1*s + A_0;
num = sym2poly(num);
den = sym2poly(den);

% Получаем передаточную функцию и строим переходный процесс
disp('Передаточная функция с полученными коэффициентами')
W_sys = tf(num, den)
figure;
step(-1*W_sys)
grid on
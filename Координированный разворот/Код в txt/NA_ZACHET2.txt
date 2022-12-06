clc
clear

syms k_psi_e k_omx k_sigm s

% Зададим параметры в системе
H = 5;
M = 0.6;
V = 192;
Zb = 0.2;
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

% Запишем неравенства для достаточного условия устойчивости
lambd_star = 2.15; % Желаемое значение лямбды
eq1d = A0^(-1) * A1 * A2 * A3^(-1) >= lambd_star;
eq2d = A1^(-1) * A2 * A3 * A4^(-1) >= lambd_star;
eq3d = A2^(-1) * A3 * A4 * A5^(-1) >= lambd_star;


% Преобразуем неравенства из символьного типа данных в функции
eq1d = matlabFunction(eq1d);
eq2d = matlabFunction(eq2d);
eq3d = matlabFunction(eq3d);

% Задание сетки
ko = [-1:0.01:4];
ks = [-1:0.01:15];
kp = [30:50:300];
[KO, KS] = meshgrid(ko, ks);
figure;
hold on;
grid on;
colors = {'green', 'blue', 'yellow', 'blue', "magenta", 'red', 'black'};
ind_col = 1;
for i=1:length(kp)
    eq1 = eq1d(k_omx, kp(i), k_sigm);
    eq2 = eq2d(k_omx, k_sigm);
    eq3 = eq3d(k_omx);
    eq1 = matlabFunction(eq1);
    eq2 = matlabFunction(eq2);
    eq3 = matlabFunction(eq3);
    c1 = eq1(KO, KS);
    c2 = eq2(KO, KS);
    c3 = eq3(KO);
    col = colors(ind_col);
    col = string(col);
    ind_col = ind_col + 1;
    contourf(ko, ks, c1 & c2 & c3, [1 1 1], FaceAlpha=0.2, FaceColor=col)
    colormap lines
end
hold off
xlabel('Komx')
ylabel('Ksigm')
title('Область устойчивости по достаточным условиям для kpsie=30, 80, 130, 180, 230, 280')

% Запишем неравенства для достаточного условия устойчивости
eq1r = A0 > 0;
eq2r = A1 > 0;
eq3r = A2 > 0;
eq4r = A3 > 0;
eq5r = A4 > 0;
eq6r = A5 > 0;
eq7r = simplify((A1*A2-A0*A3)*(A3*A4-A2*A5) - (A1*A4-A0*A5)^2) > 0;


% Преобразуем неравенства из символьного типа данных в функции
eq1r = matlabFunction(eq1r);
eq2r = matlabFunction(eq2r);
eq3r = matlabFunction(eq3r);
eq4r = matlabFunction(eq4r);
eq5r = matlabFunction(eq5r);
eq6r = matlabFunction(eq6r);
eq7r = matlabFunction(eq7r);

% Задание сетки
ko = [-1:0.01:4];
ks = [-1:0.01:18];
kp = [30:50:300];
[KO, KS] = meshgrid(ko, ks);
figure;
hold on;
grid on;
colors = {'green', 'blue', 'yellow', "magenta", 'red', 'black'};
ind_col = 1;
for i=1:length(kp)
    eq2 = eq2r(k_sigm);
    eq3 = eq3r(k_omx);
    eq7 = eq7r(k_omx, kp(i), k_sigm);
    eq2 = matlabFunction(eq2);
    eq3 = matlabFunction(eq3);
    eq7 = matlabFunction(eq7);
    c1 = eq2(KS);
    c2 = eq3(KO);
    c3 = eq7(KO, KS);
    col = colors(ind_col);
    col = string(col);
    ind_col = ind_col + 1;
    contourf(ko, ks, c1 & c2 & c3, [1 1 1], FaceAlpha=0.2, FaceColor=col)
end
hold off
xlabel('Komx')
ylabel('Ksigm')
title('Область устойчивости по критерию Рауса для kpsie=30, 80, 130, 180, 230, 280')

% Рабочая точка
komx = 0.9225;
ksigm = 3.5337;
kpsie = 109.4343;

figure;
grid on;
hold on;
ko = [-1:0.01:4];
ks = [-1:0.01:15];
[KO, KS] = meshgrid(ko, ks);
eq1 = eq1d(k_omx, kpsie, k_sigm);
eq2 = eq2d(k_omx, k_sigm);
eq3 = eq3d(k_omx);
eq1 = matlabFunction(eq1);
eq2 = matlabFunction(eq2);
eq3 = matlabFunction(eq3);
c1 = eq1(KO, KS);
c2 = eq2(KO, KS);
c3 = eq3(KO);
contourf(ko, ks, c1 & c2 & c3, [1 1 1], FaceAlpha=0.2, FaceColor='red')

ko = [-1:0.01:4];
ks = [-1:0.01:18];
kp = [30:50:300];
[KO, KS] = meshgrid(ko, ks);
eq2 = eq2r(k_sigm);
eq3 = eq3r(k_omx);
eq7 = eq7r(k_omx, kpsie, k_sigm);
eq2 = matlabFunction(eq2);
eq3 = matlabFunction(eq3);
eq7 = matlabFunction(eq7);
c1 = eq2(KS);
c2 = eq3(KO);
c3 = eq7(KO, KS);
contourf(ko, ks, c1 & c2 & c3, [1 1 1], FaceAlpha=0.2, FaceColor='green')

scatter(komx, ksigm, 'filled')
hold off
xlabel('Komx')
ylabel('Ksigm')
title('Область устойчивости по двум критериям для kpsie=109.4343')
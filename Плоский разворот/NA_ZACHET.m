clc
clear

syms k_psi k_omegay s

% Задаём параметры системы
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


% Запишем неравенства для достаточных условий устойчивости через лямбды
lambd_star = 2.15;
eq1d = A0^(-1) * A1 * A2 * A3^(-1) >= lambd_star;
eq2d = A1^(-1) * A2 * A3 * A4^(-1) >= lambd_star;
eq3d = A2^(-1) * A3 * A4 * A5^(-1) >= lambd_star;


% Преобразуем неравенства из символьного типа данных в функции
eq1d = matlabFunction(eq1d);
eq2d = matlabFunction(eq2d);
eq3d = matlabFunction(eq3d);

% Запишем неравенства для критерия Рауса
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
l1 = [-2:0.01:20];
l2 = [-2:0.01:45];
[ko, kp] = meshgrid(l1, l2);


% Проверка устойчивости узлов сетки на достаточные условия устойчивости
c1d = eq1d(kp, ko);
c2d = eq2d(kp, ko);
c3d = eq3d(ko);

% Проверка устойчивости узлов сетки по критерию Рауса
c1r = eq1r(kp);
c2r = eq2r(kp, ko);
c3r = eq3r(ko);
c4r = eq4r();
c5r = eq5r();
c6r = eq6r();
c7r = eq7r(kp, ko);

figure;
grid on
hold on
% Контур для области устойчивости по критерию Рауса
contourf(l1, l2, c1r & c2r & c3r & c4r & c5r & c6r & c7r, [1 1 1 1 1 1 1], FaceAlpha=0.5, FaceColor="red")
% Контур для области устойчивости по достаточным условиям устойчивости
contourf(l1, l2, c1d & c2d & c3d, [1 1 1],  FaceAlpha=0.5, FaceColor="blue")
title('Красная область - достаточные условия устойчивости. Синяя область - критерий Рауса')

% Синтезированные коэффициенты
k_omegay = 2.7468;
k_psi = 9.0425;
pnt1 = scatter(k_omegay, k_psi,'r','filled');
% Подобранные коэффициенты
k_omegay = 2.8172;
k_psi = 9.0425;
pnt2 = scatter(k_omegay, k_psi, 'y','filled');
legend([pnt1, pnt2],"Синтезированные коэффициенты", "Подобранные коэффициенты")
hold off
xlabel("Komy")
ylabel("Kpsi")
grid on
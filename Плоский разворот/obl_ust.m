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

% Передаточная функция в нормальном виде (старший коэффициент знаменателя
% равен 1
num = (B1*s + B0);
den = (A5*s^5 + A4*s^4 + A3*s^3 + A2*s^2 + A1*s + A0);
PF = (num / den);
disp('Числитель ПФ')
pretty(num)
disp('Знаменатель ПФ')
pretty(den)

% Запишем неравенства ждя критерия Рауса
eq1 = A0 > 0;
eq2 = A1 > 0;
eq3 = A2 > 0;
eq4 = A3 > 0;
eq5 = A4 > 0;
eq6 = A5 > 0;
eq7 = simplify((A1*A2-A0*A3)*(A3*A4-A2*A5) - (A1*A4-A0*A5)^2) > 0;

% Преобразуем неравенства из символьного типа данных в функции
eq1 = matlabFunction(eq1);
eq2 = matlabFunction(eq2);
eq3 = matlabFunction(eq3);
eq4 = matlabFunction(eq4);
eq5 = matlabFunction(eq5);
eq6 = matlabFunction(eq6);
eq7 = matlabFunction(eq7);

% Задание сетки
l1 = [-2:0.05:50];
l2 = [-2:0.05:50];
[ko, kp] = meshgrid(l1, l2);

% Проверка устойчивости узлов сетки
c1 = eq1(kp);
c2 = eq2(kp, ko);
c3 = eq3(ko);
c4 = eq4();
c5 = eq5();
c6 = eq6();
c7 = eq7(kp, ko);

% Построение области устойчивости
figure;
contourf(l1, l2, c1 & c2 & c3 & c4 & c5 & c6 & c7, [1 1 1 1 1 1 1])
colormap lines
hold on
% Синтезированные коэффициенты
k_omegay = 2.7468;
k_psi = 9.0425;
pnt1 = scatter(k_omegay, k_psi,'r','filled');
% Подобранные коэффициенты
k_omegay = 2.8172;
k_psi = 9.0425;
pnt2 = scatter(k_omegay, k_psi, 'b','filled');
legend([pnt1, pnt2],"Синтезированные коэффициенты", "Подобранные коэффициенты")
hold off
xlabel("Komy")
ylabel("Kpsi")
grid on
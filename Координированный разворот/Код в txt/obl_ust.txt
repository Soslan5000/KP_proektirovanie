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
ko = [-1:0.05:4];
ks = [-1:0.05:20];
kp = [-1:0.5:300];

% Построение трехмерного графика
% Задаём k_omega
% Строим плоскую область для параметров k_sigm, k_psi
% Прибавляет к параметру k_omega шаг 0.05
% Повторяем действия
% Таким образом - получится много плоских рисунков,
% которые в сумме будут давать пространство устойчивости
hold on
grid on
view(3)
for i=1:length(ko)
    X = [];
    Y = [];
    Z = [];
    ind = 1;
    for j=1:length(ks)
        for k=1:length(kp)
            c1 = eq1(kp(k));
            c2 = eq2(ks(j));
            c3 = eq3(ko(i));
            c4 = eq4();
            c5 = eq5();
            c6 = eq6();
            c7 = eq7(ko(i), kp(k), ks(j));
            if  c1 & c2 & c3 & c4 & c5 & c6 & c7
                X(ind) = ko(i);
                Y(ind) = ks(j);
                Z(ind) = kp(k);
                ind = ind + 1;
            end
        end
    end
    plot3(X, Y, Z)
end
xlabel('Komx')
ylabel('Ksigm')
zlabel('Kpsie')
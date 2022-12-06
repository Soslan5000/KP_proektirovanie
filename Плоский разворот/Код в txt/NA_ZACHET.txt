clc
clear

syms k_sigm k_omx s

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


W_pr = om_pr^2 / (s^2 + 2*om_pr*kzi_pr*s + om_pr^2); % ПФ привода
W1 = collect(M_x_sigme / (s - M_x_omx)); % ПФ по угловой скорости
W2 = collect(W1 * W_pr); % Связываем W1 и W2
W3 = collect(W2 / (1 - W2*k_omx)); % Положительная ОС по угловой скорости
W4 = collect(k_sigm * W3 / s); % Разомкнутая цепь по углу крена
W5 = collect(W4 / (1 - W4)); % Замкнутая цепь по углу крена
disp('ПФ замкнутой цепи по углу крена с учётом привода')
pretty(W5)


[num, den] = numden(W5);
disp('Числитель передаточной функции')
num = vpa(collect(num), 4) 
disp('Знаменатель передаточной функции')
den = vpa(collect(den), 4)


% Выделим коэффициенты в числителе и знаменателе
% Функция coeffs запишет их от младшей степени к старшей
koefs_num = coeffs(num, 's'); % Массив коэффициентов числителя
koefs_den = coeffs(den, 's'); % Массив коэффициентов знаменателя


% Перевернём массивы на 180 градусов, чтобы получить коэффициенты 
% от старшего ко младшему
koefs_num = rot90(rot90(koefs_num));
koefs_den = rot90(rot90(koefs_den));


% Запишем соответствующие коэффициенты в соответствующие переменные
B0 = koefs_num(1);
A4 = koefs_den(1);
A3 = koefs_den(2);
A2 = koefs_den(3);
A1 = koefs_den(4);
A0 = koefs_den(5);


% Запишем неравенства для достаточных условий устойчивости через лямбды
lambd_star = 2.15;
eq1d = A0^(-1) * A1 * A2 * A3^(-1) >= lambd_star;
eq2d = A1^(-1) * A2 * A3 * A4^(-1) >= lambd_star;


% Преобразуем неравенства из символьного типа данных в функции
eq1d = matlabFunction(eq1d);
eq2d = matlabFunction(eq2d);

% Запишем неравенства жля критерия Рауса
eq1r = A0 > 0;
eq2r = A1 > 0;
eq3r = A2 > 0;
eq4r = A3 > 0;
eq5r = A4 > 0;
eq6r = simplify(A1*A2*A3 - A1^2*A4 - A0*A3^2) > 0;

% Преобразуем неравенства из символьного типа данных в функции
eq1r = matlabFunction(eq1r);
eq2r = matlabFunction(eq2r);
eq3r = matlabFunction(eq3r);
eq4r = matlabFunction(eq4r);
eq5r = matlabFunction(eq5r);
eq6r = matlabFunction(eq6r);


% Задание сетки
l1 = [-2:0.005:5];
l2 = [-2:0.005:20];
[ko, ks] = meshgrid(l1, l2);


% Проверка устойчивости узлов сетки на достаточные условия устойчивости
c1d = eq1d(ko, ks);
c2d = eq2d(ko);

% Проверка устойчивости узлов сетки по критерию Рауса
c1r = eq1r(ks);
c2r = eq2r(ko);
c3r = eq3r();
c4r = eq4r();
c5r = eq5r();
c6r = eq6r(ko, ks);

figure;
grid on
hold on
% Контур для области устойчивости по критерию Рауса
contourf(l1, l2, c1r & c2r & c3r & c4r & c5r & c6r, [1 1 1 1 1 1], FaceAlpha=0.5, FaceColor="red")
% Контур для области устойчивости по достаточным условиям устойчивости
contourf(l1, l2, c1d & c2d, [1 1],  FaceAlpha=0.5, FaceColor="blue")
title('Красная область - достаточные условия устойчивости. Синяя область - критерий Рауса')

% Синтезированные коэффициенты
k_omx = 1;
k_sigm = 3.6957;
pnt1 = scatter(k_omx, k_sigm, 'y','filled');
% Подобранные коэффициенты
k_omx = 1.2597;
k_sigm = 5.4231;
pnt2 = scatter(k_omx, k_sigm, 'b','filled');
legend([pnt1, pnt2],"Синтезированные коэффициенты", "Подобранные коэффициенты")
xlabel('Komx')
ylabel('Ksigm')
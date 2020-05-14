%{
 { Вариант с учетом ограничений.
 { Задается некоторое макс. вращение w_max, и поиск ведется с его учетом.
 { Начальная точка (пара коэф. управления) для процедуры fminsearch
 { была получена в ф-ции get_good_start, как обеспечивающая минимальный
 { максимум на некотором дискретном промежутке коэффициентов.
 %}
global zs;
global m;
global b;
global g;
global w;
global t0;
global w_max;

zs = 20;
m = 2;
b = 3;
g = 9.8;
%{
 { w is 0.989949493661167
 %}
w = sqrt((m*g)/(4*b));
t0 = 4*b*w*w;

global i;
i = 0;

main1()

%{
 { function [] = main()
 {     global w_max;
 { 
 {     w_max = 3.0;
 {     k_p_d = fminsearch(@(k)(calc_integral(k)),[1 1]);
 {     kpx = k_p_d(1)
 {     kdx = k_p_d(2)
 { 
 {     lt = 2.3;
 {     rt = 3.0;
 {     while lt < rt
 {         mid = (lt + rt) * 0.5;
 {         w_max = mid
 {         k_p_d = fminsearch(@(k)(calc_integral(k)),[1 1]);
 {         kp = k_p_d(1)
 {         kd = k_p_d(2)
 { 
 {         if ((kp == 1) & (kd == 1))
 {             lt = mid;
 {         else 
 {             if ((kp == kpx) & (kd == kdx))
 {                 rt = mid;
 {             else
 {                 break
 {             end
 {         end
 {     end
 { 
 {     [z, dz] = get_solution_handle(kp, kd);
 {     t = 0:001:50;
 {     figure()
 {     plot(t, arrayfun(z, t))
 { end
 %}

function [] = main1()
	global w_max;
	global b;
	global w;

	w0 = w
	w_max = 1.3 
	%{
	 { k_p_d = fminsearch(@(k)(calc_integral(k)),[1 1]);
	 %}
	k_p_d = fminsearch(@(k)(calc_integral(k)),[0.1 9.1]);

	kp = k_p_d(1)
	kd = k_p_d(2)
	integr = calc_integral([kp kd])
	if integr == 100000000
		fprintf('solution not found\n')
	end

	[z, dz] = get_solution_handle(kp, kd);
	t = 0:001:50;
	figure()
	plot(t, arrayfun(z, t))

%{
 {     u = get_u(z, dz, kp, kd);
 {     t = 0:0.0001:50;
 {     y = u(t);
 {     [val, idx] = max(y);
 { 
 {     w_max = sqrt(val / (4 * b))
 %}
end

function [] = get_good_start()
	mn = 1000000000;
	kpa = 0;
	kda = 0;
	for kp = 0.1:1:10
		for kd = 0.1:1:10
			[z, dz] = get_solution_handle(kp, kd);
			u = get_u(z, dz, kp, kd);
			t = 0:0.0001:50;
			y = u(t);
			[val, idx] = max(y);
			if val < mn
				mn = val;
				kpa = kp;
				kpd = kd;
			end
		end
	end
	kpa
	kpd
end

function res = calc_integral(k)
	global i;
	global w_max;
	global b;

	kp = k(1);
	kd = k(2);

	%{
	 {  i = i + 1;
	 {  if (i < 8)
	 {     [z, dz] = get_solution_handle(kp, kd);
	 {     t = 0:1:50;
	 {     figure()
	 {     plot(t, arrayfun(z, t))
	 { end
	 %}

	[z, dz] = get_solution_handle(kp, kd);
	u = get_u(z, dz, kp, kd);
	 t = 0:0.0001:50;
	 y = u(t);
	 [val, idx] = max(y);
	 if val > 4 * b * w_max * w_max
		 res = 100000000;
		 return
	end

	%{
	 { [val, idx] = min(y);
	 { if val < 0
	 {      res = 100000000;
	 {      return
	 { end
	 %}

	%{
	 { xx1 = 1
	 %}

	integral_handle = get_integral_handle(kp, kd);
	res = integral(integral_handle, 0, 50,'ArrayValued', true) + penalty(kp, kd);
end

function handle = get_integral_handle(kp, kd)
	global zs;
	[z, dz] = get_solution_handle(kp, kd);
	u = get_u(z, dz, kp, kd);
	handle = @(t)double((z(t) - zs) ^ 2 + dz(t) ^ 2 + u(t) ^ 2);
end

function handle = get_u(z, dz, kp, kd)
	global t0;
	global zs;
	handle = @(t)(kp * (z(t) - zs) + kd * dz(t) + t0);
end

function [z_h, dz_h] = get_solution_handle(kp, kd)
	global m;
	global zs;
	syms z(t);
	eqn = m * diff(z,t,2) == -kp * (z - zs) - kd * diff(z,t,1);
	Dz = diff(z, t);

	cond = [z(0)==0, Dz(0)==0];
	zSol(t) = dsolve(eqn, cond);
	z_h = matlabFunction(zSol);
	dz_h = matlabFunction(diff(zSol));
end

function res = penalty(kp, kd)
	if (kp <= 0 || kd <= 0)
		res = 100000000;
	else
		res = 0;
	end
end

%{
 { Программа вывод коэф. процесса управления, первые несколько графиков
 { оптимизацинного процесса, и последний график для минимуме функционала
 %}
global zs;
global m;
global b;
global g;
global w;
global t0;

zs = 20;
m = 2;
b = 1;
g = 9.8;
w = sqrt((m*g)/(4*b));
t0 = 4*b*w*w;

global i;
i = 0;

main();

function [] = main()
	k_p_d = fminsearch(@(k)(calc_integral(k)),[1 1]);
	kp = k_p_d(1)
	kd = k_p_d(2)

	[z, dz] = get_solution_handle(kp, kd);
	t = 0:001:50;
	figure()
	plot(t, arrayfun(z, t))
end

function res = calc_integral(k)
	global i;
	kp = k(1);
	kd = k(2);

	 i = i + 1;
	 if (i < 8)
		[z, dz] = get_solution_handle(kp, kd);
		t = 0:1:50;
		figure()
		plot(t, arrayfun(z, t))
	end

	integral_handle = get_integral_handle(kp, kd);
	res = integral(integral_handle, 0, 50,'ArrayValued', true) + penalty(kp, kd);
end

function handle = get_integral_handle(kp, kd)
	global zs;
	[z, dz] = get_solution_handle(kp, kd);
	u = get_u(z, dz, kp, kd);
	handle = @(t)double((z(t) - zs) ^ 2 + dz(t) ^ 2 + u(t));
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

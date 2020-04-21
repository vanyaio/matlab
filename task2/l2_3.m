%{
 { Разбор данного модельного ур-ния приведен в рекомендуемых материалах,
 { основные выводы оттуда, что явный метод имеет малую область устойчивости, которая
 { обеспечивается step <= 2/abs(-10000), а неявный А-устойчив, т.е. устойчив
 { при всех Re(lamba) < 0
 %}
global h;
global steps;
steps = 900000;

h = 0.00009;
mainreal3();

h = 0.00019;
mainreal3();

h = 0.00021;
mainreal3();

function [] = mainreal3()
	global h;
	global steps;

	cond = 100000;
	[x, y] = euler_for3(cond, 0, h, steps);
	figure;
	plot(x, y(1, :), '-o');
	name = sprintf('Forward Euler, step=%f', h);
	legend(name);

	[x, y] = euler_back_for3(cond, 0, h, steps);
	figure;
	plot(x, y(1, :), '-o');
	name = sprintf('Backward Euler, step=%f', h);
	legend(name);
end

function [x, y] = euler_for3(y0, x0, h, steps)
	x = [];
	y = [];

	y(1) = y0;
	x(1) = x0;
	for j = 2:steps
		y(j) = (1 - h * 10000) * y(j-1);
		x(j) = x(j-1) + h;
		if (j < 10)
			y(j)
		end
	end
end
function [x, y] = euler_back_for3(y0, x0, h, steps)
	x = [];
	y = [];

	y(1) = y0;
	x(1) = x0;
	for j = 2:steps
		y(j) = (1 / (1 + h * 10000)) * y(j-1);
		x(j) = x(j-1) + h;
	end
end

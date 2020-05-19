%{
 { 1. Из подкоренного выражения в правой части виднo, что ОДУ не удовлетворяет
 { условию теоремы Пеано o непрерывности при b > arcsin(1) = pi/2
 %}
%{
 { 2.Для данного уравнения выполняется теорема существования и
 { единственности решения задачи Коши, но локально, т.е. в
 { некоторой окрестности. Соответсвенно для данных начальных условий
 { tan(x) есть решение, но оно является таковым только при x < pi/2.
 %}
%{
 { 3. Для дифф. ур-ние не выполняется условие Липшица
 { (док-во: запишем его; перенесем слагаемые в одну сторону; заметим
 { неограниченность Ly1 - y1^1/3 и зафиксируем вторую переменную).
 { Соотвественно здесь нет единственности решения.
 %}
	
main1();

function [x, y] = euler(y0, x0, f, h, steps)
	x = [];
	y = [];
	n = length(y0);
	for j = 1:steps
		for i = 1:n
			if (j == 1)
				y(i,1) = y0(i);
				x(1) = x0;
				x(2) = x0 + h;
				continue
			end
			y(i,j) = y(i,j-1) + h * f{i}(x(j-1), y(:, j-1));
			x(j) = x(j-1) + h;
		end
	end
end

function [x, y] = runge(y0, x0, f, h, steps)
	x = [];
	y = [];
	n = length(y0);

	for j = 1:steps

		if (j == 1)
			x(j) = x0;
			y(:,j) = y0;
			continue
		end

		%{
		 { k_i is row, y_prev is row
		 %}
		k1 = [];
		k2 = [];
		k3 = [];
		k4 = [];
		for i = 1:n
			curf = f{i};
			y_prev = transpose(y(:, j-1));

			k1(i) = curf(x(j-1), y_prev);
		end
		for i = 1:n
			curf = f{i};
			y_prev = transpose(y(:, j-1));

			k2(i) = curf(x(j-1) + 0.5*h, y_prev + k1 * (0.5*h));
		end
		for i = 1:n
			curf = f{i};
			y_prev = transpose(y(:, j-1));

			k3(i) = curf(x(j-1) + 0.5*h, y_prev + k2 * (0.5*h));
		end
		for i = 1:n
			curf = f{i};
			y_prev = transpose(y(:, j-1));

			k4(i) = curf(x(j-1) + h, y_prev + k3 * h);
		end
		k1 = transpose(k1);
		k2 = transpose(k2);
		k3 = transpose(k3);
		k4 = transpose(k4);
		y(:, j) = y(:, j-1) + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
		x(j) = x(j-1) + h;
	end
end


function res = p_1(x, y)
	res = sqrt(1 - (y(1)^2));
end
function res = p_2(x, y)
	res =  y(1) ^ 2 + 1;
end
function res = p_3(x, y)
	res =  y(1) ^ (1/3);
end

function [] = main1()
	h = 0.0001;

	for i=1:3
		f = cell(1);
		if i == 1
			f{1} = @p_1;
			b = pi;
		end
		if i == 2
			f{1} = @p_2;
			b = pi/2;
		end
		if i == 3
			h = 0.01;
			b = 3;
			f{1} = @p_3;

			syms y(t)
			eqn = diff(y,t) == y^(1/3);
			Dy = diff(y,t);

			cond = [y(0)==0];
			ySol(t) = dsolve(eqn,cond);
		end

		steps = ceil(b / h);
		[x, y] = euler([0], 0, f, h, steps);
		[x1, y1] = runge([0], 0, f, h, steps);
		subplot(3,1,i)
		hold on
		plot(x, y(1, :), 'green');
		hold on
		plot(x1, y1(1, :), 'red');
		if i == 1
			hold on
			plot(x, sin(x), 'blue');
			legend('euler', 'runge', 'sin')
		end
		if i == 2
			hold on
			plot(x, tan(x), 'blue');
			legend('euler', 'runge', 'tan')
		end
		if i == 3
			f = @(x,y) y^(1/3);
			x = 0:h:b;
			[x2,y2] = ode45(f, x, 0);
			hold on
			plot(x2, y2, 'blue')

			hold on
			y3 = [];
			i = 1;
			for p = x
				ans = ySol(p);
				y3(i) = ans(1);
				i = i + 1;
			end
			plot(x, y3, 'black')
			legend('euler', 'runge', 'ode', 'syms')
		end
	end
end


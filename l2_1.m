eq1();
eq2();

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

function [x, y] = euler_back(y0, x0, f, h, steps)
	x = [];
	y = [];
	n = length(y0);

	for j = 1:steps

		if (j ~= 1)
			for i = 1:n
				ys(i) = y(i,j-1) + h * f{i}(x(j-1), y(:, j-1));
			end
		end

		for i = 1:n
			if (j == 1)
				y(i,1) = y0(i);
				x(1) = x0;
				x(2) = x0 + h;
				continue
			end
			f1 = f{i}(x(j-1), y(:, j-1));
			x(j) = x(j-1) + h;
			f2 = f{i}(x(j), ys);
			y(i,j) = y(i,j-1) + h * (f1 + f2) * 0.5;
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
	res = y(2);
end
function res = p_2(x, y)
	l = 16.81;
	res = l * y(1);
end
function [] = eq1()
	syms y(t)
	eqn = diff(y,t,2) == 16.81 * y;
	Dy = diff(y,t);

	ySol(t) = dsolve(eqn)

	cond = [y(0)==1, Dy(0)==-4.1];
	ySol(t) = dsolve(eqn,cond);

	h = 0.01;
	steps = 1000;
	x = 0:h:h*steps;
	figure;
	plot(x, arrayfun(matlabFunction(ySol), x), '-o');
	legend('1st eq, analytic solution');

	f = cell(2,1);
	f{1} = @p_1;
	f{2} = @p_2;

	[x, y] = euler([1 -4.1], 0, f, h, steps);
	figure;
	plot(x, y(1, :), '-o');
	legend('1st eq, euler solution');

	[x, y] = runge([1 -4.1], 0, f, h, steps);
	figure;
	plot(x, y(1, :), '-o');
	legend('1st eq, runge solution');
end

function [] = eq2()
	syms y(t)
	eqn = diff(y,t,2) == -16.81 * y - 8.2 * diff(y,t);
	Dy = diff(y,t);

	ySol2(t) = dsolve(eqn)

	cond = [y(0)==1, Dy(0)==-4.1];
	ySol2(t) = dsolve(eqn,cond);

	h = 0.01;
	steps = 1000;

	x = 0:h:h*steps;
	figure;
	plot(x, arrayfun(matlabFunction(ySol2), x), '-o');
	legend('2st eq, analytic solution');

	f = cell(2,1);
	f{1} = @g_1;
	f{2} = @g_2;

	[x, y] = euler([1 -4.1], 0, f, h, steps);
	figure;
	plot(x, y(1, :), '-o');
	legend('2st eq, euler solution');

	[x, y] = runge([1 -4.1], 0, f, h, steps);
	figure;
	plot(x, y(1, :), '-o');
	legend('2st eq, runge solution');
end
function res = g_1(x, y)
	res = y(2);
end
function res = g_2(x, y)
	a = -16.81;
	b = -8.2;
	res = a * y(1) + b * y(2);
end

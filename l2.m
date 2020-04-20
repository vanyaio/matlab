main()

function [x, y] = euler(y0, x0, f, h, steps)
	x = []
	y = []
	n = length(y0)
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

function res = f_1(x, y)
	y
	y(2)
	res = y(2);
end
function res = f_2(x, y)
	res = -y(1);
end

function [] = main()
	syms y(t)
	eqn = diff(y,t,2) == -y;
	Dy = diff(y,t);

	ySol(t) = dsolve(eqn);

	cond = [y(0)==0, Dy(0)==1];
	ySol(t) = dsolve(eqn,cond);

	h = 0.1;
	i = 1;
	y = [];
	for x = 0:h:2 * pi
		y(i) = ySol(x);
		i = i + 1;
	end

	figure;
	plot(0:h:2 * pi, y, '-o');

	f = cell(2,1);
	f{1} = @f_1;
	f{2} = @f_2;
	h = 0.01;
	steps = 700;
	[x, y] = euler_back([0 1], 0, f, h, steps);
	%{
	 { [x, y] = euler([0 1], 0, f, h, steps);
	 %}

	figure;
	plot(x, y(1, :), '-o');
end

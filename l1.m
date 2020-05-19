1. Из подкоренного выражения в правой части видно не неприрывность
при b > arcsin(1) = pi/2
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

function res = p_1(x, y)
	res = sqrt(1 - (y(1)^2));
end
function res = p_2(x, y)
	res =  y(1) ^ 2 + 1;
end
function res = p_3(x, y)
	res =  y(1) ^ (1.0/3.0);
end

function [] = main1()
	h = 0.001;
	steps = 100;
	for i=1:3
		f = cell(1);
		if i == 1
			f{1} = @p_1;
			%{
			 { b = pi/3;
			 %}
		end
		if i == 2
			f{1} = @p_2;
		end
		if i == 3
			f{1} = @p_3;
		end

		steps = ceil(b / h);
		[x, y] = euler([0], 0, f, h, steps);
		[x1, y1] = runge([0], 0, f, h, steps);
		hold on
		subplot(3,1,i)
		plot(x, y(1, :), 'green');
		hold on
		plot(x1, y1(1, :), 'red');
		if i == 1
			hold on
			plot(x, sin(x), 'blue');
			legend('euler', 'runge', 'sin')
			%{
			 { legend('euler', 'runge')
			 %}
		else
			legend('euler', 'runge')
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


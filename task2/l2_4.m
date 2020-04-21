x0 = 0.5;
dx0 = 10;

%{
 { hold on;
 %}
figure();
analytic(x0, dx0, 1);
%{
 { odesol(x0, dx0, 1);
 %}

figure();
f = cell(2,1);
f{1} = @g_1;
f{2} = @g_2;
h = 0.0001;
steps = round(1 / h)
[x, y] = euler([x0 dx0], 0, f, h, steps);
plot(x, y(1, :), '-');


figure();
stiff_sol(x0, dx0, 1);


%{
 { legend('analytic')
 %}
%{
 { legend('analytic', 'euler')
 %}
%{
 { legend('analytic', 'euler', 'stiff');
 %}


function res = g_1(x, y)
	res = y(2);
end
function res = g_2(x, y)
	a = -100000.001;
	b = -100;
	res = a * y(2) + b * y(1);
end

function sol = analytic(x0, dx0, end_t)
	syms x(t)
	eqn = diff(x,t,2) == -100 * x - 100000.001 * diff(x,t,1)
	Dx = diff(x,t);

	cond = [x(0)==x0, Dx(0)==dx0];
	sol(t) = dsolve(eqn,cond)

	x = [];
	i = 1;
	h = 0.001;
	for t = 0:h:end_t
		x(i) = sol(t);
		i = i + 1;
	end

	plot(0:h:end_t, x, '-')
end

function [] = odesol(x0, dx0, end_t)
	syms y(t)
	x_end = end_t;

	[V] = odeToVectorField(diff(y,2) == -100 * y - 100000.001 * diff(y))
	M = matlabFunction(V,'vars', {'t','Y'})
	sol = ode45(M,[0 x_end],[x0 dx0]);

	fplot(@(x)deval(sol,x,1), [0, x_end])
end

function [] = stiff_sol(x0, dx0, end_t)
	syms y(t)
	x_end = end_t;

	[V] = odeToVectorField(diff(y,2) == -100 * y - 100000.001 * diff(y))
	M = matlabFunction(V,'vars', {'t','Y'})
	sol = ode15s(M,[0 x_end],[x0 dx0]);

	fplot(@(x)deval(sol,x,1), [0, x_end])
end



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

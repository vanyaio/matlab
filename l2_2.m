%{
 { main_21();
 %}
main_22();

function [] = main_21
	t2_1(0.5);
	t2_3(0.5);
	t2_1(1.5);
	t2_3(1.5);
	t2_1(2);
	t2_3(2);
	t2_1(3);
	t2_3(3);
	t2_1(5);
	t2_3(5);
end
function [] = main_22
	t2_2(0.5);
	t2_4(0.5);
	t2_2(1.5);
	t2_4(1.5);
	t2_2(2);
	t2_4(2);
	t2_2(3);
	t2_4(3);
	t2_2(5);
	t2_4(5);
end

function [] = t2_1(v)
	syms theta(t) theta_t(t);
	eqs = [diff(theta)   == theta_t;
		   diff(theta_t) == -sin(theta)];

	vars = [theta, theta_t];
	[M,F] = massMatrixForm(eqs,vars);
	f = M\F;
	f = odeFunction(f, vars);

	x0 = [0; v];
	tInit  = 0;
	tFinal = 7;
	ols = ode45(f,[tInit tFinal],x0);

	name = sprintf('nonlinear no-friction, init speed=%d', v);
	figure('Name', name);

	yyaxis left;
	plot(ols.x, ols.y(1,:), '-o');
	ylabel('\theta (rad)');
	xlabel('time');
end

function [] = t2_2(v)
	syms y(t)
	x_end = 20

	[V] = odeToVectorField(diff(y, 2) == -0.5 * diff(y) - sin(y))
	M = matlabFunction(V,'vars', {'t','Y'})
	sol = ode45(M,[0 x_end],[0 v]);

	name = sprintf('nonlinear with-friction, init speed=%d', v);
	figure('Name', name);

	fplot(@(x)deval(sol,x,1), [0, x_end])
end

function [] = t2_3(v)
	syms y(t)
	eqn = diff(y,t,2) == -y;
	Dy = diff(y,t);

	cond = [y(0)==0, Dy(0)==v];
	ySol(t) = dsolve(eqn,cond);

	y = [];
	i = 1;
	h = 0.3;
	for x = 0:h:7
		y(i) = ySol(x);
		i = i + 1;
	end

	name = sprintf('linear no-friction, init speed=%d', v);
	figure('Name', name);

	plot(0:h:7, y, '-o');
end

function [] = t2_4(v)
	syms y(t)
	eqn = diff(y,t,2) == -y - 0.5 * diff(y,t,1);
	Dy = diff(y,t);

	cond = [y(0)==0, Dy(0)==v];
	ySol(t) = dsolve(eqn,cond);

	y = [];
	i = 1;
	h = 0.3;
	for x = 0:h:6*pi
		y(i) = ySol(x);
		i = i + 1;
	end

	name = sprintf('linear with-friction, speed=%d', v);
	figure('Name', name);

	plot(0:h:6 * pi, y, '-o');
end

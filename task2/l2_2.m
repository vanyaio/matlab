%{
 { В моделях без трения при малых начальных скоростях лин. и
 { нелин. вариант ведут схожим образом (углы малы - лин.вариант адекватен),
 { при увеличении начальной скорости постепенно различия
 { нарастают (сравните скорости 0.5 и 1.5),
 { но сам вид кривой весьма схож - однако это свойство теряется
 { при верхнем положении устойчивости (скорость=2), когда нелинейный
 { вариант перестает иметь синусообразный вид, который всегда остается
 { при лин. варианте просто из-за соответвующего явного решения. Отметим,
 { что верхняя точка численно неустойчива и все равно маятник из-за
 { погрешности вычислений начинает сваливаться. В моделях с трением
 { справедливы аналогичные выводы про их схожесть в зависимости от нач.
 { скоростей и соответственно малости углов, но вдобаков отметим, что
 { лин. вариант со скоростями (например 5), которые дают прохождение
 { верхней точки в нелин. модели, этот перелет не обнаруживают.
 %}

main_21();
main_22();

function [] = main_21
	figure(); hold on;
	v = 0.5
	t2_1(v);
	t2_3(v);
	name = sprintf('nonlinear no-friction, init speed=%f', v);
	name1 = sprintf('linear no-friction, init speed=%f', v);
	legend(name, name1);

	figure(); hold on;
	v = 1.5
	t2_1(v);
	t2_3(v);
	name = sprintf('nonlinear no-friction, init speed=%f', v);
	name1 = sprintf('linear no-friction, init speed=%f', v);
	legend(name, name1);

	figure(); hold on;
	v = 2
	t2_1(v);
	t2_3(v);
	name = sprintf('nonlinear no-friction, init speed=%f', v);
	name1 = sprintf('linear no-friction, init speed=%f', v);
	legend(name, name1);

	figure(); hold on;
	v = 3
	t2_1(v);
	t2_3(v);
	name = sprintf('nonlinear no-friction, init speed=%f', v);
	name1 = sprintf('linear no-friction, init speed=%f', v);
	legend(name, name1);

	figure(); hold on;
	v = 5
	t2_1(v);
	t2_3(v);
	name = sprintf('nonlinear no-friction, init speed=%f', v);
	name1 = sprintf('linear no-friction, init speed=%f', v);
	legend(name, name1);
end
function [] = main_22
	figure(); hold on;
	v = 0.5
	t2_2(v);
	t2_4(v);
	name = sprintf('nonlinear with-friction, init speed=%f', v);
	name1 = sprintf('linear with-friction, init speed=%f', v);
	legend(name, name1);

	figure(); hold on;
	v = 1.5
	t2_2(v);
	t2_4(v);
	name = sprintf('nonlinear with-friction, init speed=%f', v);
	name1 = sprintf('linear with-friction, init speed=%f', v);
	legend(name, name1);

	figure(); hold on;
	v = 2
	t2_2(v);
	t2_4(v);
	name = sprintf('nonlinear with-friction, init speed=%f', v);
	name1 = sprintf('linear with-friction, init speed=%f', v);
	legend(name, name1);

	figure(); hold on;
	v = 3
	t2_2(v);
	t2_4(v);
	name = sprintf('nonlinear with-friction, init speed=%f', v);
	name1 = sprintf('linear with-friction, init speed=%f', v);
	legend(name, name1);

	figure(); hold on;
	v = 5
	t2_2(v);
	t2_4(v);
	name = sprintf('nonlinear with-friction, init speed=%f', v);
	name1 = sprintf('linear with-friction, init speed=%f', v);
	legend(name, name1);
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

	yyaxis left;
	plot(ols.x, ols.y(1,:), 'g');
end

function [] = t2_2(v)
	syms y(t)
	x_end = 20

	[V] = odeToVectorField(diff(y, 2) == -0.5 * diff(y) - sin(y))
	M = matlabFunction(V,'vars', {'t','Y'})
	sol = ode45(M,[0 x_end],[0 v]);


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


	plot(0:h:7, y, 'r');
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


	plot(0:h:6 * pi, y, 'r');
end

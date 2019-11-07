%xn is row

%solve_meth
main()

function x_val = sys(xn, fun, x)
	n = numel(xn)

	A = ones(n, 1)
	for i = 2:n
		A(:,i) = transpose(xn) .* A(:, i - 1)
	end

	b = arrayfun(fun, transpose(xn))
	coeff = linsolve(A, b)
	
	x_val = []
	for ptr = x
		powers = [1]
		for i = 2:n
			powers = [powers; powers(i - 1) * ptr]
		end
		x_val = [x_val, sum(powers .* coeff)]
	end
end

function x_val = lagr(xn, fun, x)
	x_val = []
	for ptr = x
		res = 0
		for node1 = xn
			prod1 = 1
			prod2 = 1
			for node2 = xn
				if (node1 == node2)
					continue
				end
				prod1 = prod1 * (ptr - node2)
				prod2 = prod2 * (node1 - node2)
			end
			res = res + fun(node1) *  (prod1 / prod2)
		end
		x_val = [x_val res]
	end
end

%ptr_meth
function xn = cheb(a, b, n)
	for i = 0:n-1
		xn(1, i + 1) = 0.5 * ((b - a) * cos( ((2*i + 1) / (2*n)) * pi ) + a + b)
	end
end

function err = get_error(a, b, xn, solve_meth, fun)
	x = a:1:b
	true_val = arrayfun(fun, x)
	calc_val = solve_meth(xn, fun, x)
	err = max(abs(true_val - calc_val))
end

function ret = analys(a, b, ptr_meth, solve_meth, fun)
	ret = []
	eps = 0.001
	for n = 25 
		xn = ptr_meth(a, b, n)
		ret(n, 1) = n
		ret(n, 2) = get_error(a, b, xn, solve_meth, fun)
		if (ret(n, 2) < eps)
			break
		end
	end
end

function [] = main_fun(a, b, fun)
	cheb_sys = analys(a, b, @cheb, @sys, fun)
	cheb_lagr = analys(a, b, @cheb, @lagr, fun)

	lin_lagr = analys(a, b, @linspace, @lagr, fun)
	lin_sys = analys(a, b, @linspace, @sys, fun)
end

function [] = main()
	a = - 5 * pi
	b = 5 * pi
	main_fun(a, b, @sin)

	main_fun(-1, 1, @sec_fun)
end

function res = sec_fun(x)
	res = 1 / (1 + 12 * x * x)
end

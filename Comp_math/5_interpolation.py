#This programm create interpolate Newton's polinomial and cubic spline for set of points in file in.txt
import math
import sympy
import numpy as np
import matplotlib.pyplot as plt

def f(list):
	if len(list) > 1:
		return (f(list[1:]) - f(list[:-1]))/(list[-1][0] - list[0][0])
	else:
		return list[0][1]

#Function creates interpolation polinomial P
def newton(list):
	t = sympy.symbols('t')
	P, s, j = 0, 1, 0
	list_arg = []
	for i in list:
		j += 1
		P = P + f(list[:j]) * s
		s = s * (t - i[0])
	return (P,t)

#Function for coefficients P_(N, i) in cubic spline
def get_coef(P, x):
	ans = []
	P_x = sympy.diff(P, sympy.symbols('t'))
	for i in x:
		ans.append(P_x.subs(sympy.symbols('t'), i))
	return ans

#Function create cubic spline on segments in list with polinimial P(argument P_X is his derivative)
def cubic_spline(P, list):
	x, f, ans = [], [], []
	for i in list:
		x.append(i[0])
		f.append(i[1])
	C = get_coef(P, x)
	t = sympy.Symbol('t')
	for i in range(len(list) - 1):
		a3 =    ((C[i+1] * (x[i+1] - x[i]) -
			2 * (f[i+1] - f[i]) +
			C[i] * (x[i+1] - x[i])) /
			(x[i+1] - x[i])**3)
		a2 =    ((-C[i+1] * (x[i+1] - x[i]) * (x[i+1] + 2 * x[i])+
			3 * (f[i+1] - f[i]) * (x[i+1] + x[i]) -
			C[i] * (x[i+1] - x[i]) * (x[i] + 2 * x[i+1])) /
			(x[i+1] - x[i])**3)
		a1 =    ((C[i+1] * x[i] * (2* x[i+1] + x[i]) * (x[i+1] - x[i]) -
			6 * (f[i+1] - f[i]) * x[i] * x[i+1] +
			C[i] * x[i+1] * (x[i+1] + 2 * x[i]) * (x[i+1] - x[i]))/
			(x[i+1] - x[i])**3)
		a0 =    ((-C[i+1] * x[i]**2 * x[i+1] * (x[i+1] - x[i]) +
			f[i+1] * x[i]**2 * (3 * x[i+1] - x[i])+
			f[i] * x[i+1]**2 * (x[i+1] - 3 * x[i]) -
			C[i] * x[i] * x[i+1]**2 * (x[i+1] - x[i])) /
			(x[i+1] - x[i])**3)
		ans.append((a3 * t**3 + a2 * t**2 + a1 * t + a0, [x[i], x[i+1]]))
	return ans

def calc_spline(spl, x):
	for i in spl:
		if x >= i[1][0] and x <= i[1][1]:
			return i[0].subs(sympy.symbols('t'), x)
	return 0

def make_graph(size, f, s, e, n, num, name):
	x = []
	for i in np.linspace(s, e, n):
		x.append([i, f(i)])
	P = newton(x)
	spline3 = cubic_spline(P[0], x)
	
	arg = np.linspace(s - 1, e + 1, 100)
	plt.subplot(2, size, num)
	plt.plot([i[0] for i in x], [i[1] for i in x], 'go')
	plt.plot(arg, [f(i) for i in arg], 'b')
	plt.plot(arg, [P[0].subs(P[1], i) for i in arg], 'r')
	plt.title(r'Newtons polinomial for %s'%(name))
	plt.legend(['Interpolted points', r'Original function', 'Newtons polinomial'])
	plt.xlabel(r'x')
	plt.ylabel(r'y')
	plt.grid(True)
	
	arg = np.linspace(s, e, 100)
	value = [calc_spline(spline3, float(i)) for i in arg]
	plt.subplot(2, size, num + size)
	plt.plot([i[0] for i in x], [i[1] for i in x], 'go')
	plt.plot(arg, [f(i) for i in arg], 'b')
	plt.plot(arg, value, 'r')
	plt.legend(['Interpolted points', r'Original function', 'Cubic spline'])
	plt.xlabel(r'x')
	plt.ylabel(r'y')
	plt.title(r'Cubic spline for %s'%(name))
	plt.grid(True)


f1 = lambda x: math.exp(x)
f2 = lambda x: (x - 0.4) * (x - 0.5) * (x - 0.3) * (x - 0.7) * (x - 0.12)
f3 = lambda x: math.sin(x * math.pi)

make_graph(3, f1, -3, 3, 5, 1, r'exponent $e^x$, 5 points')
make_graph(3, f2, -10, 10, 5, 2, 'polinomial of 5 degree, 5 points')
make_graph(3, f3, -3, 3, 10, 3, r'sinus $\sin(\pi x)$,10 points')

plt.subplots_adjust(left=0.04, bottom=0.05, right=0.99, top=0.97, wspace=0.15, hspace=0.20)
plt.show()

import math
import sympy


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
	return P

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

file = open('var_1', 'r')
x = []
for line in file:
	print(line)
	x.append([float(i) for i in line.split()])

P = newton(x).expand()
print(P)
spline3 = cubic_spline(P, x)
for i in spline3:
	print(i[1],': P(t)=',i[0])
print('Input x to get value of cubic spline in point x')
print('Print E to end programm')
while (1):
	inp = input()
	if inp == 'E':
		break
	print(calc_spline(spline3, float(inp)))

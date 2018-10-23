#This programm include two parts:
#1. Show one disadvantages of Newton's Method
#2. Show how many iterations needs for to find all roots of a random polynomial of constant degree
#3. Show how work dichotomy for different degrees of polinomial with same number of iteration
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.misc import derivative

#Function for creating vector (1, x, x^2 ... x^n)^T
def vector_pow(x, N):
	ans = np.zeros((N,))
	for i in range(N):
		ans[i] = x**i
	return ans

#This function give knowledge about convergence of Newton's Method
def newton_convergence(f, x, stop):
	x_0, x_prev = x, 0
	N = 0
	f_val, x_val = [], []
	while N < stop:
		f_val.append(abs(f(x_0)))
		x_val.append(abs(x_prev - x_0))
		N += 1
		df = derivative(f, x_0, dx = 0.0000001)
		if df != 0:
			x_prev, x_0 = x_0, x_0 - (df)**(-1) * f(x_0)
		else:
			return (f_val, x_val)
	return (f_val, x_val)


#Arguments -function that's roots must be found, interval, number of subinterval
#Return number of found roots
def dichotomy(f, interval, N):
	s, e = interval
	step = (e - s)/N
	num_of_roots = 0
	f_1, f_2 = f(s), f(s + step)
	for i in range(N):
		if f_1 * f_2 < 0:
			num_of_roots += 1
		f_1, f_2 = f_2, f(s + (i + 2) * step)
	return num_of_roots


#The first part: newton's method
#Number of iterations for Newton's Method
NI = 15
#f(x) = x^3 - 2x + 2
#This example will be in cycle with start point x = 0
a = [2, -2, 0, 1]
list_x = [i * 10 for i in range(0, 3, 1)]
colors = ['b', 'r', 'y']
legend = []
j = 0
for i in list_x:
	#polinomial with coefficient's vector a
	poly = lambda x: np.dot(a, vector_pow(x, 4))
	q_bad_poly = newton_convergence(poly, i, NI)
	plt.subplot(2, 2, 1)
	plt.plot(range(NI), q_bad_poly[0], colors[j])
	plt.subplot(2, 2, 2)
	plt.plot(range(NI), q_bad_poly[1], colors[j])
	legend.append('Start points %s'%(i))
	j+=1
print("Newton Method was done success")
plt.subplot(2, 2, 1)
plt.grid(True)
plt.xlabel(r'Number of iteration $n$')
plt.ylabel(r'$f(x_n) = x_n^3-2x_n+2$')
plt.legend(legend)
plt.title(r'Convergence of Newton Method in bad case')
plt.subplot(2, 2, 2)
plt.grid(True)
plt.xlabel(r'Number of iteration $n$')
plt.ylabel(r'$x_n - x_{n-1}$')
plt.legend(legend)
plt.title(r'Stabilization of Newton Method bad case')

#The second part: dichotomy
degree = 4
num, num_of_roots = [], []
for i in range(1, 101, 2):
	NR = []
	for j in range(500):
		roots = np.random.uniform(-1, 1, (degree,))
		f = lambda x: np.prod(np.array([x for i in range(degree)]) - roots)
		N = dichotomy(f, (-1, 1), i)
		NR.append(N)
	num_of_roots.append(np.array(NR).mean())
	if int(i - 1) % 10 == 0:
		print("Dichotomy was done success: %d %%"%(i-1))
print("Dichotomy was done success: 100 %%")

plt.subplot(2, 2, 3)
plt.plot(range(1, 101, 2), num_of_roots, 'r')
plt.grid(True)
plt.xlabel(r'Number of interval')
plt.ylabel(r'Num of roots')
plt.title(r'Dychotomy')


#The third part: different degrees
degree = 4
num, num_of_roots = [], [0,]
i = 100
for degree in range(1, 51, 1):
	NR = []
	for j in range(100):
		roots = np.random.uniform(-1, 1, (degree * 10,))
		f = lambda x: np.prod(np.array([x for i in range(degree * 10)]) - roots)
		N = dichotomy(f, (-1, 1), i)
		NR.append(N)
	num_of_roots.append(np.array(NR).mean())
	if int(degree) % 5 == 0:
		print("Dichotomy for different degrees of polinomial was done success: %d %%"%(degree * 2))

plt.subplot(2, 2, 4)
plt.plot([i * 10 for i in range(0, 51, 1)], num_of_roots, 'r')
plt.grid(True)
plt.xlabel(r'Degree of polinomial')
plt.ylabel(r'Num of roots')
plt.title(r'100 intervals for several degree')
plt.show()

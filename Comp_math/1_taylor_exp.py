import math
from matplotlib.pylab import plot, savefig, show, grid
import matplotlib.pyplot as plt
import numpy as np

def Delta(u_find, u):
	return abs(u_find - u)

a = float(input())
tuple = np.linspace(start = -a, stop = a, num = a * 10)
exp_x, exp = [], []

for x in tuple:
	sum, term, i = 1.0, 1.0, 1.0
	
	exp.append(math.exp(x))
	while term != 0:
		term_prev, term = term, term * x / i
		sum += term
		i+= 1

	if term_prev < 0 and x < 0:
		sum = sum - term_prev

	exp_x.append(sum)

plot(tuple, exp_x)
grid(True)
plt.xlabel(r'$x$')
plt.ylabel(r'$\exp(x)$')
a = [Delta(exp_x[i], exp[i]) for i in range(len(exp))]
plt.title(r'Exponent, maximum of Delta %s'%(max(a)))


print('Maximum of delta ',max(a))
show()

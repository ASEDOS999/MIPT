import math
import random
import numpy as np
import matplotlib.pylab as plt
def Monte_Carlo_cube(f, N, num_arg):
	x = np.random.uniform(0, 1, (N, num_arg))
	sum = 0
	for i in range(N):
		sum += f(x[i])
	return sum / N

def rectangle(f, N, start, end):
	step = (end - start) / N
	sum = 0
	for i in range(N):
		sum += f(start + (i + 1/2) * step) * step
	return sum

def trapezium(f, N, start, end):
	step = (end - start) / N
	sum = 0
	for i in range(N):
		sum += step * (f(start + i * step) + f(start + (i + 1) * step)) / 2
	return sum

def sympson(f, N, start, end):
	step = (end - start) / N
	sum = 0
	for i in range(N):
		sum += step * (f(start + i * step) + 4 * f (start + (i + 1/ 2) * step)+ f(start + (i + 1) * step)) / 6
	return sum


size_test = 100
f_test_1 = lambda x: math.exp(x)
f_test_2 = lambda x: math.pi * math.sin(math.pi * x) - 0.5
f_test_3 = lambda x: 3 * math.exp(3 * x)
TI = [(math.e - 1), 1.5, (math.e**3 - 1)]
rec_value1, rec_value2, rec_value3 = [], [], []
tr_value1, tr_value2, tr_value3 = [], [], []
sym_value1, sym_value2, sym_value3 = [], [], []
for i in range(size_test):
	rec_value1.append(math.log(abs(rectangle(f_test_1, i + 1, 0, 1) - TI[0]), 10))
	rec_value2.append(math.log(abs(rectangle(f_test_2, i + 1, 0, 1) - TI[1]), 10))
	rec_value3.append(math.log(abs(rectangle(f_test_3, i + 1, 0, 1) - TI[2]), 10))
	sym_value1.append(math.log(abs(sympson(f_test_1, i + 1, 0, 1) - TI[0]), 10))
	sym_value2.append(math.log(abs(sympson(f_test_2, i + 1, 0, 1) - TI[1]), 10))
	sym_value3.append(math.log(abs(sympson(f_test_3, i + 1, 0, 1) - TI[2]), 10))
	tr_value1.append(math.log(abs(trapezium(f_test_1, i + 1, 0, 1) - TI[0]), 10))
	tr_value2.append(math.log(abs(trapezium(f_test_2, i + 1, 0, 1) - TI[1]), 10))
	tr_value3.append(math.log(abs(trapezium(f_test_3, i + 1, 0, 1) - TI[2]), 10))
plt.subplot(2, 3, 1)
plt.plot(range(1, size_test + 1, 1), rec_value1, 'r')
plt.plot(range(1, size_test + 1, 1), tr_value1, 'b')
plt.legend(['Rectangle', 'Trapezium'])
plt.xlabel('Number of parts')
plt.ylabel('Differences')
plt.grid(True)
plt.title(r'$e^x$')
plt.subplot(2, 3, 2)
plt.plot(range(1, size_test + 1, 1), rec_value2, 'r')
plt.plot(range(1, size_test + 1, 1), tr_value2, 'b')
plt.legend(['Rectangle', 'Trapezium'])
plt.xlabel('Number of parts')
plt.ylabel('Differences')
plt.title(r'$\pi \sin(\pi x)$')
plt.grid(True)
plt.subplot(2, 3, 3)
plt.plot(range(1, size_test + 1, 1), rec_value3, 'r')
plt.plot(range(1, size_test + 1, 1), tr_value3, 'b')
plt.legend(['Rectangle', 'Trapezium'])
plt.xlabel('Number of parts')
plt.ylabel('Differences')
plt.title(r'$3e^{3x}$')
plt.grid(True)

plt.subplot(2, 3, 4)
plt.plot(range(1, size_test + 1, 1), sym_value1, 'r')
plt.plot(range(1, size_test + 1, 1), tr_value1, 'b')
plt.legend(['Simpson', 'Trapezium'])
plt.xlabel('Number of parts')
plt.ylabel('Differences')
plt.grid(True)
plt.title(r'$e^x$')
plt.subplot(2, 3, 5)
plt.plot(range(1, size_test + 1, 1), sym_value2, 'r')
plt.plot(range(1, size_test + 1, 1), tr_value2, 'b')
plt.legend(['Simpson', 'Trapezium'])
plt.xlabel('Number of parts')
plt.ylabel('Differences')
plt.title(r'$\pi \sin(\pi x)$')
plt.grid(True)
plt.subplot(2, 3, 6)
plt.plot(range(1, size_test + 1, 1), sym_value3, 'r')
plt.plot(range(1, size_test + 1, 1), tr_value3, 'b')
plt.legend(['Simpson', 'Trapezium'])
plt.xlabel('Number of parts')
plt.ylabel('Differences')
plt.title(r'$3e^{3x}$')
plt.grid(True)
plt.subplots_adjust(left=0.04, bottom=0.05, right=0.99, top=0.97, wspace=0.15, hspace=0.20)
plt.show()

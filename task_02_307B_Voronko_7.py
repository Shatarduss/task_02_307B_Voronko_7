from scipy.special import spherical_jn as jn
from scipy.special import spherical_yn as yn
import matplotlib.pyplot as plt
import urllib.request as url
from numpy import pi
from re import split
import numpy as np
import os

# определение составляющих ряда
def hn(l, z):
    return jn(l, z) + 1j * yn(l, z)

def an(l, z):
    return jn(l, z) / hn(l, z)

def bn(l, z):
    return (z * jn(l - 1, z) - l * jn(l, z)) \
           / (z * hn(l - 1, z) - l * hn(l, z))

# определение исходных данных
URL = 'https://jenyay.net/uploads/Student/Modelling/task_02.txt'
file = url.urlopen(URL)
list = file.readlines()
my_string = list[6].decode("utf-8") # переводим bytes в string
values = split("[=\r;]", my_string) # выделяем значения
D = float(values[1])
fmin = float(values[3])
fmax = float(values[5])

# рассчет ЭПР
N = 1000    # число точек на отрезке f
r = 0.5 * D
f = np.linspace(fmin, fmax, N)  # создаем решетку
Lambda = 3e8 / f
k = 2 * pi / Lambda # волновое число
Sum_arr = [((-1) ** n) * (n + 0.5) * (an(n, k * r) - bn(n, k * r)) \
           for n in range(1, 50)]
Sum = np.sum(Sum_arr, axis=0)   # считаем "ряд"
Sigma = (Lambda ** 2) / pi * (np.abs(Sum) ** 2) # находим ЭПР

# построение графика
plt.plot(f/10e6, Sigma)
plt.xlabel('$f, МГц$')
plt.ylabel('$\sigma, м^2$')
plt.grid()
plt.show()

# создаем папку выхода с результатом
try:
    os.mkdir('results')
except OSError:
    pass
complete_file = os.path.join('results', 'task_02_307B_Voronko_7.txt')
ftl = f.tolist() # преобразуем nparray в list, иначе ошибка
Stl = Sigma.tolist()
f = open(complete_file, 'w')
f.write('f    Sigma\n')
for i in range(N):
    f.write(str(ftl[i])+'    '+str(Stl[i])+"\n")
f.close()
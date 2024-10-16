#Численные методы задание №1.1 Интерполяция

import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
import numdifftools as nd
np.random.seed(42)

#построение полинома Лагранжа
def lagrange_interpolation(x_values, y_values):
    x = sp.symbols('x')
    def create_l_i(x_values, i): #подсчет базовых полиномов
        def l_i(x):
            l_i = 1
            for j in range(len(x_values)):
                if j != i:
                    l_i *= (x - x_values[j])/(x_values[i] - x_values[j])
            return l_i
        return l_i

    l_i = []
    for i in range(len(x_values)):
        l_i.append(create_l_i(x_values, i))

    def lagrange_polynomial(x):
        result = 0
        for i in range(len(y_values)):
            result += y_values[i]* l_i[i](x)
        return result
      # Формируем символьный полином
    polynomial_expr = sum(y_values[i] * create_l_i(x_values, i)(x) for i in range(len(y_values)))
    
    # Округляем коэффициенты полинома
    rounded_coefficients = [round(coef.evalf(), 2) for coef in polynomial_expr.as_coefficients_dict().values()]
    
    # Создаем новый полином с округленными коэффициентами
    rounded_polynomial_expr = sum(coef * x**i for i, coef in enumerate(rounded_coefficients))

    return lagrange_polynomial, rounded_polynomial_expr

#опредление функции
def f(x):
    return 1/3 + np.cos(10 + 2.3 ** np.abs(x))

a = -1
b = 2.5
k = 4 #степень полинома (Степень интерполяционного многочлена на единицу меньше числа узлов.)
xi = np.random.uniform(a, b) #выбираем рандомную точку на нашем промежутке. Благодаря строчке №7 при перезапуске кода не изменяется
#определение узлов
x_nodes = np.linspace(a, b, k + 1) #k+1 кол-во узлов
y_nodes = f(x_nodes)
lagrange_poly_func, lagrange_poly_expr = lagrange_interpolation(x_nodes, y_nodes) #Обозначаем полином Лагрнажа и его вид (№1)

#распределение точек
x_values = np.linspace(a, b, 500)
y_values = lagrange_poly_func(x_values)#значения полинома в точках

# Вычисляем производную с помощью numdifftools и находим его значение в выбранной точке xi
f_diff = nd.Derivative(f, n = 5)
f_diff_xi = f_diff(xi)

# R_k(x) = {f'(xi)/(n+1)!}*произведение по i от 0 до к(xi - x_nodes)- погрешность сверху в узлах 
product_term = 1
R_kk = []
for i in range(k + 1):
    product_term *= (xi - x_nodes[i])
    R_kk. append(abs((f_diff_xi/sp.factorial(k+1))*product_term))
R_k = np.min(abs((f_diff_xi/sp.factorial(k+1))*product_term)) #оценка сверху, так как остаточный член равен супремумм, то берем мин(№2)


y_f = f(x_values)#значения изначальной функции в точках
difference = abs(y_f - y_values)
error_fact = np.max(difference)#фактическая погрешность (№3)

print(lagrange_poly_expr)
print(f'Погрешность сверху: {round(R_k, 2)}, при значении: {round(xi, 2)}')
print(f'Фактическая ошибка: {round(error_fact, 2)}')

plt.figure(figsize=(12, 8))
plt.plot(x_values, y_f, label='f(x)', color='red', linestyle='-')
plt.plot(x_values, y_values,label='Полином Лагранжа', color='green', linestyle='-' )
plt.xlabel('x')
plt.ylabel('y')
plt.axhline(0, color='black', linewidth=0.5, ls='--')  # Горизонтальная линия на уровне y=0
plt.grid(color='gray', linestyle='--', linewidth=0.5)
plt.legend()

# Показать график
plt.show()
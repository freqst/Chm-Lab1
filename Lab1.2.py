#Численные методы задание №1.2 Интерполяция

import numpy as np
import sympy as sp
import numdifftools as nd
import matplotlib.pyplot as plt


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

a = -1
b = 2.5

def f(x):
    return 1/3 + np.cos(10 + 2.3 ** abs(x))

def chebyshev_roots(n):
    # Вычисление корней по формуле
    k = np.arange(1, n + 1)
    roots = np.cos((2 * k - 1) * np.pi / (2 * n))
    
    return roots

#roots т.к степень полинома Лагранжа равна n - 1, где n - степень полинома Чебышева, 
#для построения полинома нужной нам степени мы повышаем степень полинома Чебышева при построении узлов
degree = [3, 4, 5, 6]
x_roots = []
y_roots = []
lagrange_funcs = []
lagrange_expr = []
x_plot = np.linspace(-1, 2.5, 100)
max_errors = []
plt.figure(figsize=(12, 6))
for n in degree:
    x_n = chebyshev_roots(n + 1)
    y_n = f(x_n)
    x_roots.append(x_n)
    y_roots.append(y_n)
    lagrange_poly_func, lagrange_poly_exprr = lagrange_interpolation(x_n, y_n)
    lagrange_funcs.append(lagrange_poly_func) 
    lagrange_expr.append(lagrange_poly_exprr) # Сохраняем функцию в список
    # Вычисляем погрешность интерполяции
    error = abs(f(x_plot) - lagrange_poly_func(x_plot))
    plt.subplot(1, 2, 1)
    plt.plot(x_plot, error, label=f'Ошибка интерполяции (n={n})')
    # Нахождение максимальной погрешности
    max_error = np.max(np.abs(error))
    max_errors.append(max_error)


# График фактической погрешности интерполяции
plt.subplot(1, 2, 1)
plt.axhline(0, color='black', lw=0.5, ls='--')
plt.title('Фактическая погрешность интерполяции')
plt.xlabel('x')
plt.ylabel('Погрешность')
plt.legend()
plt.grid()

# График максимальной погрешности в зависимости от степени полинома
plt.subplot(1, 2, 2)
plt.plot(degree, max_errors, marker='o')
plt.title('Максимальная погрешность интерполяции')
plt.xlabel('Степень полинома')
plt.ylabel('Максимальная погрешность')
plt.grid()
plt.xticks(degree)

# Показываем графики
plt.tight_layout()
plt.show()

# Печать максимальных значений погрешности
for n, max_error in zip(degree, max_errors):
    print(f"Максимальная погрешность интерполяции для n={n}: {max_error:.4f}")
for n, lagrange_expr in zip(degree, lagrange_expr):
    print(f'Полином степени n = {n}: {lagrange_expr}')
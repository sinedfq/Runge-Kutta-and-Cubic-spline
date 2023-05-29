import math

# Класс кубического сплайна
class CubicSpline:
    def __init__(self, x_values, y_values):
        self.xData = x_values
        self.yData = y_values
        n = len(self.xData)

        self.h = [self.xData[i + 1] - self.xData[i] for i in range(n - 1)]  # Разности между значениями x
        self.alpha = [0.0] * (n - 1)  # Список для хранения промежуточных значений alpha
        self.l = [0.0] * n  # Список для хранения промежуточных значений l
        self.u = [0.0] * (n - 1)  # Список для хранения промежуточных значений u
        self.z = [0.0] * n  # Список для хранения промежуточных значений z
        self.c = [0.0] * n  # Список для хранения коэффициентов c
        self.b = [0.0] * (n - 1)  # Список для хранения коэффициентов b
        self.d = [0.0] * (n - 1)  # Список для хранения коэффициентов d

        for i in range(1, n - 1):
            self.alpha[i] = 3.0 * ((self.yData[i + 1] - self.yData[i]) / self.h[i] - (self.yData[i] - self.yData[i - 1]) / self.h[i - 1])
            #  Формула естесвенного сплайна
        self.l[0] = 1.0
        self.u[0] = 0.0
        self.z[0] = 0.0

        for i in range(1, n - 1):
            self.l[i] = 2.0 * (self.xData[i + 1] - self.xData[i - 1]) - self.h[i - 1] * self.u[i - 1]
            self.u[i] = self.h[i] / self.l[i]
            self.z[i] = (self.alpha[i] - self.h[i - 1] * self.z[i - 1]) / self.l[i]

        self.l[n - 1] = 1.0
        self.z[n - 1] = 0.0
        self.c[n - 1] = 0.0

        for j in range(n - 2, -1, -1):
            self.c[j] = self.z[j] - self.u[j] * self.c[j + 1]
            self.b[j] = (self.yData[j + 1] - self.yData[j]) / self.h[j] - self.h[j] * (self.c[j + 1] + 2.0 * self.c[j]) / 3.0
            self.d[j] = (self.c[j + 1] - self.c[j]) / (3.0 * self.h[j])

    # Функция интерполяции
    def interpolate(self, x):
        n = len(self.xData)
        index = 0

        # Находим интервал для интерполяции
        for i in range(n - 1):
            if x >= self.xData[i] and x <= self.xData[i + 1]:
                index = i
                break

        delta_x = x - self.xData[index]
        interpolated_value = (                      #полином степени
            self.yData[index]
            + self.b[index] * delta_x
            + self.c[index] * delta_x ** 2
            + self.d[index] * delta_x ** 3
        )

        return interpolated_value

# Функция f(x, y) для методов Рунге-Кутты
def f(x, y):
    return (math.exp(x) + y + (y - y * x)) / 3.0

# Метод Рунге-Кутта 2-го порядка
def runge_kutta2(x0, y0, x, h):
    n = int((x - x0) / h)
    y = y0

    print("Метод Рунге-Кутта 2-го порядка:")
    print(f"y({x0}) = {y}")

    for i in range(1, n + 1):
        k1 = h * f(x0, y)
        k2 = h * f(x0 + h, y + k1)
        y += 0.5 * (k1 + k2)
        x0 += h

        print("y({:.6f}) = {:.6f}".format(x0, y))

# Метод Рунге-Кутта 4-го порядка
def runge_kutta4(x0, y0, x, h):
    n = int((x - x0) / h)
    y = y0

    print("Метод Рунге-Кутта 4-го порядка:")
    print(f"y({x0}) = {y}")
    for i in range(1, n + 1):
        k1 = h * f(x0, y)
        k2 = h * f(x0 + 0.5 * h, y + 0.5 * k1)
        k3 = h * f(x0 + 0.5 * h, y + 0.5 * k2)
        k4 = h * f(x0 + h, y + k3)

        y += (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        x0 += h

        print("y({:.6f}) = {:.6f}".format(x0, y))

x0 = 0.0  # Начальное значение x
y0 = 1.0  # Начальное значение y
x = 1.0  # Конечное значение x
h = 0.1  # Шаг

runge_kutta2(x0, y0, x, h)
print()
runge_kutta4(x0, y0, x, h)
print()

xData = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]  # Значения x
yData = [y0]  # Начальное значение y

# Вычисляем значения y с использованием метода Рунге-Кутты 4-го порядка
xi = x0
yi = y0
n = int((x - x0) / h)

for i in range(1, n + 1):
    k1 = h * f(xi, yi)
    k2 = h * f(xi + 0.5 * h, yi + 0.5 * k1)
    k3 = h * f(xi + 0.5 * h, yi + 0.5 * k2)
    k4 = h * f(xi + h, yi + k3)

    yi += (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
    xi += h

    yData.append(yi)

print("Интерполяция кубическим сплайном:")
spline = CubicSpline(xData, yData)
for x in xData:
    interpolatedValue = spline.interpolate(x)
    print("P({:.6f}) = {:.6f}".format(x, interpolatedValue))
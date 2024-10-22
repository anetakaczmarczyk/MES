import math


class GaussTable:
    def __init__(self, n):
        self.n = n
        self.points = []
        self.weights = []
        if self.n == 1:
            self.points.append(0)
            self.weights.append(2)
        elif self.n == 2:
            self.points.append(-(1 / math.sqrt(3)))
            self.weights.append(1)
            self.points.append(1 / math.sqrt(3))
            self.weights.append(1)
        elif self.n == 3:
            self.points.append(-(math.sqrt(3 / 5)))
            self.weights.append(5 / 9)
            self.points.append(0)
            self.weights.append(8 / 9)
            self.points.append(math.sqrt(3 / 5))
            self.weights.append(5 / 9)
        elif self.n == 4:
            self.points.append(-0.861135)
            self.weights.append(0.347855)
            self.points.append(-0.339981)
            self.weights.append(0.652145)
            self.points.append(0.339981)
            self.weights.append(0.652145)
            self.points.append(0.861135)
            self.weights.append(0.347855)
        elif self.n == 5:
            self.points.append(-0.906180)
            self.weights.append(0.236927)
            self.points.append(-0.538469)
            self.weights.append(0.478629)
            self.points.append(0)
            self.weights.append(0.568889)
            self.points.append(0.538469)
            self.weights.append(0.478629)
            self.points.append(0.906180)
            self.weights.append(0.236927)
        else:
            print("ERROR")


def solve1d(fun, n):
    sum = 0
    gauss = GaussTable(n)
    for i in range(n):
        sum += gauss.weights[i] * fun(gauss.points[i])
    return sum


def solve2d(fun, n):
    sum = 0
    gauss = GaussTable(n)
    for i in range(n):
        for j in range(n):
            sum += gauss.weights[i] * gauss.weights[j] * fun(gauss.points[i], gauss.points[j])
    return sum


def fun1d(x):
    return 5 * x ** 2 + 3 * x + 6


def fun2d(x, y):
    return 5 * x ** 2 * y ** 2 + 3 * x * y + 6


print(solve1d(fun1d, 4))  # 15.333

print(solve2d(fun2d, 4))  # 26.222

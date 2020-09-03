import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

data = pd.read_csv('03_Владивосток.csv')
columns = data.columns[1:]
date_first, date_last = data['Год'].min(), data['Год'].max()

def approx_polynom(x, y, k):
    A = np.zeros(shape=(k + 1, k + 1))

    deg_sum = []
    for i in range(2 * k + 1):
        deg_sum.append(sum(j ** i for j in x))

    for i in range(k + 1):
        for j in range(k + 1):
            A[i][j] = deg_sum[i + j]

    b = np.array([sum(x[i] ** j * y[i] for i in range(len(x))) for j in range(k + 1)])

    return np.dot(np.linalg.inv(A), b)


def compute_polynom(coefs, p):
    res = 0
    for i, coef in enumerate(coefs):
        res += coef * p ** i
    return res


def interpolate_poly(month, deg):
    x = list(range(len(data[month])))
    y = list(data[month])
    if y[-1] == 999.9:
        y.pop()
        x.pop()
    coefs = approx_polynom(x, y, deg)

    plt.figure(figsize=(10, 7))
    plt.plot(np.array(x) + date_first, [compute_polynom(coefs, p) for p in x])
    plt.plot(np.array(x) + date_first, y, marker='o', linestyle='')
    plt.show()


if __name__ == '__main__':
    month, deg = input("введите месяц, степень многочлена\n").split()
    interpolate_poly(month, int(deg))
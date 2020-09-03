
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

data = pd.read_csv('03_Владивосток.csv')
columns = data.columns[1:]


def make_lagrange_polynom(x, y, p):
    res = 0
    for i in range(len(x)):
        frac = 1
        for j in range(len(x)):
            if j != i and j < len(x) and y[j] != 999.9:
                frac *= (p - x[j])/(x[i] - x[j])
        res += frac*y[i]
    return res

def compute_binom(q, n):
    res = 1
    for i in range(n):
        res *= (q - i) / (i + 1)
    return res


def fin_diff(y, n):
    fin_diffs = [y]
    for i in range(n):
        fin_diffs.append([])
        for j in range(n - i):
            fin_diffs[-1].append(fin_diffs[-2][j + 1] - fin_diffs[-2][j])
    return fin_diffs


def make_neuton_polynom_forward(x, y, p):
    h = x[1] - x[0]
    q = (p - x[0]) / h
    n = len(x) - 1
    fin_diffs = fin_diff(y, n)

    res = y[0]
    for i in range(1, n + 1):
        res += compute_binom(q, i) * fin_diffs[i][0]
    return res


def make_neuton_polynom_backward(x, y, p):
    h = x[1] - x[0]
    q = (p - x[-1]) / h
    n = len(x) - 1
    fin_diffs = fin_diff(y, n)

    res = y[-1]
    for i in range(1, n + 1):
        res += compute_binom(q + i - 1, i) * fin_diffs[i][-1]
    return res


def interpolate_weather(data, year, month, method, period):
    date_first, date_last = data['Год'].min(), data['Год'].max()
    start = year - date_first
    end = start + period
    grid = np.arange(year, year + period + 0.1, 0.1)
    x, y = data['Год'][start:end + 1].values, data[month][start:end + 1].values
    vals = list(map(lambda p: method(x, y, p), grid))

    plt.figure(figsize=(10, 7))
    plt.plot(grid, vals, label=month)
    plt.plot(x, y, marker='o', linestyle='dashed')
    plt.legend()
    plt.show()
'''
def interpolate_weather(data, year, month):
    date_first, date_last = data['Год'].min(), data['Год'].max()
    start = year - date_first
    end = start + 12
    grid = np.arange(year, year+12.1, 0.1)
    x, y = data['Год'][start:end+1], data[month][start:end+1]
    vals = list(map(lambda p: make_lagrange_polynom(x, y, p), grid))
    
    plt.figure(figsize=(10, 7))
    plt.plot(grid, vals, label=month)
    plt.plot(x.values, y.values, marker='o', linestyle='dashed')
    plt.legend()
'''


if __name__ == '__main__':
    year, month, period, method = input("введите год, месяц, период, метод\n").split()
    year = int(year)
    if method == 'Lagrange':
        interpolate_weather(data, year, month, make_lagrange_polynom, int(period))
    elif method == 'neuton_forward':
        interpolate_weather(data, year, month, make_neuton_polynom_forward, int(period))
    elif method == 'neuton_back':
        interpolate_weather(data, year, month, make_neuton_polynom_backward, int(period))





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
    plt.show()


if __name__ == '__main__':
    year, month = input("введите год и месяц\n").split()
    year = int(year)
    interpolate_weather(data, year, month)






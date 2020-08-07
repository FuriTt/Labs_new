import numpy as np
from numpy.linalg import inv
from copy import copy


class TeylorMethod():
    def compute(self, f, a, b, y0, partials, grid=None, n_init=10, eps=0.01):
        if grid is None:
            grid = self._get_optimal_grid(f, a, b, y0, partials, n_init=n_init, eps=eps)
        return self._teylor_method(f, grid, y0, partials)

    def _teylor_method(self, f, grid, y0, partials):
        result = []
        for i, x in enumerate(grid):
            if i == 0:
                result.append(y0)
            else:
                h = grid[i] - grid[i - 1]
                x, y = grid[i - 1], result[-1]
                new_y = result[-1] + h * (f(x, y) + 0.5 * h * (partials['x'](x, y) + partials['y'](x, y) * f(x, y)))
                result.append(new_y)
        return np.array(result)

    def _get_optimal_grid(self, f, a, b, y0, partials, n_init=10, eps=0.01):
        delta = np.inf
        h = (b - a) / n_init
        while delta > eps:
            y_old = self._teylor_method(f, np.arange(a, b, h), y0, partials)
            y_new = self._teylor_method(f, np.arange(a, b, h / 2), y0, partials)
            delta = max(np.abs(y_new[0::2] - y_old))
            h /= 2
        return np.arange(a, b, h)


class NeutonRafson():
        @staticmethod
        def neuton_rafson(f, nabla, gesse, p_init=None, eps=0.00001):
            dim = len(gesse)
            if p_init is None:
                p_init = NeutonRafson._swenn(f, dim)
            grad = np.array([part(*p_init) for part in nabla])
            p = np.array(p_init)
            gesse_val = np.matrix([[gesse[i][j](*p) for j in range(dim)] for i in range(dim)])

            while np.sqrt(sum(grad ** 2)) > eps:
                p -= np.array(inv(gesse_val).dot(grad))[0]
                grad = np.array([part(*p) for part in nabla])
                gesse_val = np.matrix([[gesse[i][j](*p) for j in range(dim)] for i in range(dim)])
            return p

        @staticmethod
        def _swenn(f, dim):
            p_init = np.random.random(size=dim)
            volume = []
            for coord in range(dim):
                volume.append(NeutonRafson._swenn_line(f, p_init, coord))
            return volume

        @staticmethod
        def _swenn_line(f, p_init, coord, h=0.001, max_iter=1000):
            p0, p1 = copy(p_init), copy(p_init)
            p1[coord] += h

            if f(*p0) < f(*p1):
                h *= -1
                p0, p1 = p1, p0

            p2 = copy(p1)
            p2[coord] += h

            iter_count = 0
            while f(*p1) > f(*p2):
                if iter_count > max_iter:
                    raise ValueError("Swenn method exceed max_iteration number")

                p_new = copy(p2)
                p_new[coord] += h
                p0, p1, p2 = p1, p2, p_new
                h *= 2
                iter_count += 1

            return np.array([p0[coord], p2[coord]])





if __name__ == '__main__':
    f = lambda x, y: -x * y + np.sin(x)
    a, b = 1, 10
    y0 = 0.5

    partials = {
        'x': lambda x, y: -y + np.cos(x),
        'y': lambda x, y: -x
    }

    teylor = TeylorMethod()
    print(teylor.compute(f, a, b, y0, partials))
    print('/////////////')
    f = lambda x, y: x ** 2 + y ** 2 + x ** 3

    nabla = np.array([lambda x, y: 2 * x + 3 * x ** 2,
                      lambda x, y: 2 * y])

    gesse = [[lambda x, y: 2 + 6 * x, lambda x, y: 0],
             [lambda x, y: 0, lambda x, y: 2]]

    p_init = (0.3, 0.9)
    neuton = NeutonRafson()
    print(NeutonRafson.neuton_rafson(f, nabla, gesse))

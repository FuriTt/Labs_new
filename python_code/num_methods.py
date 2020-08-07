import numpy as np
from numpy.linalg import inv


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
        def neuton_rafson(self,f, nabla, gesse, p_init, eps=0.00001):
            dim = len(gesse)
            grad = np.array([part(*p_init) for part in nabla])
            p = np.array(p_init)
            gesse_val = np.matrix([[gesse[i][j](*p) for j in range(dim)] for i in range(dim)])

            while np.sqrt(sum(grad ** 2)) > eps:
                p -= np.array(inv(gesse_val).dot(grad))[0]
                grad = np.array([part(*p) for part in nabla])
                gesse_val = np.matrix([[gesse[i][j](*p) for j in range(dim)] for i in range(dim)])
            return p


if __name__ == '__main__':
    f = lambda x, y: -x * y + np.sin(x)
    a, b = 1, 10
    y0 = 0.5

    partials = {
        'x': lambda x, y: -y + np.cos(x),
        'y': lambda x, y: -x
    }

    teylor = TeylorMethod()
    print(teylor.compute(f, a, b, y0, partials)
    print(/////////////)
    f = lambda x, y: x ** 2 + y ** 2 + x ** 3

    nabla = np.array([lambda x, y: 2 * x + 3 * x ** 2,
                      lambda x, y: 2 * y])

    gesse = [[lambda x, y: 2 + 6 * x, lambda x, y: 0],
             [lambda x, y: 0, lambda x, y: 2]]

    p_init = (0.3, 0.9)
    neuton = NeutonRafson()
    print(neuton.neuton_rafson(f, nabla, gesse, p_init))

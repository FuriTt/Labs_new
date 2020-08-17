#include <iostream>
#include "math.h"
#include <unordered_map>
#include <string>
#include <functional>
#include <vector>
#include<bits/stdc++.h>


class TeylorMethod {
public:
    TeylorMethod(int f_order): f_order(f_order) {};

    std::vector<std::vector<float>> compute(
            std::function<float(std::vector<float>&)>& f,
            float a,
            float b,
            std::vector<float>& y0,
            std::vector<std::function<float(std::vector<float>&)>>& partials,
            std::vector<float>& grid,
            int n_init=10,
            double eps=0.01
            ) {
        if (grid.empty()) {
            grid = _get_optimal_grid(f, a, b, y0, partials, n_init, eps);
        }
        return _teylor_method(f, grid, y0, partials);
    }


private:
    int f_order;

    std::vector<float> _get_optimal_grid(
             std::function<float(std::vector<float>&)>& f,
             float a,
             float b,
             std::vector<float>& y0,
             std::vector<std::function<float(std::vector<float>&)>>& partials,
             int n_init=10,
             double eps=0.01) const {
        auto delta = FLT_MAX;
        double h = (b - a) / float(n_init);
        auto grid = _make_grid(a, b, h);

        while (delta > eps) {
            grid = _make_grid(a, b, h);
            auto y_old = _teylor_method(f, grid, y0, partials);

            auto new_grid = _make_grid(a, b, h / 2);
            auto y_new = _teylor_method(f, new_grid, y0, partials);

            delta = 0;
            for (int j = 0; j < y_old.size(); ++j) {
                for (int i = 0; i < y_old[0].size() && 2 * i < y_new[0].size(); ++i) {
                    delta = std::max(delta, std::abs(y_new[j][2 * i] - y_old[j][i]));
                }
            }
            h /= 2;
        }
        return grid;
    }

    std::vector<std::vector<float>> _teylor_method(
            std::function<float(std::vector<float>&)>& f,
            std::vector<float>& grid,
            std::vector<float>& y0,
            std::vector<std::function<float(std::vector<float>&)>>& partials
            ) const {

        std::vector<std::vector<float>> result(f_order + 1, std::vector<float>(grid.size()));
        for (int j = 0; j < result.size(); ++j) {
            result[j][0] = y0[j];
        }

        for (int i = 1; i < grid.size(); ++i) {
            std::vector<float> prev_val;

            prev_val.push_back(grid[i - 1]);
            for (int j = 0; j < result.size(); ++j) {
                prev_val.push_back(result[i-1][j]);
            }

            result[result.size() - 1][i] = f(prev_val);
            if (result.size() < 2) {
                return result;
            }

            int der_ord = result.size() - 2;
            float step = grid[i] - grid[i - 1];
            result[der_ord][i] = result[der_ord][i-1] + step * f(prev_val);
            for (int j = 0; j < result.size(); ++j) {
                result[der_ord][i] += (step*step/2) * partials[j](prev_val) * (j == 0 ? 1 : result[j][i - 1]);
            }

            for (int j = 0; j + 2 < result.size(); ++j) {
                result[j][i] = result[j][i-1] + step * result[j + 1][i - 1] + (step * step / 2) * result[j + 2][i - 1];
            }
        }

        return result;
    }

    std::vector<float> _make_grid(float a, float b, float h) const {
        std::vector<float> grid;
        for (int i = 0; a + h*i < h; ++i) {
            grid.push_back(a + h*i);
        }
        return grid;
    }


};




float f(std::vector<float>& args) {
    return args[1] - args[0]*std::exp(args[0]);
}


int main() {
    TeylorMethod teylor(2);
    float a = 1, b = 10;
    float y0 = 0.5;

    std::vector<std::function<float(std::vector<float>&)>> partials = {
            [](std::vector<float>& args)->float { return -args[0]*std::exp(args[0]) - std::exp(args[0]);},
            [](std::vector<float>& args)->float { return -args[0]*std::exp(args[0]) - std::exp(args[0]);}
    };

    std::vector<std::vector<float>> result = teylor.compute(f, a, b, y0, partials);

    for (int i = 0; i < result.size(); ++i) {
        std::cout << i << " derivation: ";
        for (float val: result[i]) {
            std::cout << val << " ";
        }
        std::cout << "\n";
    }

}

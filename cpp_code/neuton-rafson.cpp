#include<vector>
#include<functional>
#include "Eigen/Dense"
#include <fstream>
#include <iostream>
#include <cmath>

using namespace Eigen;

class NeutonRafson {
public:
    std::vector<VectorXd> trajectory;

    NeutonRafson() {
        trajectory = std::vector<VectorXd>();
    }

    void compute(std::function<float(VectorXd&)>& f,
            std::vector<std::function<float(VectorXd&)>>& nabla,
            std::vector<std::vector<std::function<float(VectorXd&)>>>& gesse,
            VectorXd& p_init,
            float eps=0.0001,
            float learning_rate=0.1) {

        VectorXd p(p_init);
        float norm = std::pow(10, 9);

        while (norm > eps) {
            std::cout << "Norm: " << norm << "\n";
            trajectory.push_back(p);

            auto grad = _materialize_vector(nabla, p);
            auto gesse_val = _materialize_matrix(gesse, p);

            VectorXd direction = -gesse_val.inverse() * grad;
            std::cout << "Direction: " << direction << "\n";
            std::cout << "Gradient: " << grad << "\n";

            auto step = _get_step(f, p, direction, learning_rate);
            std::cout << "Step: " << step << "\n";
            p += step * learning_rate * direction;
            norm = grad.norm();
        }


    }

private:
    static float _get_step(
            std::function<float(VectorXd&)>& f,
            VectorXd& p,
            VectorXd& direction,
            float learning_rate) {

        std::function<float(float)> func = [&](float s)->float {
            VectorXd new_point = p + s * learning_rate * direction;
            return f(new_point);
        };

        /*
        float step = 0.01;
        float min_val = func(step);
        for (int i = 1; i < 30; ++i) {
            if (func(i*0.01) < min_val) {
                min_val = func(i*0.01);
                step = i*0.01;
            }
        }
        std::cout << "Step " << step << "\n";
        return step;
        */
        auto unimodal_interval = _get_unimodal_interval(func, learning_rate);
        return _get_min_uniform(func, unimodal_interval, learning_rate) / learning_rate;
    }

    static std::pair<float, float> _get_unimodal_interval(
            std::function<float(float)>& f,
            float lrate,
            float thr = 100) {
        float start = -lrate;
        float step = 1;
        if (f(start) < f(0)) {
            start *= -1;
            step *= -1;
        }
        float p_next = step * lrate;
        while ((f(p_next) < f(0)) & (step < thr) ) {
            step += 1;
            p_next = step * lrate;
        }

        std::cout << "Interval: " << start << " " << p_next << "\n";

        if (start < p_next) {
            return {lrate, p_next};
        }
        return {lrate, start};
    }

    static float _get_min_uniform(
            std::function<float(float)>& f,
            std::pair<float, float>& interval,
            float lrate) {
        float min_val = std::pow(10, 9);;
        float point;
        for (int i = 0; interval.first + i * lrate <= interval.second; ++i) {
            float new_point = interval.first + i * lrate;
            if (f(new_point) < min_val) {
                min_val = f(new_point);
                point = new_point;
            }
        }

        std::cout << "Point: " << point << "\n";

        return point;
    }


    static VectorXd _materialize_vector(std::vector<std::function<float(VectorXd&)>>& nabla, VectorXd& p) {
        VectorXd grad(nabla.size());
        for (int i = 0; i < nabla.size(); ++i) {
            grad(i) = nabla[i](p);
        }
        return grad;
    }

    static MatrixXd _materialize_matrix(
            std::vector<std::vector<std::function<float(VectorXd&)>>>& gesse,
            VectorXd& p) {
        auto gesse_val = MatrixXd(gesse.size(), gesse.size());
        for (int i = 0; i < gesse.size(); ++i) {
            for (int j = 0; j < gesse.size(); ++j) {
                gesse_val(i, j) = gesse[i][j](p);
            }
        }
        return gesse_val;
    }
};


float f(VectorXd& args) {
    return std::exp(3.5)*std::exp(args[0]*0.5)*(args[0] + std::pow(args[1], 2) - 8*args[1] + 23);
}


int main() {
    std::vector<std::function<float(VectorXd&)>> nabla = {
            [](VectorXd& args)->float {
                return 0.5*std::exp(3.5)*std::exp(args[0]*0.5)*(args[0] + std::pow(args[1], 2) - 8*args[1] + 25);
            },
            [](VectorXd& args)->float {
                return std::exp(3.5)*std::exp(args[0]*0.5)*(2*args[1] - 8);
            }
    };

    std::vector<std::vector<std::function<float(VectorXd&)>>> gesse = {
            {
                [](VectorXd& args)->float {
                    return 0.25*std::exp(3.5)*std::exp(args[0]*0.5)*(args[0] + std::pow(args[1], 2) - 8*args[1] + 27);
                },
                [](VectorXd& args)->float {
                    return 0.5*std::exp(3.5)*std::exp(args[0]*0.5)*(2*args[1] - 8);
                },
            },
            {
                [](VectorXd& args)->float {
                    return 0.5*std::exp(3.5)*std::exp(args[0]*0.5)*(2*args[1] - 8);
                },
                [](VectorXd& args)->float {
                    return 2*std::exp(3.5)*std::exp(args[0]*0.5);
                },

            }
    };

    VectorXd p_init(2);
    p_init << 8, 8;
    std::function<float(VectorXd&)> func(f);

    auto neuton_rafgon = NeutonRafson();
    neuton_rafgon.compute(func, nabla, gesse, p_init, 0.001, 0.01);

    std::ofstream output("neuton_rafson.csv");
    for (auto& p: neuton_rafgon.trajectory) {
        for (int i = 0; i < p.size(); ++i) {
            output << p[i] << ",";
        }
        output << f(p) << "\n";
    }
}
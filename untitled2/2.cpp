#include <vector>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;

vector<vector<float>> parse2DCsvFile(string inputFileName) {

    vector<vector<float> > data;
    ifstream inputFile(inputFileName);
    int l = 0;

    while (inputFile) {
        l++;
        string s;
        if (!getline(inputFile, s)) break;
        if (s[0] != '#') {
            istringstream ss(s);
            vector<float> record;

            while (ss) {
                string line;
                if (!getline(ss, line, ','))
                    break;
                try {
                    record.push_back(stof(line));
                }
                catch (const std::invalid_argument e) {
                    cout << "NaN found in file " << inputFileName << " line " << l
                         << endl;
                    e.what();
                }
            }

            data.push_back(record);
        }
    }

    if (!inputFile.eof()) {
        cerr << "Could not read file " << inputFileName << "\n";
        __throw_invalid_argument("File not found.");
    }

    return data;
}

/*
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
 */

VectorXd approx_polynom(std::vector<float>& x, std::vector<float>& y, int k) {
    MatrixXd A(k+1, k+1);
    vector<float> deg_sum(2*k+1);
    for (int i = 0; i < deg_sum.size(); ++i) {
        for (auto j: x) {
            deg_sum[i] += pow(j, i);
        }
    }

    for (int i = 0; i < k + 1; ++i) {
        for (int j = 0; j < k + 1; ++j) {
            A(i, j) = deg_sum[i + j];
        }
    }

    VectorXd b(k+1);
    for (int j = 0; j < k+1; ++j) {
        b(j) = 0;
        for (int i = 0; i < x.size(); ++i) {
            b(j) += pow(x[i], j)*y[i];
        }
    }

    return A.inverse()*b;
}

/*
 def compute_polynom(coefs, p):
    res = 0
    for i, coef in enumerate(coefs):
        res += coef * p ** i
    return res
 */


float compute_polynom(std::vector<float>& coefs, float p) {
    float res = 0;
    for (int i = 0; i < coefs.size(); ++i) {
        res += coefs[i] * pow(p, i);
    }
    return res;
}

std::vector<std::vector<float>> interpolate_poly(std::vector<std::vector<float>>& data, int month, int deg) {
    vector<float> x(data.size());
    for (int i = 0; i < x.size(); ++i) {
        x[i] = i;
    }

    vector<float> y(data.size());
    for (int i = 0; i < y.size(); ++i) {
        y[i] = data[i][month];
    }

    if (y.back() == 999.9) {
        y.pop_back();
        x.pop_back();
    }

    auto coefs_ = approx_polynom(x, y, deg);
    vector<float> coefs;
    for (int j = 0; j <= deg; ++j) {
        coefs.push_back(coefs_(j));
    }

    vector<float> predicted;
    for (auto el: x) {
        predicted.push_back(compute_polynom(coefs, el));
    }

    return {x, y, predicted};
}

int main() {
    vector<vector<float>> data = parse2DCsvFile("1.csv");

    ofstream learn_data("learn_data2.csv");
    ofstream predict_data("predict_data2.csv");

    auto result = interpolate_poly(data, 1, 3);

    auto x = result[0];
    auto y = result[2];
    for (int i = 0; i < x.size(); ++i) {
        predict_data << x[i] << "," << y[i] << "\n";
    }

    x = result[0];
    y = result[1];
    for (int i = 0; i < x.size(); ++i) {
        learn_data << x[i] << "," << y[i] << "\n";
    }

    predict_data.close();
    learn_data.close();


    return 0;
}


#include <vector>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;

float make_lagrange_polynom(std::vector<float>& x, std::vector<float>& y, float p) {
    float res = 0;
    for (int i = 0; i < x.size(); ++i) {
        float frac = 1;
        for (int j = 0; j < x.size(); ++j) {
            if (j != i && j < x.size() && y[j] != 999.9) {
                frac *= (p - x[j])/(x[i] - x[j]);
            }
        }
        res += frac*y[i];
    }
    return res;
}

float compute_binom(float q, int n) {
    float res = 1;
    for (int i = 0; i < n; ++i) {
        res *= float(q - i) / float(i + 1);
    }
    return res;
}


std::vector<std::vector<float>> fin_diff(std::vector<float>& y, int n) {
    std::vector<std::vector<float>> res(n+1, std::vector<float>());
    res[0] = y;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n - i; ++j) {
            res[i + 1].push_back(res[i][j+1] - res[i][j]);
        }
    }
    return res;
}

float make_neuton_polynom_forward(std::vector<float>& x, std::vector<float>& y, float p) {
    float h = x[1] - x[0];
    float q = (p - x[0]) / h;
    int n = x.size() - 1;

    auto fin_diffs = fin_diff(y, n);
    float res = y[0];

    for (int i = 1; i < n + 1; i++) {
        res += compute_binom(q, i) * fin_diffs[i][0];
    }
    return res;
}

float make_neuton_polynom_backward(std::vector<float>& x, std::vector<float>& y, float p) {
    float h = x[1] - x[0];
    float q = (p - x.back()) / h;
    int n = x.size() - 1;

    auto fin_diffs = fin_diff(y, n);
    float res = y.back();

    for (int i = 1; i < n + 1; i++) {
        res += compute_binom(q + i - 1, i) * fin_diffs[i].back();
    }
    return res;
}


std::vector<std::vector<float>> interpolate_weather(std::vector<std::vector<float>>& data,
        int year, int month, std::function<float(std::vector<float>&, std::vector<float>&, float)>& method,
        int period) {
    int date_first = 1872;
    int start = year - date_first;
    int end = start + period;
    std::vector<float> grid;
    for (float i = start; i <= start + period + 0.1; i += 0.1) {
        grid.push_back(i);
    }
    std::vector<float> x;
    std::vector<float> y;

    for (int i = start; i <= start + period; ++i) {
        x.push_back(i);
        y.push_back(data[i][month]);
    }

    std::vector<float> vals;
    for (auto p: grid) {
        vals.push_back(method(x, y, p));
    }

    return {x, y, grid, vals};
}



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


int main()
{
    vector<vector<float>> data = parse2DCsvFile("1.csv");

    ofstream learn_data("learn_data.csv");
    ofstream predict_data("predict_data.csv");

    auto lagrange = std::function<float(std::vector<float>&, std::vector<float>&, float)>(make_lagrange_polynom);
    auto neuton_fow = std::function<float(std::vector<float>&, std::vector<float>&, float)>(make_neuton_polynom_forward);
    auto neuton_back = std::function<float(std::vector<float>&, std::vector<float>&, float)>(make_neuton_polynom_backward);

    auto result = interpolate_weather(data, 1880, 1, neuton_back, 12);

    auto x = result[2];
    auto y = result[3];
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


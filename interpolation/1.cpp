#include <vector>
#include <fstream>

/*
def make_neuton_polynom_forward(x, y, p):
    h = x[1] - x[0]
    q = (p - x[0]) / h
    n = len(x) - 1
    fin_diffs = fin_diff(y, n)

    res = y[0]
    for i in range(1, n + 1):
        res += compute_binom(q, i) * fin_diffs[i][0]
    return res
*/

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
    return res''
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


class CSVRow
{
public:
    std::string_view operator[](std::size_t index) const
    {
        return std::string_view(&m_line[m_data[index] + 1], m_data[index + 1] -  (m_data[index] + 1));
    }
    std::size_t size() const
    {
        return m_data.size() - 1;
    }
    void readNextRow(std::istream& str)
    {
        std::getline(str, m_line);

        m_data.clear();
        m_data.emplace_back(-1);
        std::string::size_type pos = 0;
        while((pos = m_line.find(',', pos)) != std::string::npos)
        {
            m_data.emplace_back(pos);
            ++pos;
        }
        // This checks for a trailing comma with no data after it.
        pos   = m_line.size();
        m_data.emplace_back(pos);
    }
private:
    std::string         m_line;
    std::vector<int>    m_data;
};

std::istream& operator>>(std::istream& str, CSVRow& data)
{
    data.readNextRow(str);
    return str;
}






int main()
{
    std::ifstream  file("03_Владивосток.csv");
    CSVRow  row;
    while(file >> row)
    {
        std::cout << "4th Element(" << row[3] << ")\n";
    }
}


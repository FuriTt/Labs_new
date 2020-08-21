#include<vector>
#include<functional>

class NeutonRafson {
public:
    NeutonRafson() = default;

    std::vector<float> compute(std::functional<float(std::vector<float>&)>& f,
            std::vector<std::functional<float(std::vector<float>&)>>& nabla,
            std::vector<std::vector<std::functional<float(std::vector<float>&)>>>& gesse,
            vector<float>& p_init,
            float eps=0.00001) {

        auto grad = std::vector<float>;
        for (auto& part: nabla) {
            grad.push_back(part(p_init));
        }

        std::vector<std::vector<float>> gesse_val(gesse.size(), std::vector<float>(gesse.size()));
        for (int i = 0; i < gesse.size(); ++i) {
            for (int j = 0; j < gesse.size(); ++j) {
                gesse_val[i][j] = gesse[i][j](p_init);
            }
        }

        while (_get_norm(grad) > eps) {

        }

        auto inv_gesse = _inv(gesse_val);
        p -= _dot(inv_gesse, grad);

        /*
        grad = np.array([part(*p) for part in nabla])
        gesse_val = np.matrix([[gesse[i][j](*p) for j in range(dim)] for i in range(dim)])
         */
        //TODO

    }

private:
    float _get_norm(std::vector<float>& v) {
        //TODO
    }

    std::vector<std::vector<float>> _inv(std::vector<std::vector<float>>& matrix) {
        //TODO
    }

    std::vector<float> _dot(std::vector<std::vector<float>>& matrix,
            std::vector<float>& v) {
        //TODO
    }
};
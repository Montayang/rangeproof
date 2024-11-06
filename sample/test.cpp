#include <mcl/bls12_381.hpp>
#include "subgroups.h"
#include "kzg.h"
#include <vector>

constexpr int MAX_N = 1 << 5;

using namespace mcl::bls12;

Fr tmp1[MAX_N], tmp2[MAX_N]; 

void INTT(Fr* f, __int128 n, const subgroup& group) {
    if (n == 1) return;
    for (__int128 i = 0; i < n; ++i) tmp1[i] = f[i];
    for (__int128 i = 0; i < n; ++i) {
        if (i & 1)
            f[n / 2 + i / 2] = tmp1[i];
        else
            f[i / 2] = tmp1[i];
    }
    Fr *g = f, *h = f + n / 2;
    INTT(g, n / 2, group), INTT(h, n / 2, group);
    Fr cur = 1, step = -group.getGenerator();
    for (__int128 k = 0; k < n / 2; ++k) {
        tmp1[k] = g[k] + cur * h[k];
        tmp1[k + n / 2] = g[k] - cur * h[k];
        cur *= step;
    }
    for (__int128 i = 0; i < n; ++i) f[i] = tmp1[i];
}

void INTT_d2(Fr (&f)[MAX_N][MAX_N], __int128 n1, __int128 n2, const subgroup& group1, const subgroup& group2) {
    for (__int128 i = 0; i < n1; ++i) INTT(f[i], n2, group2);
    for (__int128 j = 0; j < n2; ++j) {
        for (__int128 i = 0; i < n1; ++i) tmp2[i] = f[i][j];
        INTT(tmp2, n1, group1);
        for (__int128 i = 0; i < n1; ++i) f[i][j] = tmp2[i];
    }
}

int main() {
    initPairing(mcl::BLS12_381);

    //INTT test
    __int128 l = 2, k = 4;
    subgroup H_l(l), H_k(k);

    Fr y[MAX_N][MAX_N];
    y[0][0] = 2;
    y[0][1] = -1;
    y[0][2] = -2;
    y[0][3] = 7;
    y[1][0] = 4;
    y[1][1] = -2;
    y[1][2] = -3;
    y[1][3] = 1;

    INTT_d2(y, l, k, H_l, H_k);
    std::vector<std::vector<Fr>> coeffs;
    for (__int128 i = 0;i< l;i++) {
        std::vector<Fr> vec;
        coeffs.push_back(vec);
        for (__int128 j = 0;j< k;j++) {
            coeffs[i].push_back(y[i][j]/l/k);
            //std::cout << y[i][j]<<std::endl;
        }
    }
    Fr gx = 1;
    for (__int128 i = 0;i< l;i++) {
        Fr gy = 1;
        for (__int128 j = 0;j< k;j++) {
            //std::cout << evaluate_polynomial_d2(coeffs, gx, gy)<<std::endl;
            gy *= H_k.getGenerator();
        }
        gx *= H_l.getGenerator();
    }

    //Polynomial division test
    // std::vector<Fr> f, d;
    // f.push_back(3);
    // f.push_back(4);
    // f.push_back(3);
    // f.push_back(1);
    // d.push_back(1);
    // d.push_back(1);
    // d.push_back(1);
    // std::pair<std::vector<Fr>, std::vector<Fr>> res = polynomial_division(f, d);
    // std::cout << "quotient: " << std::endl;
    // for (auto i : res.first) std::cout << i << std::endl;
    // std::cout << "remainder: " << std::endl;
    // for (auto i : res.second) std::cout << i << std::endl;

    //Polynomial division d2 test
    std::vector<Fr> x0={2,2,1}, x1={2,2,1}, x2={1,1,0};
    std::vector<std::vector<Fr>> ff;
    ff.push_back(x0);
    ff.push_back(x1);
    ff.push_back(x2);
    std::vector<Fr> d1={1,1}, d2={1,1};
    std::pair<std::vector<std::vector<Fr>>, std::vector<std::vector<Fr>>> res = polynomial_division_d2(ff,d1,d2);

    std::cout << "1:" <<std::endl;
    for (auto i : res.first) {
        for (auto j : i) std::cout << j <<"  ";
        std::cout << std::endl;
    }
    std::cout << "2:" <<std::endl;
    for (auto i : res.second) {
        for (auto j : i) std::cout << j <<"  ";
        std::cout << std::endl;
    }

    return 0;
}
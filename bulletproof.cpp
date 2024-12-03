#include <iostream>
#include <vector>
#include <mcl/bn.hpp>
#include <random>
#include <cstring>
#include <cmath>
#include <algorithm>
#include "kzg.h"
#include "subgroups.h"
#include <chrono>

using namespace mcl::bn;

constexpr int n = 4;
constexpr int R = 1 << n;
constexpr int MAX_N = 1 << 25;
Fr* tmp1 = new Fr[MAX_N]; //for NTT innerly
Fr* y_1 = new Fr[MAX_N];
Fr* y_2 = new Fr[MAX_N]; // for poly_product and polynomial_extension

void print_int128(__int128 x) {
    if (x == 0) std::cout << "0";

    bool isNegative = (x < 0);
    if (isNegative) x = -x;  

    std::string result;
    while (x > 0) {
        result = char('0' + (x % 10)) + result;
        x /= 10;
    }

    if (isNegative) result = '-' + result; 

    std::cout << result << std::endl;
}

std::string generate_eff(__int128 x) {
    // 使用std::random_device和std::mt19937_64生成随机数
    std::random_device rd;  // 用于生成随机种子
    std::mt19937_64 gen(rd());  // 使用64位的Mersenne Twister引擎
    std::uniform_int_distribution<uint64_t> dist(0, x - 1);  // 分布范围为[0, x - 1]

    uint64_t randomValue = dist(gen);  // 生成随机值
    return std::to_string(randomValue);  // 转为字符串并返回
}

__int128 Fr_to_int(Fr x) {
    uint8_t bf[32]={0};
    int size = x.getLittleEndian(bf,32);	
    __int128 V = 0;	
    for(int j = size - 1; j >= 0; j--) {	
        V = (V << 8) + bf[j];	
    }
    return V;	
}

std::vector<Fr> generate_polynomial(const int &l) {
    std::vector<Fr> coeffs;
    Fr r;
    for (int i = 0; i < l; i++) {
        r.setStr(generate_eff(R));
        coeffs.push_back(r);
    }
    return coeffs;
}

void NTT(Fr* f, int n, const subgroup& group, int reverse) {
    if (n == 1) return;
    for (int i = 0; i < n; ++i) tmp1[i] = f[i];
    for (int i = 0; i < n; ++i) {
        if (i & 1)
            f[n / 2 + i / 2] = tmp1[i];
        else
            f[i / 2] = tmp1[i];
    }
    Fr *g = f, *h = f + n / 2;
    subgroup newgroup(n / 2);
    NTT(g, n / 2, newgroup, reverse), NTT(h, n / 2, newgroup, reverse);
    Fr cur = 1;
    Fr step = group.getGenerator();
    if (reverse == -1) step = 1 / step;
    for (int k = 0; k < n / 2; ++k) {
        tmp1[k] = g[k] + cur * h[k];
        tmp1[k + n / 2] = g[k] - cur * h[k];
        cur *= step;
    }
    for (int i = 0; i < n; ++i) f[i] = tmp1[i];
}

std::vector<Fr> poly_product(const std::vector<Fr>& poly1, const std::vector<Fr>& poly2) {
    int n = 1;
    while (n < poly1.size() + poly2.size()) {
        n <<= 1;
    }
    subgroup group(n);

    for (int i = 0; i < n; i++) {
        if (i < poly1.size()) y_1[i] = poly1[i];
        else y_1[i] = 0;
        if (i < poly2.size()) y_2[i] = poly2[i];
        else y_2[i] = 0;
    }
    NTT(y_1, n, group, 1);
    NTT(y_2, n, group, 1);

    for (int i = 0; i < n; i++) {
        y_1[i] = y_1[i] * y_2[i];
    }
    NTT(y_1, n, group, -1);
    // NTT(y_2, n, group, -1);
    // for (int i = 0; i < n; i++) {std::cout<< y_1[i] << std::endl;}
    // for (int i = 0; i < n; i++) {std::cout<< y_2[i] << std::endl;}
    std::vector<Fr> result;
    for (int i = 0; i < n; i++) result.push_back(y_1[i] / n);
    return result;
}

Fr vec_product(std::vector<Fr> vec1, std::vector<Fr> vec2) {
    Fr res = 0;
    int size = vec1.size();
    for (int i = 0; i < size; i++) res += vec1[i] * vec2[i];
    return res;
}

int main() {
    initPairing(mcl::BN_SNARK1);
    G1 g1;
    hashAndMapToG1(g1, "ab", 2);
    G2 g2;
    hashAndMapToG2(g2, "abc", 3);
    std::chrono::duration<double> prover_time(0);
    std::chrono::duration<double> verifier_time(0);
    size_t proof_size = 0;

    //input
    Fr v;
    v.setStr(generate_eff(R));
        //v.setStr("1048577"); //for testing over range
    subgroup H(n);
    Fr g = H.getGenerator();
    Fr h, h_seed;
    h_seed.setStr(generate_eff(n));
    Fr::pow(h, g, h_seed);
    Fr V, gamma;
    gamma.setRand();
    Fr tmp_V;
    Fr::pow(tmp_V, g, v);
    Fr::pow(V, h, gamma);
    V *= tmp_V;
    
    //bit decomposition
    std::vector<__int128> G;
    __int128 base = 1;
    for (int i = 0; i < n; i++) {
        G.push_back(Fr_to_int(base));
        base *= 2;
    }
    std::vector<Fr> a_L,a_R;
    __int128 v_int = Fr_to_int(v);
    for (int i = n - 1; i >= 0; i--) {
        a_L.push_back(v_int / G[i]);
        v_int = v_int % G[i];
    }
    for (int i = 0; i < n; i++) a_R[i] = a_L[i] - 1;

    Fr alpha;
    alpha.setRand();
    Fr A = 0;
    for (int i = 0; i < n; i++) A += a_R[i];
    A += alpha;
    Fr::pow(A, h, A);
    Fr tmp_A = 0;
    for (int i = 0; i < n; i++) tmp_A += a_L[i];
    Fr::pow(tmp_A, g, tmp_A);
    A *= tmp_A;
    
    std::vector<Fr> s_L, s_R;
    for (int i = 0; i < n; i++) {
        Fr s_L_seed, s_R_seed;
        s_L_seed.setRand();
        s_R_seed.setRand();
        s_L.push_back(s_L_seed);
        s_R.push_back(s_R_seed);
    }

    Fr rou;
    rou.setRand();
    Fr S = 0;
    for (int i = 0; i < n; i++) S += s_R[i];
    S += rou;
    Fr::pow(S, h, S);
    Fr tmp_S = 0;
    for (int i = 0; i < n; i++) tmp_S += s_L[i];
    Fr::pow(tmp_S, g, tmp_S);
    S *= tmp_S;

    //A, S sent
    Fr y, z;
    y.setRand();
    z.setRand();

    std::vector<std::vector<Fr>> l_poly, r_poly; //vector polynomial
    std::vector<Fr> t_poly;
    std::vector<Fr> tmp_l_poly(n, 0);
    for (int i = 0; i < n; i++) tmp_l_poly[i] = a_L[i] - z;
    l_poly.push_back(tmp_l_poly);
    l_poly.push_back(s_L);
    std::vector<Fr> tmp_r_poly(n, 0);
    Fr cur = 1;
    for (int i = 0; i < n; i++) {
        tmp_r_poly[i] = a_R[i] + z;
        tmp_r_poly[i] *= cur;
        cur *= y;
    }
    cur = 1;
    for (int i = 0; i < n; i++) {
        tmp_r_poly[i] += z * z * cur;
        cur *= 2;
    }
    r_poly.push_back(tmp_r_poly);
    cur = 1;
    for (int i = 0; i < n; i++) {
        tmp_r_poly[i] = s_R[i] * cur;
        cur *= y;
    }
    r_poly.push_back(tmp_r_poly);
    Fr t_0, t_1, t_2;
    t_0 = vec_product(l_poly[0], r_poly[0]);
    t_1 = vec_product(l_poly[0], r_poly[1]) + vec_product(l_poly[1], r_poly[0]);
    t_2 = vec_product(l_poly[1], r_poly[1]);
    t_poly.push_back(t_0);
    t_poly.push_back(t_1);
    t_poly.push_back(t_2);

    Fr tao_1, tao_2;
    tao_1.setRand();
    tao_2.setRand();
    Fr T_1, T_2;
    Fr tmp_T;
    Fr::pow(tmp_T, h, tao_1);
    Fr::pow(T_1, g, t_1);
    T_1 *= tmp_T;
    Fr::pow(tmp_T, h, tao_2);
    Fr::pow(T_2, g, t_2);
    T_2 *= tmp_T;

    //T_1, T_2 sent
    Fr x;
    x.setRand();

    std::vector<Fr> l = l_poly[0];
    for (int i = 0; i < n; i++) l[i] += l_poly[1][i] * x;
    std::vector<Fr> r = r_poly[0];
    for (int i = 0; i < n; i++) r[i] += r_poly[1][i] * x;
    Fr t_hat = vec_product(l, r);
    Fr tao_x = tao_2 * x * x + tao_1 * x + z * z * gamma;
    Fr miu = alpha + rou * x;

    

    std::cout << "1" << std::endl;

    std::cout << "Prover Time:" << prover_time.count() << std::endl;
    std::cout << "Verifier Time:" << verifier_time.count() << std::endl;
    std::cout << "Proof Size:" << ((double) proof_size) / 1024.0 << " KB" << std::endl;
    return 0;
}
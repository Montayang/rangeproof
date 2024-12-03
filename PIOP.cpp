<<<<<<< HEAD
#include <iostream>
#include <vector>
#include <mcl/bls12_381.hpp>
#include <random>
#include <cstring>
#include <cmath>
#include "kzg.h"
#include "subgroups.h"

using namespace mcl::bls12;

constexpr int MAX_N = 1 << 6;
Fr tmp1[MAX_N], tmp2[MAX_N]; //for INTT

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

__int128 generate_eff() {
    // 随机数引擎
    std::random_device rd;   // 使用真正的随机数种子
    std::mt19937_64 gen(rd()); // 使用64位的Mersenne Twister引擎
    std::uniform_int_distribution<uint64_t> dis(1, std::numeric_limits<uint32_t>::max());

    // 生成 1 到 2^32 - 1 的随机数并转换为 __int128 类型
    return static_cast<__int128>(dis(gen));
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

std::vector<Fr> generate_polynomial(const __int128 &l) {
    std::vector<Fr> coeffs;
    Fr r;
    for (int i = 0; i < l; i++) {
        r = generate_eff();
        //r = (Fr_to_int(r));
        coeffs.push_back(r);
    }
    return coeffs;
}

void NTT(Fr* f, __int128 n, const subgroup& group, __int128 reverse) {
    if (n == 1) return;
    for (__int128 i = 0; i < n; ++i) tmp1[i] = f[i];
    for (__int128 i = 0; i < n; ++i) {
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
    for (__int128 k = 0; k < n / 2; ++k) {
        tmp1[k] = g[k] + cur * h[k];
        tmp1[k + n / 2] = g[k] - cur * h[k];
        cur *= step;
    }
    for (__int128 i = 0; i < n; ++i) f[i] = tmp1[i];
}

void NTT_d2(Fr (&f)[MAX_N][MAX_N], __int128 n1, __int128 n2, const subgroup& group1, const subgroup& group2, __int128 reverse) {
    for (__int128 i = 0; i < n1; ++i) NTT(f[i], n2, group2, reverse);
    for (__int128 j = 0; j < n2; ++j) {
        for (__int128 i = 0; i < n1; ++i) tmp2[i] = f[i][j];
        NTT(tmp2, n1, group1, reverse);
        for (__int128 i = 0; i < n1; ++i) f[i][j] = tmp2[i];
    }
}

std::vector<Fr> poly_product(const std::vector<Fr>& poly1, const std::vector<Fr>& poly2) {
    __int128 n = 1;
    while (n < poly1.size() + poly2.size()) {
        n <<= 1;
    }
    subgroup group(n);

    Fr y_1[MAX_N], y_2[MAX_N];
    for (__int128 i = 0; i < n; i++) {
        if (i < poly1.size()) y_1[i] = poly1[i];
        else y_1[i] = 0;
        if (i < poly2.size()) y_2[i] = poly2[i];
        else y_2[i] = 0;
    }
    NTT(y_1, n, group, 1);
    NTT(y_2, n, group, 1);
    for (__int128 i = 0; i < n; i++) {
        y_1[i] = y_1[i] * y_2[i];
    }
    NTT(y_1, n, group, -1);
    // NTT(y_2, n, group, -1);
    // for (__int128 i = 0; i < n; i++) {std::cout<< y_1[i] << std::endl;}
    // for (__int128 i = 0; i < n; i++) {std::cout<< y_2[i] << std::endl;}
    std::vector<Fr> result;
    for (__int128 i = 0; i < n; i++) result.push_back(y_1[i]/n);
    return result;
}

std::vector<std::vector<Fr>> poly_product_d2(const std::vector<std::vector<Fr>>& poly1, const std::vector<std::vector<Fr>>& poly2) {
    __int128 n1 = 1, n2 = 1;
    while (n1 < poly1.size() + poly2.size()) {
        n1 <<= 1;
    }
    while (n2 < poly1[0].size() + poly2[0].size()) {
        n2 <<= 1;
    }
    subgroup group1(n1), group2(n2);

    Fr yy_1[MAX_N][MAX_N], yy_2[MAX_N][MAX_N];
    for (__int128 i = 0; i < n1; i++) {
        for (__int128 j = 0; j < n2; j++) {
            if (i < poly1.size() && j < poly1[0].size()) yy_1[i][j] = poly1[i][j];
            else yy_1[i][j] = 0;
            if (i < poly2.size() && j < poly2[0].size()) yy_2[i][j] = poly2[i][j];
            else yy_2[i][j] = 0;
        }
    }
    NTT_d2(yy_1, n1, n2, group1, group2, 1);
    NTT_d2(yy_2, n1, n2, group1, group2, 1);
    for (__int128 i = 0; i < n1; i++) {
        for (__int128 j = 0; j < n2; j++) {
            yy_1[i][j] *= yy_2[i][j];
        }
    }
    NTT_d2(yy_1, n1, n2, group1, group2, -1);
    std::vector<std::vector<Fr>> result;
    for (__int128 i = 0; i < n1; i++) {
        std::vector<Fr> tmp;
        result.push_back(tmp);
        for (__int128 j = 0; j < n2; j++) {
            result[i].push_back(yy_1[i][j] / n1 / n2);
        }
    }
    return result;
}

std::vector<Fr> sumset_representation(const Fr& R, const Fr& b) {
    __int128 R_j = Fr_to_int(R), b_int = Fr_to_int(b);
    std::vector<Fr> G;
    __int128 G_j = (R_j + 1) / b_int;
    while (R_j > b_int) {
        G.push_back(G_j);
        R_j = R_j - (b_int - 1) * G_j;
        G_j = (R_j + 1) / b_int;
    }
    if (G[G.size() - 1] != 1) G.push_back(1);
    return G;
}

int main() {
    initPairing(mcl::BLS12_381);
    G1 g1;
    hashAndMapToG1(g1, "ab", 2);
    G2 g2;
    hashAndMapToG2(g2, "abc", 3);

    //input
    Fr R;
    __int128 b, k, l;
    R.setStr("18446744073709551615"); //2^64 - 1
    b = 16, k = 16, l = 32;
    std::vector<Fr> a = generate_polynomial(l);
    CRS crs = kzg_setup(g1, g2, 2 * l + 2, 2 * k + 2);

    //a[6].setStr("18446744073709551616"); //for testing over range

    std::vector<Fr> G = sumset_representation(R, b);
    subgroup H_l(l), H_k(k), H_b(b);

    // bit decomposition
    std::vector<__int128> GG;
    for (auto i : G) {
        GG.push_back(Fr_to_int(i));
    }
    std::vector<std::vector<Fr>> d;
    for (__int128 i = 0; i < l; i++) {
        std::vector<Fr> tmp;
        d.push_back(tmp);
        __int128 ai = Fr_to_int(a[i]);
        for (__int128 j = 0; j < k; j++) {
            d[i].push_back(ai / GG[j]);
            ai = ai % GG[j];
        }
    }
    std::vector<Fr> m;
    for (__int128 i = 0; i < b; i++) m.push_back(0);
    for (auto i : d) {
        for (auto j : i) {
            __int128 index = Fr_to_int(j);
            if (index >= 0 && index < b) m[index] += 1;
        }
    }
    
    //Step 1:
    Fr y[MAX_N], yy[MAX_N][MAX_N];
    for (__int128 i = 0; i < l; i++) {// d_extension
        for (__int128 j = 0; j < k; j++) yy[i][j] = d[i][j];
    }
    NTT_d2(yy, l, k, H_l, H_k, -1);
    std::vector<std::vector<Fr>> d_poly;
    for (__int128 i = 0; i < l; i++) {
        std::vector<Fr> vec;
        d_poly.push_back(vec);
        for (__int128 j = 0; j < k; j++) {
            d_poly[i].push_back(yy[i][j]/l/k);
        }
    }

    for (__int128 i = 0; i < l; i++) {// a_extension
        y[i] = a[i];
    }
    NTT(y, l, H_l, -1);
    std::vector<Fr> a_poly;
    for (__int128 i = 0; i < l; i++) {
        a_poly.push_back(y[i]/l);
    }

    for (__int128 i = 0; i < k; i++) {// c_extension
        y[i] = G[i];
    }
    NTT(y, k, H_k, -1);
    std::vector<Fr> c_poly;
    for (__int128 i = 0; i < k; i++) {
        c_poly.push_back(y[i]/k);
    }

    for (__int128 i = 0; i < b; i++) {// I_extension
        y[i] = i;
    }
    NTT(y, b, H_b, -1);
    std::vector<Fr> I_poly;
    for (__int128 i = 0; i < b; i++) {
        I_poly.push_back(y[i] / b);
    }

    for (__int128 i = 0; i < b; i++) {// m_extension
        y[i] = m[i];
    }
    NTT(y, b, H_b, -1);
    std::vector<Fr> m_poly;
    for (__int128 i = 0; i < b; i++) {
        m_poly.push_back(y[i] / b);
    }
    G1 d_commitment = kzg_commit_d2(d_poly, crs);
    //assert(verify_kzg_open_d2(d_commitment, d_poly, crs));
    G1 a_commitment = kzg_commit(a_poly, crs);
    //assert(verify_kzg_open(a_commitment, a_poly, crs));
    G1 c_commitment = kzg_commit(c_poly, crs);
    //assert(verify_kzg_open(c_commitment, c_poly, crs));
    G1 I_commitment = kzg_commit(I_poly, crs);
    //assert(verify_kzg_open(I_commitment, I_poly, crs));
    G1 m_commitment = kzg_commit(m_poly, crs);
    //assert(verify_kzg_open(m_commitment, m_poly, crs));

    //Step2:
    Fr alpha, r_0;
    alpha.setRand();
    r_0.setRand();
    std::vector<std::vector<Fr>> T;
    for (__int128 i = 0; i < l; i++) {
        std::vector<Fr> tmp;
        T.push_back(tmp);
        for (__int128 j = 0; j < k; j++) {
            T[i].push_back(1 / (alpha + d[i][j]));
        }
    }
    std::vector<Fr> e;
    for (__int128 i = 0; i < b; i++) e.push_back(m[i] / (alpha + i));
    for (__int128 i = 0; i < l; i++) {// T_extension
        for (__int128 j = 0; j < k; j++) yy[i][j] = T[i][j];
    }
    NTT_d2(yy, l, k, H_l, H_k, -1);
    std::vector<std::vector<Fr>> T_poly;
    for (__int128 i = 0; i < l; i++) {
        std::vector<Fr> vec;
        T_poly.push_back(vec);
        for (__int128 j = 0; j < k; j++) {
            T_poly[i].push_back(yy[i][j]/l/k);
        }
    }

    for (__int128 i = 0; i < b; i++) {// e_extension
        y[i] = e[i];
    }
    NTT(y, b, H_b, -1);
    std::vector<Fr> e_poly;
    for (__int128 i = 0; i < b; i++) {
        e_poly.push_back(y[i]/b);
    }
    G1 T_commitment = kzg_commit_d2(T_poly, crs);
    //assert(verify_kzg_open_d2(T_commitment, T_poly, crs));
    G1 e_commitment = kzg_commit(e_poly, crs);
    //assert(verify_kzg_open(e_commitment, e_poly, crs));

    //Step3:
    std::vector<Fr> F1;
    std::vector<Fr> d_r0_poly;
    for (__int128 j = 0; j < d_poly[0].size(); j++) {
        std::vector<Fr> tmp;
        for (__int128 i = 0; i < d_poly.size(); i++) tmp.push_back(d_poly[i][j]);
        d_r0_poly.push_back(evaluate_polynomial(tmp, r_0));
    }
    F1 = poly_product(c_poly, d_r0_poly);
    F1[0] -= evaluate_polynomial(a_poly, r_0) / k;
    std::vector<Fr> Z_H_k(k + 1, 0);
    Z_H_k[0] = -1, Z_H_k[k] = 1;
    std::pair<std::vector<Fr>, std::vector<Fr>> tmp_pair = polynomial_division(F1, Z_H_k);
    std::vector<Fr> h_1_poly = tmp_pair.first;
    std::vector<Fr> g_1_poly = tmp_pair.second;
    assert(g_1_poly[0] == 0); //check remainder have no constant when divide F1
    for (__int128 i = 0; i < g_1_poly.size() - 1; i++) g_1_poly[i] = g_1_poly[i + 1];
    g_1_poly[g_1_poly.size() - 1] = 0;

    std::vector<std::vector<Fr>> F2 = d_poly;
    F2[0][0] += alpha;
    F2 = poly_product_d2(F2, T_poly);
    F2[0][0] -= 1;
    std::vector<Fr> Z_H_l(l + 1, 0);
    Z_H_l[0] = -1, Z_H_l[l] = 1;
    std::pair<std::vector<std::vector<Fr>>, std::vector<std::vector<Fr>>> tmp_pair2 = polynomial_division_d2(F2, Z_H_l, Z_H_k);
    std::vector<std::vector<Fr>> p_poly = tmp_pair2.first;
    std::vector<std::vector<Fr>> q_poly = tmp_pair2.second;

    std::vector<Fr> F3 = I_poly;
    F3[0] += alpha;
    F3 = poly_product(F3, e_poly);
    for (__int128 i = 0; i < b; i++) {
        F3[i] -= m_poly[i];
    }
    std::vector<Fr> Z_H_b(b + 1, 0);
    Z_H_b[0] = -1, Z_H_b[b] = 1;
    std::pair<std::vector<Fr>, std::vector<Fr>> tmp_pair3 = polynomial_division(F3, Z_H_b);
    for (auto i : tmp_pair3.second) assert(i == 0); //check no remainder when divide F3
    std::vector<Fr> h_2_poly = tmp_pair3.first;
    G1 h_1_commitment = kzg_commit(h_1_poly, crs);
    //assert(verify_kzg_open(h_1_commitment, h_1_poly, crs));
    G1 g_1_commitment = kzg_commit(g_1_poly, crs);
    //assert(verify_kzg_open(g_1_commitment, g_1_poly, crs));
    G1 p_commitment = kzg_commit_d2(p_poly, crs);
    //assert(verify_kzg_open_d2(p_commitment, p_poly, crs));
    G1 q_commitment = kzg_commit_d2(q_poly, crs);
    //assert(verify_kzg_open_d2(q_commitment, q_poly, crs));
    G1 h_2_commitment = kzg_commit(h_2_poly, crs);
    //assert(verify_kzg_open(h_2_commitment, h_2_poly, crs));

    //Step3 check:
    std::pair<Fr, G1> e_oracle = kzg_createWitness(e_poly, crs, r_0);
    assert(kzg_verifyEval(e_commitment, crs, r_0, e_oracle.first, e_oracle.second));
    std::pair<Fr, G1> I_oracle = kzg_createWitness(I_poly, crs, r_0);
    assert(kzg_verifyEval(I_commitment, crs, r_0, I_oracle.first, I_oracle.second));
    std::pair<Fr, G1> m_oracle = kzg_createWitness(m_poly, crs, r_0);
    assert(kzg_verifyEval(m_commitment, crs, r_0, m_oracle.first, m_oracle.second));
    std::pair<Fr, G1> h_2_oracle = kzg_createWitness(h_2_poly, crs, r_0);
    assert(kzg_verifyEval(h_2_commitment, crs, r_0, h_2_oracle.first, h_2_oracle.second));
    Fr Z_H_b_eval = evaluate_polynomial(Z_H_b, r_0);
    assert(e_oracle.first * (alpha + I_oracle.first) - m_oracle.first == h_2_oracle.first * Z_H_b_eval);

    std::pair<Fr, std::pair<G1, G1>> T_oracle = kzg_createWitness_d2(T_poly, crs, 0, 0);
    assert(kzg_verifyEval_d2(T_commitment, crs, 0, 0, T_oracle.first, T_oracle.second));
    e_oracle = kzg_createWitness(e_poly, crs, 0);
    assert(kzg_verifyEval(e_commitment, crs, 0, e_oracle.first, e_oracle.second));
    assert(T_oracle.first * l * k == e_oracle.first * b);

    //Step4: 
    Fr r_1;
    r_1.setRand();

    std::pair<Fr, G1> a_oracle = kzg_createWitness(a_poly, crs, r_0);
    assert(kzg_verifyEval(a_commitment, crs, r_0, a_oracle.first, a_oracle.second));
    std::pair<Fr, G1> c_oracle = kzg_createWitness(c_poly, crs, r_1);
    assert(kzg_verifyEval(c_commitment, crs, r_1, c_oracle.first, c_oracle.second));
    std::pair<Fr, std::pair<G1, G1>> d_oracle = kzg_createWitness_d2(d_poly, crs, r_0, r_1);
    assert(kzg_verifyEval_d2(d_commitment, crs, r_0, r_1, d_oracle.first, d_oracle.second));
    std::pair<Fr, G1> h_1_oracle = kzg_createWitness(h_1_poly, crs, r_1);
    assert(kzg_verifyEval(h_1_commitment, crs, r_1, h_1_oracle.first, h_1_oracle.second));
    std::pair<Fr, G1> g_1_oracle = kzg_createWitness(g_1_poly, crs, r_1);
    assert(kzg_verifyEval(g_1_commitment, crs, r_1, g_1_oracle.first, g_1_oracle.second));
    Fr Z_H_k_eval = evaluate_polynomial(Z_H_k, r_1);
    assert(c_oracle.first * d_oracle.first - a_oracle.first / k == h_1_oracle.first * Z_H_k_eval + r_1 * g_1_oracle.first);

    T_oracle = kzg_createWitness_d2(T_poly, crs, r_0, r_1);
    assert(kzg_verifyEval_d2(T_commitment, crs, r_0, r_1, T_oracle.first, T_oracle.second));
    std::pair<Fr, std::pair<G1, G1>> p_oracle = kzg_createWitness_d2(p_poly, crs, r_0, r_1);
    assert(kzg_verifyEval_d2(p_commitment, crs, r_0, r_1, p_oracle.first, p_oracle.second));
    std::pair<Fr, std::pair<G1, G1>> q_oracle = kzg_createWitness_d2(q_poly, crs, r_0, r_1);
    assert(kzg_verifyEval_d2(q_commitment, crs, r_0, r_1, q_oracle.first, q_oracle.second));
    Fr Z_H_l_eval = evaluate_polynomial(Z_H_l, r_0);
    assert(T_oracle.first * (alpha + d_oracle.first) - 1 == p_oracle.first * Z_H_l_eval + q_oracle.first * Z_H_k_eval);

    std::cout<<"1"<<std::endl;

    // __int128 uu = 8;
    // subgroup group(uu);
    // Fr generator = group.getGenerator();
    // Fr y[8];
    // y[0] = 0;
    // y[1] = 1;
    // y[2] = 2;
    // y[3] = 3;
    // y[4] = -1;
    // y[5] = 2;
    // y[6] = 4;
    // y[7] = -2;
    // NTT(y, uu, group, 1);
    // Fr cur = 1;
    // std::cout<< group.getGenerator() << std::endl;
    // for (int i = 0; i < uu; i++) {
    //     std::cout<< cur << std::endl;
    //     cur *= generator;
    // }
    // NTT(y, uu, group, -1);
    // for (int i = 0; i < uu; i++) std::cout<< y[i]/uu << std::endl;
    // std::vector<std::vector<Fr>> poly1 = {{0, 1}, {0, 1}, {0, 1}}, poly2 = {{1, 0}, {1, 0}};
    // std::vector<std::vector<Fr>> product = poly_product_d2(poly1, poly2);
    // for (auto i : product) { 
    //     for (auto j : i) std::cout<< j << "  ";
    //     std::cout<< std::endl;
    // }
    // Fr aa = 200;
    // std::cout<< sizeof(size_t)<<std::endl;
    // __int128 xx = (Fr_to_int(R)+1);
    // print_int128(H_l.getOrder());
    // std::cout<< H_l.getGenerator()<<std::endl;
    // for (size_t i=0; i<H_l.order; i++) std::cout<< H_l[i]<<std::endl;
    // __int128 xx = 1152921504606846976;
    // print_int128(xx *16);
    //for (auto i : G) std::cout<< i << std::endl;

    return 0;
=======
#include <iostream>
#include <vector>
#include <mcl/bls12_381.hpp>
#include <random>
#include <cstring>
#include <cmath>
#include "kzg.h"
#include "subgroups.h"

using namespace mcl::bls12;

constexpr int MAX_N = 1 << 6;
Fr tmp1[MAX_N], tmp2[MAX_N]; //for INTT

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

__int128 generate_eff() {
    // 随机数引擎
    std::random_device rd;   // 使用真正的随机数种子
    std::mt19937_64 gen(rd()); // 使用64位的Mersenne Twister引擎
    std::uniform_int_distribution<uint64_t> dis(1, std::numeric_limits<uint32_t>::max());

    // 生成 1 到 2^32 - 1 的随机数并转换为 __int128 类型
    return static_cast<__int128>(dis(gen));
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

std::vector<Fr> generate_polynomial(const __int128 &l) {
    std::vector<Fr> coeffs;
    Fr r;
    for (int i = 0; i < l; i++) {
        r = generate_eff();
        //r = (Fr_to_int(r));
        coeffs.push_back(r);
    }
    return coeffs;
}

void NTT(Fr* f, __int128 n, const subgroup& group, __int128 reverse) {
    if (n == 1) return;
    for (__int128 i = 0; i < n; ++i) tmp1[i] = f[i];
    for (__int128 i = 0; i < n; ++i) {
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
    for (__int128 k = 0; k < n / 2; ++k) {
        tmp1[k] = g[k] + cur * h[k];
        tmp1[k + n / 2] = g[k] - cur * h[k];
        cur *= step;
    }
    for (__int128 i = 0; i < n; ++i) f[i] = tmp1[i];
}

void NTT_d2(Fr (&f)[MAX_N][MAX_N], __int128 n1, __int128 n2, const subgroup& group1, const subgroup& group2, __int128 reverse) {
    for (__int128 i = 0; i < n1; ++i) NTT(f[i], n2, group2, reverse);
    for (__int128 j = 0; j < n2; ++j) {
        for (__int128 i = 0; i < n1; ++i) tmp2[i] = f[i][j];
        NTT(tmp2, n1, group1, reverse);
        for (__int128 i = 0; i < n1; ++i) f[i][j] = tmp2[i];
    }
}

std::vector<Fr> poly_product(const std::vector<Fr>& poly1, const std::vector<Fr>& poly2) {
    __int128 n = 1;
    while (n < poly1.size() + poly2.size()) {
        n <<= 1;
    }
    subgroup group(n);

    Fr y_1[MAX_N], y_2[MAX_N];
    for (__int128 i = 0; i < n; i++) {
        if (i < poly1.size()) y_1[i] = poly1[i];
        else y_1[i] = 0;
        if (i < poly2.size()) y_2[i] = poly2[i];
        else y_2[i] = 0;
    }
    NTT(y_1, n, group, 1);
    NTT(y_2, n, group, 1);
    for (__int128 i = 0; i < n; i++) {
        y_1[i] = y_1[i] * y_2[i];
    }
    NTT(y_1, n, group, -1);
    // NTT(y_2, n, group, -1);
    // for (__int128 i = 0; i < n; i++) {std::cout<< y_1[i] << std::endl;}
    // for (__int128 i = 0; i < n; i++) {std::cout<< y_2[i] << std::endl;}
    std::vector<Fr> result;
    for (__int128 i = 0; i < n; i++) result.push_back(y_1[i]/n);
    return result;
}

std::vector<std::vector<Fr>> poly_product_d2(const std::vector<std::vector<Fr>>& poly1, const std::vector<std::vector<Fr>>& poly2) {
    __int128 n1 = 1, n2 = 1;
    while (n1 < poly1.size() + poly2.size()) {
        n1 <<= 1;
    }
    while (n2 < poly1[0].size() + poly2[0].size()) {
        n2 <<= 1;
    }
    subgroup group1(n1), group2(n2);

    Fr yy_1[MAX_N][MAX_N], yy_2[MAX_N][MAX_N];
    for (__int128 i = 0; i < n1; i++) {
        for (__int128 j = 0; j < n2; j++) {
            if (i < poly1.size() && j < poly1[0].size()) yy_1[i][j] = poly1[i][j];
            else yy_1[i][j] = 0;
            if (i < poly2.size() && j < poly2[0].size()) yy_2[i][j] = poly2[i][j];
            else yy_2[i][j] = 0;
        }
    }
    NTT_d2(yy_1, n1, n2, group1, group2, 1);
    NTT_d2(yy_2, n1, n2, group1, group2, 1);
    for (__int128 i = 0; i < n1; i++) {
        for (__int128 j = 0; j < n2; j++) {
            yy_1[i][j] *= yy_2[i][j];
        }
    }
    NTT_d2(yy_1, n1, n2, group1, group2, -1);
    std::vector<std::vector<Fr>> result;
    for (__int128 i = 0; i < n1; i++) {
        std::vector<Fr> tmp;
        result.push_back(tmp);
        for (__int128 j = 0; j < n2; j++) {
            result[i].push_back(yy_1[i][j] / n1 / n2);
        }
    }
    return result;
}

std::vector<Fr> sumset_representation(const Fr& R, const Fr& b) {
    __int128 R_j = Fr_to_int(R), b_int = Fr_to_int(b);
    std::vector<Fr> G;
    __int128 G_j = (R_j + 1) / b_int;
    while (R_j > b_int) {
        G.push_back(G_j);
        R_j = R_j - (b_int - 1) * G_j;
        G_j = (R_j + 1) / b_int;
    }
    if (G[G.size() - 1] != 1) G.push_back(1);
    return G;
}

int main() {
    initPairing(mcl::BLS12_381);
    G1 g1;
    hashAndMapToG1(g1, "ab", 2);
    G2 g2;
    hashAndMapToG2(g2, "abc", 3);

    //input
    Fr R;
    __int128 b, k, l;
    R.setStr("18446744073709551615"); //2^64 - 1
    b = 16, k = 16, l = 32;
    std::vector<Fr> a = generate_polynomial(l);
    CRS crs = kzg_setup(g1, g2, 2 * l + 2, 2 * k + 2);

    //a[6].setStr("18446744073709551616"); //for testing over range

    std::vector<Fr> G = sumset_representation(R, b);
    subgroup H_l(l), H_k(k), H_b(b);

    // bit decomposition
    std::vector<__int128> GG;
    for (auto i : G) {
        GG.push_back(Fr_to_int(i));
    }
    std::vector<std::vector<Fr>> d;
    for (__int128 i = 0; i < l; i++) {
        std::vector<Fr> tmp;
        d.push_back(tmp);
        __int128 ai = Fr_to_int(a[i]);
        for (__int128 j = 0; j < k; j++) {
            d[i].push_back(ai / GG[j]);
            ai = ai % GG[j];
        }
    }
    std::vector<Fr> m;
    for (__int128 i = 0; i < b; i++) m.push_back(0);
    for (auto i : d) {
        for (auto j : i) {
            __int128 index = Fr_to_int(j);
            if (index >= 0 && index < b) m[index] += 1;
        }
    }
    
    //Step 1:
    Fr y[MAX_N], yy[MAX_N][MAX_N];
    for (__int128 i = 0; i < l; i++) {// d_extension
        for (__int128 j = 0; j < k; j++) yy[i][j] = d[i][j];
    }
    NTT_d2(yy, l, k, H_l, H_k, -1);
    std::vector<std::vector<Fr>> d_poly;
    for (__int128 i = 0; i < l; i++) {
        std::vector<Fr> vec;
        d_poly.push_back(vec);
        for (__int128 j = 0; j < k; j++) {
            d_poly[i].push_back(yy[i][j]/l/k);
        }
    }

    for (__int128 i = 0; i < l; i++) {// a_extension
        y[i] = a[i];
    }
    NTT(y, l, H_l, -1);
    std::vector<Fr> a_poly;
    for (__int128 i = 0; i < l; i++) {
        a_poly.push_back(y[i]/l);
    }

    for (__int128 i = 0; i < k; i++) {// c_extension
        y[i] = G[i];
    }
    NTT(y, k, H_k, -1);
    std::vector<Fr> c_poly;
    for (__int128 i = 0; i < k; i++) {
        c_poly.push_back(y[i]/k);
    }

    for (__int128 i = 0; i < b; i++) {// I_extension
        y[i] = i;
    }
    NTT(y, b, H_b, -1);
    std::vector<Fr> I_poly;
    for (__int128 i = 0; i < b; i++) {
        I_poly.push_back(y[i] / b);
    }

    for (__int128 i = 0; i < b; i++) {// m_extension
        y[i] = m[i];
    }
    NTT(y, b, H_b, -1);
    std::vector<Fr> m_poly;
    for (__int128 i = 0; i < b; i++) {
        m_poly.push_back(y[i] / b);
    }
    G1 d_commitment = kzg_commit_d2(d_poly, crs);
    assert(verify_kzg_open_d2(d_commitment, d_poly, crs));
    G1 a_commitment = kzg_commit(a_poly, crs);
    assert(verify_kzg_open(a_commitment, a_poly, crs));
    G1 c_commitment = kzg_commit(c_poly, crs);
    assert(verify_kzg_open(c_commitment, c_poly, crs));
    G1 I_commitment = kzg_commit(I_poly, crs);
    assert(verify_kzg_open(I_commitment, I_poly, crs));
    G1 m_commitment = kzg_commit(m_poly, crs);
    assert(verify_kzg_open(m_commitment, m_poly, crs));

    //Step2:
    Fr alpha, r_0;
    alpha.setRand();
    r_0.setRand();
    std::vector<std::vector<Fr>> T;
    for (__int128 i = 0; i < l; i++) {
        std::vector<Fr> tmp;
        T.push_back(tmp);
        for (__int128 j = 0; j < k; j++) {
            T[i].push_back(1 / (alpha + d[i][j]));
        }
    }
    std::vector<Fr> e;
    for (__int128 i = 0; i < b; i++) e.push_back(m[i] / (alpha + i));
    for (__int128 i = 0; i < l; i++) {// T_extension
        for (__int128 j = 0; j < k; j++) yy[i][j] = T[i][j];
    }
    NTT_d2(yy, l, k, H_l, H_k, -1);
    std::vector<std::vector<Fr>> T_poly;
    for (__int128 i = 0; i < l; i++) {
        std::vector<Fr> vec;
        T_poly.push_back(vec);
        for (__int128 j = 0; j < k; j++) {
            T_poly[i].push_back(yy[i][j]/l/k);
        }
    }

    for (__int128 i = 0; i < b; i++) {// e_extension
        y[i] = e[i];
    }
    NTT(y, b, H_b, -1);
    std::vector<Fr> e_poly;
    for (__int128 i = 0; i < b; i++) {
        e_poly.push_back(y[i]/b);
    }
    G1 T_commitment = kzg_commit_d2(T_poly, crs);
    assert(verify_kzg_open_d2(T_commitment, T_poly, crs));
    G1 e_commitment = kzg_commit(e_poly, crs);
    assert(verify_kzg_open(e_commitment, e_poly, crs));

    //Step3:
    std::vector<Fr> F1;
    std::vector<Fr> d_r0_poly;
    for (__int128 j = 0; j < d_poly[0].size(); j++) {
        std::vector<Fr> tmp;
        for (__int128 i = 0; i < d_poly.size(); i++) tmp.push_back(d_poly[i][j]);
        d_r0_poly.push_back(evaluate_polynomial(tmp, r_0));
    }
    F1 = poly_product(c_poly, d_r0_poly);
    F1[0] -= evaluate_polynomial(a_poly, r_0) / k;
    std::vector<Fr> Z_H_k(k + 1, 0);
    Z_H_k[0] = -1, Z_H_k[k] = 1;
    std::pair<std::vector<Fr>, std::vector<Fr>> tmp_pair = polynomial_division(F1, Z_H_k);
    std::vector<Fr> h_1_poly = tmp_pair.first;
    std::vector<Fr> g_1_poly = tmp_pair.second;
    assert(g_1_poly[0] == 0); //check remainder have no constant when divide F1
    for (__int128 i = 0; i < g_1_poly.size() - 1; i++) g_1_poly[i] = g_1_poly[i + 1];
    g_1_poly[g_1_poly.size() - 1] = 0;

    std::vector<std::vector<Fr>> F2 = d_poly;
    F2[0][0] += alpha;
    F2 = poly_product_d2(F2, T_poly);
    F2[0][0] -= 1;
    std::vector<Fr> Z_H_l(l + 1, 0);
    Z_H_l[0] = -1, Z_H_l[l] = 1;
    std::pair<std::vector<std::vector<Fr>>, std::vector<std::vector<Fr>>> tmp_pair2 = polynomial_division_d2(F2, Z_H_l, Z_H_k);
    std::vector<std::vector<Fr>> p_poly = tmp_pair2.first;
    std::vector<std::vector<Fr>> q_poly = tmp_pair2.second;

    std::vector<Fr> F3 = I_poly;
    F3[0] += alpha;
    F3 = poly_product(F3, e_poly);
    for (__int128 i = 0; i < b; i++) {
        F3[i] -= m_poly[i];
    }
    std::vector<Fr> Z_H_b(b + 1, 0);
    Z_H_b[0] = -1, Z_H_b[b] = 1;
    std::pair<std::vector<Fr>, std::vector<Fr>> tmp_pair3 = polynomial_division(F3, Z_H_b);
    for (auto i : tmp_pair3.second) assert(i == 0); //check no remainder when divide F3
    std::vector<Fr> h_2_poly = tmp_pair3.first;
    G1 h_1_commitment = kzg_commit(h_1_poly, crs);
    assert(verify_kzg_open(h_1_commitment, h_1_poly, crs));
    G1 g_1_commitment = kzg_commit(g_1_poly, crs);
    assert(verify_kzg_open(g_1_commitment, g_1_poly, crs));
    G1 p_commitment = kzg_commit_d2(p_poly, crs);
    assert(verify_kzg_open_d2(p_commitment, p_poly, crs));
    G1 q_commitment = kzg_commit_d2(q_poly, crs);
    assert(verify_kzg_open_d2(q_commitment, q_poly, crs));
    G1 h_2_commitment = kzg_commit(h_2_poly, crs);
    assert(verify_kzg_open(h_2_commitment, h_2_poly, crs));

    //Step3 check:
    std::pair<Fr, G1> e_oracle = kzg_createWitness(e_poly, crs, r_0);
    assert(kzg_verifyEval(e_commitment, crs, r_0, e_oracle.first, e_oracle.second));
    std::pair<Fr, G1> I_oracle = kzg_createWitness(I_poly, crs, r_0);
    assert(kzg_verifyEval(I_commitment, crs, r_0, I_oracle.first, I_oracle.second));
    std::pair<Fr, G1> m_oracle = kzg_createWitness(m_poly, crs, r_0);
    assert(kzg_verifyEval(m_commitment, crs, r_0, m_oracle.first, m_oracle.second));
    std::pair<Fr, G1> h_2_oracle = kzg_createWitness(h_2_poly, crs, r_0);
    assert(kzg_verifyEval(h_2_commitment, crs, r_0, h_2_oracle.first, h_2_oracle.second));
    Fr Z_H_b_eval = evaluate_polynomial(Z_H_b, r_0);
    assert(e_oracle.first * (alpha + I_oracle.first) - m_oracle.first == h_2_oracle.first * Z_H_b_eval);

    std::pair<Fr, std::pair<G1, G1>> T_oracle = kzg_createWitness_d2(T_poly, crs, 0, 0);
    assert(kzg_verifyEval_d2(T_commitment, crs, 0, 0, T_oracle.first, T_oracle.second));
    e_oracle = kzg_createWitness(e_poly, crs, 0);
    assert(kzg_verifyEval(e_commitment, crs, 0, e_oracle.first, e_oracle.second));
    assert(T_oracle.first * l * k == e_oracle.first * b);

    //Step4: 
    Fr r_1;
    r_1.setRand();

    std::pair<Fr, G1> a_oracle = kzg_createWitness(a_poly, crs, r_0);
    assert(kzg_verifyEval(a_commitment, crs, r_0, a_oracle.first, a_oracle.second));
    std::pair<Fr, G1> c_oracle = kzg_createWitness(c_poly, crs, r_1);
    assert(kzg_verifyEval(c_commitment, crs, r_1, c_oracle.first, c_oracle.second));
    std::pair<Fr, std::pair<G1, G1>> d_oracle = kzg_createWitness_d2(d_poly, crs, r_0, r_1);
    assert(kzg_verifyEval_d2(d_commitment, crs, r_0, r_1, d_oracle.first, d_oracle.second));
    std::pair<Fr, G1> h_1_oracle = kzg_createWitness(h_1_poly, crs, r_1);
    assert(kzg_verifyEval(h_1_commitment, crs, r_1, h_1_oracle.first, h_1_oracle.second));
    std::pair<Fr, G1> g_1_oracle = kzg_createWitness(g_1_poly, crs, r_1);
    assert(kzg_verifyEval(g_1_commitment, crs, r_1, g_1_oracle.first, g_1_oracle.second));
    Fr Z_H_k_eval = evaluate_polynomial(Z_H_k, r_1);
    assert(c_oracle.first * d_oracle.first - a_oracle.first / k == h_1_oracle.first * Z_H_k_eval + r_1 * g_1_oracle.first);

    T_oracle = kzg_createWitness_d2(T_poly, crs, r_0, r_1);
    assert(kzg_verifyEval_d2(T_commitment, crs, r_0, r_1, T_oracle.first, T_oracle.second));
    std::pair<Fr, std::pair<G1, G1>> p_oracle = kzg_createWitness_d2(p_poly, crs, r_0, r_1);
    assert(kzg_verifyEval_d2(p_commitment, crs, r_0, r_1, p_oracle.first, p_oracle.second));
    std::pair<Fr, std::pair<G1, G1>> q_oracle = kzg_createWitness_d2(q_poly, crs, r_0, r_1);
    assert(kzg_verifyEval_d2(q_commitment, crs, r_0, r_1, q_oracle.first, q_oracle.second));
    Fr Z_H_l_eval = evaluate_polynomial(Z_H_l, r_0);
    assert(T_oracle.first * (alpha + d_oracle.first) - 1 == p_oracle.first * Z_H_l_eval + q_oracle.first * Z_H_k_eval);

    std::cout<<"1"<<std::endl;

    // __int128 uu = 8;
    // subgroup group(uu);
    // Fr generator = group.getGenerator();
    // Fr y[8];
    // y[0] = 0;
    // y[1] = 1;
    // y[2] = 2;
    // y[3] = 3;
    // y[4] = -1;
    // y[5] = 2;
    // y[6] = 4;
    // y[7] = -2;
    // NTT(y, uu, group, 1);
    // Fr cur = 1;
    // std::cout<< group.getGenerator() << std::endl;
    // for (int i = 0; i < uu; i++) {
    //     std::cout<< cur << std::endl;
    //     cur *= generator;
    // }
    // NTT(y, uu, group, -1);
    // for (int i = 0; i < uu; i++) std::cout<< y[i]/uu << std::endl;
    // std::vector<std::vector<Fr>> poly1 = {{0, 1}, {0, 1}, {0, 1}}, poly2 = {{1, 0}, {1, 0}};
    // std::vector<std::vector<Fr>> product = poly_product_d2(poly1, poly2);
    // for (auto i : product) { 
    //     for (auto j : i) std::cout<< j << "  ";
    //     std::cout<< std::endl;
    // }
    // Fr aa = 200;
    // std::cout<< sizeof(size_t)<<std::endl;
    // __int128 xx = (Fr_to_int(R)+1);
    // print_int128(H_l.getOrder());
    // std::cout<< H_l.getGenerator()<<std::endl;
    // for (size_t i=0; i<H_l.order; i++) std::cout<< H_l[i]<<std::endl;
    // __int128 xx = 1152921504606846976;
    // print_int128(xx *16);
    //for (auto i : G) std::cout<< i << std::endl;

    return 0;
>>>>>>> ea97283b7d0d33a96163f269cae4f9e92a8490e0
}
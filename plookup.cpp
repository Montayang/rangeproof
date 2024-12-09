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

constexpr int n = 1 << 20;

std::vector<int> batch_bit_reverse(int log_n) {
    int n = 1 << log_n;
    std::vector<int> res(n, 0);
    for (int i = 0; i < n; ++i) {
        res[i] = (res[i >> 1] >> 1) | ((i & 1) << (log_n - 1));
    }
    return res;
}

void NTT(std::vector<Fr>& coeff, const subgroup& group, int reverse) {
    int n = coeff.size();
    int log_n = static_cast<int>(std::log2(n));
    auto rank = batch_bit_reverse(log_n);

    for (int i = 0; i < n; ++i) {
        if (i < rank[i]) {
            std::swap(coeff[i], coeff[rank[i]]);
        }
    }

    int log_m = 0;
    Fr step = group.getGenerator();
    if (reverse == -1) step = 1 / step;
    for (int i = 0; i < log_n; ++i) {
        Fr w_m;
        Fr::pow(w_m, step, n >> (log_m + 1));
        int m = 1 << log_m;
        for (int j = 0; j < n; j += m * 2) {
            Fr w = 1;
            for (int k = 0; k < m; ++k) {
                Fr t = w * coeff[j + k + m];
                coeff[j + k + m] = coeff[j + k] - t;
                coeff[j + k] += t;
                w *= w_m;
            }
        }
        ++log_m;
    }
    if (reverse == -1) {
        for (int i = 0; i < coeff.size(); i++) coeff[i] /= n;
    }
}

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

std::string generate_eff() {
    // 使用std::random_device和std::mt19937_64生成随机数
    std::random_device rd;  // 用于生成随机种子
    std::mt19937_64 gen(rd());  // 使用64位的Mersenne Twister引擎
    std::uniform_int_distribution<uint64_t> dist(0, n - 1);  // 分布范围为[0, 2^20 - 1]

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
        r.setStr(generate_eff());
        coeffs.push_back(r);
    }
    return coeffs;
}

std::vector<Fr> poly_product(const std::vector<Fr>& poly1, const std::vector<Fr>& poly2) {
    int n = 1;
    while (n < poly1.size() + poly2.size()) {
        n <<= 1;
    }
    subgroup group(n);

    std::vector<Fr> y_1(n, 0), y_2(n, 0);
    for (int i = 0; i < poly1.size(); i++) y_1[i] = poly1[i];
    for (int i = 0; i < poly2.size(); i++) y_2[i] = poly2[i];
    NTT(y_1, group, 1);
    NTT(y_2, group, 1);

    for (int i = 0; i < n; i++) {
        y_1[i] = y_1[i] * y_2[i];
    }
    NTT(y_1, group, -1);
    // NTT(y_2, n, group, -1);
    // for (int i = 0; i < n; i++) {std::cout<< y_1[i] << std::endl;}
    // for (int i = 0; i < n; i++) {std::cout<< y_2[i] << std::endl;}
    return y_1;
}

std::vector<Fr> poly_product_force(const std::vector<Fr>& poly1, const std::vector<Fr>& poly2) {
    size_t n1 = poly1.size();
    size_t n2 = poly2.size();
    size_t result_size = n1 + n2 - 1;

    // 初始化结果多项式的系数为0
    std::vector<Fr> result(result_size, 0);

    // 暴力计算多项式乘积
    for (size_t i = 0; i < n1; ++i) {
        for (size_t j = 0; j < n2; ++j) {
            result[i + j] += poly1[i] * poly2[j];
        }
    }

    return result;
}

// void preprocessLagrange(const std::vector<Fr>& x_points, std::vector<Fr>& prefix_product, std::vector<Fr>& suffix_product) {
//     int n = x_points.size();
//     prefix_product.resize(n);
//     suffix_product.resize(n);

//     // 初始化前缀积
//     prefix_product[0] = 1;
//     for (int i = 1; i < n; ++i) {
//         prefix_product[i] = prefix_product[i - 1] * (x_points[i - 1]);
//     }

//     // 初始化后缀积
//     suffix_product[n - 1] = 1;
//     for (int i = n - 2; i >= 0; --i) {
//         suffix_product[i] = suffix_product[i + 1] * (x_points[i + 1]);
//     }
// }

// Fr computeLagrangeAt(const std::vector<Fr>& x_points, const std::vector<Fr>& prefix_product,
//                          const std::vector<Fr>& suffix_product, int i, Fr r0) {
//     int n = x_points.size();

//     // 计算分子
//     Fr numerator = (i == 0 ? 1 : prefix_product[i]) * (i == n - 1 ? 1 : suffix_product[i]);

//     // 计算分母
//     Fr denominator = 1;
//     for (int j = 0; j < n; ++j) {
//         if (j != i) {
//             denominator *= (x_points[i] - x_points[j]);
//         }
//     }

//     return numerator / denominator;
// }

int main() {
    initPairing(mcl::BN_SNARK1);
    G1 g1;
    hashAndMapToG1(g1, "ab", 2);
    G2 g2;
    hashAndMapToG2(g2, "abc", 3);
    std::chrono::duration<double> prover_time(0);
    std::chrono::duration<double> verifier_time(0);
        std::chrono::duration<double> test_time(0);
    size_t proof_size = 0;

    //input
    Fr N = n;
    int m = 16384;
    std::vector<Fr> t;
    Fr cur = 0;
    for (int i = 0; i < n; i++) {
        t.push_back(cur);
        cur += 1;
    }
    std::vector<Fr> f = generate_polynomial(m);
        //f[6].setStr("1048577"); //for testing over range
    CRS crs = kzg_setup(g1, g2, 2 * n + 3, 0);
    subgroup H(n);
    Fr g = H.getGenerator();

    //construct L_1 and L_0
    std::vector<Fr> y_1(n, 0), y_2(n, 0);//for ntt
    std::vector<Fr> L_1, L_0;
    for (int i = 1; i < n; i++) {
        y_1[i] = 0;
    }
    y_1[1] = 1;
    NTT(y_1, H, -1);
    for (int i = 0; i < n; i++) {
        L_1.push_back(y_1[i]/n);
    }
    for (int i = 1; i < n; i++) {
        y_1[i] = 0;
    }
    y_1[0] = 1;
    NTT(y_1, H, -1);
    for (int i = 0; i < n; i++) {
        L_0.push_back(y_1[i]);
    }
    
    // std::vector<Fr> x_points;
    // Fr cur_x = 1;
    // for (int i = 0; i < n; i++) {
    //     x_points.push_back(cur_x);
    //     cur_x *= g;
    // }
    // std::vector<Fr> prefix_product, suffix_product;
    // preprocessLagrange(x_points, prefix_product, suffix_product);

    std::cout << "Setup Completed." << std::endl;
    //construt s
    for (int i = m; i < n; i++) f.push_back(n - 1);
    std::vector<Fr> s = t;
    for (auto i : f) s.push_back(i);
    std::sort(s.begin(), s.end());

    // auto start = std::chrono::high_resolution_clock::now();
    // kzg_commit(L_1, crs);
    // auto end = std::chrono::high_resolution_clock::now();
    // test_time = end - start;
    // std::cout << test_time.count() << std::endl;
    
    //core of protocol
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<Fr> t_poly, f_poly;
    for (int i = 0; i < n; i++) {// construct t
        y_1[i] = t[i];
    }
    NTT(y_1, H, -1);
    for (int i = 0; i < n; i++) {
        t_poly.push_back(y_1[i]);
    }

    for (int i = 0; i < n; i++) {// construct f
        y_1[i] = f[i];
    }
    NTT(y_1, H, -1);
    for (int i = 0; i < n; i++) {
        f_poly.push_back(y_1[i]);
    }

    std::vector<Fr> h_1_poly, h_2_poly;
    for (int i = 1; i < n; i++) {// construct h_1
        y_1[i] = s[i];
    }
    y_1[0] = s[n];
    NTT(y_1, H, -1);
    for (int i = 0; i < n; i++) {
        h_1_poly.push_back(y_1[i]);
    }

    for (int i = 1; i < n; i++) {// construct h_2
        y_1[i] = s[i + n - 1];
    }
    y_1[0] = s[2 * n - 1];
    NTT(y_1, H, -1);
    for (int i = 0; i < n; i++) {
        h_2_poly.push_back(y_1[i]);
    }
    
    // G1 L_1_commitment = kzg_commit(L_1, crs);
    // G1 L_n_commitment = kzg_commit(L_n, crs);
    auto end = std::chrono::high_resolution_clock::now();
    prover_time += end - start;
    G1 t_commitment = kzg_commit(t_poly, crs);
    G1 f_commitment = kzg_commit(f_poly, crs);
    start = std::chrono::high_resolution_clock::now();
    G1 h_1_commitment = kzg_commit(h_1_poly, crs);
    G1 h_2_commitment = kzg_commit(h_2_poly, crs);
    end = std::chrono::high_resolution_clock::now();
    prover_time += end - start;
        test_time = end - start;
        std::cout << "input commited  " << test_time.count() << std::endl;

    proof_size += sizeof(f_commitment) + sizeof(h_1_commitment) + sizeof(h_2_commitment);

    //verifier construct beta, gamma
    Fr beta, gamma;
    beta.setRand();
    gamma.setRand();

    //construct Z
    start = std::chrono::high_resolution_clock::now();
    std::vector<Fr> Z_poly;
    y_1[0] = 1;
    y_1[1] = 1;
    for (int i = 2; i < n; i++) {
        // Fr Z_eval;
        // Fr::pow(Z_eval, 1 + beta, i - 1);
        // for (int l = 1; l < i; l++) {
        //     Z_eval *= ((gamma + f[l]) * (gamma * (1 + beta) + t[l] + beta * t[l + 1]));
        //     Z_eval /= ((gamma * (1 + beta) + s[l] + beta * s[l + 1]) * (gamma * (1 + beta) + s[n + l - 1] + beta * s[n + l]));
        // }
        y_1[i] = y_1[i-1] * (1 + beta) * ((gamma + f[i - 1]) * (gamma * (1 + beta) + t[i - 1] + beta * t[i]));
        y_1[i] /= ((gamma * (1 + beta) + s[i - 1] + beta * s[i]) * (gamma * (1 + beta) + s[n + i - 2] + beta * s[n + i - 1]));
    }
        // cur = 1;
        // std::vector <Fr> lhs(1, 0), rhs(1, 0);
        // for (int i = 1; i < n - 1; i++) {
        //     lhs.push_back(y_1[i]*(1 + beta)*(gamma+f[i])*(gamma*(1+beta)+t[i]+beta*t[i+1]));
        //     cur*=g;
        // }
    NTT(y_1, H, -1);
    for (int i = 0; i < n; i++) {
        Z_poly.push_back(y_1[i]);
    }
        
    G1 Z_commitment = kzg_commit(Z_poly, crs);
    end = std::chrono::high_resolution_clock::now();
    prover_time += end - start;
        test_time = end - start;
        std::cout << "Z commited  " << test_time.count() << std::endl;

    start = std::chrono::high_resolution_clock::now();
    //prepare for zerotest
    std::vector<Fr> F1, F2, F3, F4, g_4_poly(2 * n, 0);
    F1 = Z_poly;
    F1[0] -= 1;
    F1 = poly_product(L_1, F1);
    //for (int i = n; i < 2 * n - 1; i++) g_1_poly.push_back(F1[i]); //construct g_1

    F2 = Z_poly;
    F2[0] -= 1;
    F2 = poly_product(L_0, F2);
    //for (int i = n; i < 2 * n - 1; i++) g_2_poly.push_back(F2[i]); //construct g_2

    F3 = h_1_poly;
    cur = 1;
    for (int i = 0; i < n; i++) {
        F3[i] -= h_2_poly[i] * cur;
        cur *= g;
    }
    F3 = poly_product(L_0, F3);
    //for (int i = n; i < 2 * n - 1; i++) g_3_poly.push_back(F3[i]); //construct g_3

    std::vector<Fr> tmp1_F4, tmp2_F4, tmp3_F4;
    tmp1_F4 = t_poly;
    tmp1_F4[0] += gamma * (1 + beta);
    cur = 1;
    for (int i = 0; i < n; i++) {
        tmp1_F4[i] += t_poly[i] * cur * beta;
        cur *= g;
    }
    tmp2_F4 = f_poly;
    tmp2_F4[0] += gamma;
    tmp1_F4 = poly_product(tmp2_F4, tmp1_F4);
    tmp1_F4 = poly_product(Z_poly, tmp1_F4);
    tmp1_F4 = poly_product_force({1 + beta, -1 - beta - g - beta * g, g + beta * g}, tmp1_F4);
    
    tmp2_F4 = h_2_poly;
    tmp2_F4[0] += gamma * (1 + beta);
    cur = 1;
    for (int i = 0; i < n; i++) {
        tmp2_F4[i] += h_2_poly[i] * cur * beta;
        cur *= g;
    }
    tmp3_F4 = h_1_poly;
    tmp3_F4[0] += gamma * (1 + beta);
    cur = 1;
    for (int i = 0; i < n; i++) {
        tmp3_F4[i] += h_1_poly[i] * cur * beta;
        cur *= g;
    }
    tmp2_F4 = poly_product(tmp2_F4, tmp3_F4);
    cur = 1;
    for (int i = 0; i < n; i++) {
        tmp3_F4[i] = Z_poly[i] * cur;
        cur *= g;
    }
    tmp2_F4 = poly_product(tmp2_F4, tmp3_F4);
    tmp2_F4 = poly_product_force({1, -1 - g, g}, tmp2_F4);
    for (int i = 0; i < 3 * n; i++) F4.push_back(tmp1_F4[i] - tmp2_F4[i]);

    //combine F1, F2, F3, F4
    Fr alpha_1, alpha_2, alpha_3, alpha_4;
    alpha_1.setRand();
    alpha_2.setRand();
    alpha_3.setRand();
    alpha_4.setRand();
    for (int i = 0; i < 2 * n - 1; i++) F4[i] = F1[i] * alpha_1 + F2[i] * alpha_2 + F3[i] * alpha_3 + F4[i] * alpha_4;
    for (int i = 2 * n - 1; i < 3 * n; i++) F4[i] *= alpha_4;

    //construct g_4
    for (int i = 2 * n; i < 3 * n; i++) { 
        g_4_poly[i - n] = F4[i];
    }
    for (int i = 0; i < n; i++) { 
        g_4_poly[i] = -F4[i];
    }
        // cur = 1;
        // for (int i = 0; i < n;i++) {
        //     std::cout<<evaluate_polynomial(F4, cur) - evaluate_polynomial(g_4_poly, cur) * (cur - 1)<<std::endl;
        //     cur*=g;
        // }
        // cur = 1;
        // for (int i = 0; i < n;i++) {
        //     std::cout<<evaluate_polynomial(tmp2_F4, cur) - lhs[i]<<std::endl;
        //     cur*=g;
        // }
    end = std::chrono::high_resolution_clock::now();
    prover_time += end - start;
        test_time = end - start;
        std::cout << "F constructed  " << test_time.count() << std::endl;

    start = std::chrono::high_resolution_clock::now();

    // G1 g_1_commitment = kzg_commit(g_1_poly, crs);
    // G1 g_2_commitment = kzg_commit(g_2_poly, crs);
    // G1 g_3_commitment = kzg_commit(g_3_poly, crs);
    G1 g_4_commitment = kzg_commit(g_4_poly, crs);

    end = std::chrono::high_resolution_clock::now();
    prover_time += end - start;
        test_time = end - start;
        std::cout << "g commited  " << test_time.count() << std::endl;

    proof_size += sizeof(Z_commitment) + sizeof(g_4_commitment);

    //Verify Phase
    Fr r_0;
    r_0.setRand();

    Fr L_1_r0 = evaluate_polynomial(L_1, r_0);
    Fr L_0_r0 = evaluate_polynomial(L_0, r_0);
    
    start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<Fr>> funcs1;
    // funcs1.push_back(L_1);
    // funcs1.push_back(L_0);
    funcs1.push_back(t_poly);
    funcs1.push_back(f_poly);
    funcs1.push_back(h_1_poly);
    funcs1.push_back(h_2_poly);
    funcs1.push_back(Z_poly);
    // funcs1.push_back(g_1_poly);
    // funcs1.push_back(g_2_poly);
    // funcs1.push_back(g_3_poly);
    funcs1.push_back(g_4_poly);
    Batching_witness batching1 = kzg_createWitness_batching(funcs1, crs, r_0);

    std::vector<std::vector<Fr>> funcs2;
    funcs2.push_back(t_poly);
    funcs2.push_back(h_1_poly);
    funcs2.push_back(h_2_poly);
    funcs2.push_back(Z_poly);
    Batching_witness batching2 = kzg_createWitness_batching(funcs2, crs, g * r_0);

    end = std::chrono::high_resolution_clock::now();
    prover_time += end - start;
        test_time = end - start;
        std::cout << "opened  " << test_time.count() << std::endl;

    std::vector<Fr> Z_H(n + 1, 0);
    Z_H[0] = -1, Z_H[n] = 1;
    Fr Z_H_eval = evaluate_polynomial(Z_H, r_0);
    start = std::chrono::high_resolution_clock::now();
    std::vector<G1> commitments1;
    // commitments1.push_back(L_1_commitment);
    // commitments1.push_back(L_n_commitment);
    commitments1.push_back(t_commitment);
    commitments1.push_back(f_commitment);
    commitments1.push_back(h_1_commitment);
    commitments1.push_back(h_2_commitment);
    commitments1.push_back(Z_commitment);
    // commitments1.push_back(g_1_commitment);
    // commitments1.push_back(g_2_commitment);
    // commitments1.push_back(g_3_commitment);
    commitments1.push_back(g_4_commitment);
    assert(kzg_verifyEval_batching(commitments1, crs, r_0, batching1));

    std::vector<G1> commitments2;
    commitments2.push_back(t_commitment);
    commitments2.push_back(h_1_commitment);
    commitments2.push_back(h_2_commitment);
    commitments2.push_back(Z_commitment);
    assert(kzg_verifyEval_batching(commitments2, crs, g * r_0, batching2));
    
    // assert(batching1.evals[0] * (batching1.evals[6] - 1) == batching1.evals[7] * Z_H_eval); //F1
    // assert(batching1.evals[1] * (batching1.evals[6] - 1) == batching1.evals[8] * Z_H_eval); //F2
    // assert(batching1.evals[1] * (batching1.evals[4] - batching2.evals[2]) == batching1.evals[9] * Z_H_eval); //F3
    assert(alpha_1 * L_1_r0 * (batching1.evals[4] - 1) + alpha_2 * L_0_r0 * (batching1.evals[4] - 1) + alpha_3 * L_0_r0 * (batching1.evals[2] - batching2.evals[2]) +
                alpha_4 * ((r_0 - 1) * (g * r_0 - 1) * batching1.evals[4] * (1 + beta) * (gamma + batching1.evals[1]) * (gamma * (1 + beta) + batching1.evals[0] + beta * batching2.evals[0])
                -  (r_0 - 1) * (g * r_0 - 1) * batching2.evals[3] * (gamma * (1 + beta) + batching1.evals[2] + beta * batching2.evals[1]) * (gamma * (1 + beta) + batching1.evals[3] + beta * batching2.evals[2]))
                == batching1.evals[5] * Z_H_eval); //F4

    end = std::chrono::high_resolution_clock::now();
    verifier_time += end - start;

    proof_size += sizeof(batching1) + sizeof(batching2);

    std::cout << "1" << std::endl;

    std::cout << "Prover Time:" << prover_time.count() << std::endl;
    std::cout << "Verifier Time:" << verifier_time.count() << std::endl;
    std::cout << "Proof Size:" << ((double) proof_size) / 1024.0 << " KB" << std::endl;
    return 0;
}
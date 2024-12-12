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

constexpr int n = 64;
constexpr __int128 R = (__int128)1 << n;

std::chrono::duration<double> prover_time(0);
std::chrono::duration<double> verifier_time(0);
size_t proof_size = 0;

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

Fr vec_product(std::vector<Fr> vec1, std::vector<Fr> vec2) {
    Fr res = 0;
    int size = vec1.size();
    for (int i = 0; i < size; i++) res += vec1[i] * vec2[i];
    return res;
}

Fr Fr_pow(const Fr &x, const Fr &y) {
    Fr tmp;
    Fr::pow(tmp, x , y);
    return tmp;
}

bool Inner_Product_Argument(const std::vector<G1> &g, const std::vector<G1> &h, const G1 &u, const G1 &P, const std::vector<Fr> &a, const std::vector<Fr> &b) {
    int len = g.size();
    if (len == 1) {
        proof_size += sizeof(a[0]) + sizeof(b[0]);
        auto start = std::chrono::high_resolution_clock::now();
        Fr c = a[0] * b[0];
        return (P == g[0] * a[0] + h[0] * b[0] + u * c);
        auto end = std::chrono::high_resolution_clock::now();
        verifier_time += end - start;
    }

    auto start = std::chrono::high_resolution_clock::now();
    int len_ = len / 2;
    Fr c_L = 0, c_R = 0;
    for (int i = 0; i < len_; i++) {
        c_L += a[i] * b[len_ + i];
        c_R += a[len_ + i] * b[i];
    }
    G1 L, R;

    L = u * c_L;
    R = u * c_R;
    for (int i = 0; i < len_; i++) {
        L += g[len_ + i] * a[i] + h[i] * b[len_ + i];
        R += g[i] * a[len_ + i] + h[len_ + i] * b[i];
    }
    auto end = std::chrono::high_resolution_clock::now();
    prover_time += end - start;

    proof_size += 2 * G1::getSerializedByteSize();//L, R

    Fr x;
    x.setRand();

    start = std::chrono::high_resolution_clock::now();
    std::vector<G1> g_(len_), h_(len_);
    for (int i = 0; i < len_; i++) {
        g_[i] = g[i] * (1 / x) + g[len_ + i] * x;
        h_[i] = h[i] * x + h[len_ + i] * (1 / x);
    }
    G1 P_ = L * (x * x) + P + R * (1 / x / x);
    end = std::chrono::high_resolution_clock::now();
    prover_time += end - start;
    verifier_time += end - start;

    start = std::chrono::high_resolution_clock::now();
    std::vector<Fr> a_(len_), b_(len_);
    for (int i = 0; i < len_; i++) {
        a_[i] = a[i] * x + a[len_ + i] * (1 / x);
        b_[i] = b[i] * (1 / x) + b[len_ + i] * x;
    }
    end = std::chrono::high_resolution_clock::now();
    prover_time += end - start;
    
    proof_size += sizeof(a_) + sizeof(b_);

    return Inner_Product_Argument(g_, h_, u, P_, a_, b_);
}

int main() {
    initPairing(mcl::BN_SNARK1);
    G1 g;
    hashAndMapToG1(g, "ab", 2);
        
    //input
    int l = 16384;
    std::vector<Fr> v;
    v = generate_polynomial(l);
        //v[6].setStr("18446744073709551616"); //for testing over range
        //std::cout << sizeof(v) <<"  "<< sizeof(g)<<"  "<<G1::getSerializedByteSize()<<std::endl;
    Fr seed;
    G1 h;
    seed.setRand();
    h = g * seed;
    std::vector<G1> g_vec, h_vec;
    for (int i = 0; i < n * l; i++) {
        seed.setRand();
        g_vec.push_back(g * seed);
        seed.setRand();
        h_vec.push_back(g * seed);
    }
    std::vector<G1> V(l);
    std::vector<Fr> gamma(l);
    for (int j = 0; j < l; j++) {
        gamma[j].setRand();
        V[j] = g * v[j] + h * gamma[j];
    }
    
    //bit decomposition
    std::vector<__int128> G;
    __int128 base = 1;
    for (int i = 0; i < n; i++) {
        G.push_back(base);
        base *= 2;
    }
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<Fr> a_L(n * l), a_R(n * l);
    for (int j = 0; j < l; j++) {
        __int128 v_int = Fr_to_int(v[j]);
        for (int i = n - 1; i >= 0; i--) {
            a_L[j * n + i] = v_int / G[i];
            v_int = v_int % G[i];
        }
    }
    for (int i = 0; i < n * l; i++) a_R[i] = a_L[i] - 1;

    Fr alpha;
    alpha.setRand();
    G1 A;
    A = h * alpha;
    for (int i = 0; i < n * l; i++) {
        A += g_vec[i] * a_L[i] + h_vec[i] * a_R[i];
    }
    
    std::vector<Fr> s_L, s_R;
    for (int i = 0; i < n * l; i++) {
        Fr s_L_seed, s_R_seed;
        s_L_seed.setRand();
        s_R_seed.setRand();
        s_L.push_back(s_L_seed);
        s_R.push_back(s_R_seed);
    }

    Fr rou;
    rou.setRand();
    G1 S;
    S = h * rou;
    for (int i = 0; i < n * l; i++) {
        S += g_vec[i] * s_L[i] + h_vec[i] * s_R[i];
    }
    auto end = std::chrono::high_resolution_clock::now();
    prover_time += end - start;

    proof_size += 2 * G1::getSerializedByteSize(); //A, S
    //A, S sent
    Fr y, z;
    y.setRand();
    z.setRand();

    start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<Fr>> l_poly, r_poly; //vector polynomial
    std::vector<Fr> t_poly;
    std::vector<Fr> tmp_l_poly(n * l, 0);
    for (int i = 0; i < n * l; i++) tmp_l_poly[i] = a_L[i] - z;
    l_poly.push_back(tmp_l_poly);
    l_poly.push_back(s_L);
    std::vector<Fr> tmp_r_poly(n * l, 0);
    Fr cur = 1;
    for (int i = 0; i < n * l; i++) {
        tmp_r_poly[i] = a_R[i] + z;
        tmp_r_poly[i] *= cur;
        cur *= y;
    }
    
    Fr cur_z = z * z;
    for (int j = 0; j < l; j++) {
        cur = 1;
        for (int i = 0; i < n; i++) {
            tmp_r_poly[j * n + i] += cur_z * cur;
            cur *= 2;
        }
        cur_z *= z;
    }

    r_poly.push_back(tmp_r_poly);
    cur = 1;
    for (int i = 0; i < n * l; i++) {
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

    Fr delta;
    Fr y_nl, two_n;
    y_nl = Fr_pow(y, n * l);
    two_n = Fr_pow(2, n);
    delta = (y_nl - 1) / (y - 1) * (z - z * z);
    cur_z = z * z * z;
    for (int j = 0; j < l; j++) {
        delta -= (two_n - 1) * cur_z;
        cur_z *= z;
    }

        // Fr check_t = delta;
        // cur_z = z*z;
        // for (int j = 0; j < l; j++) {
        //     check_t += cur_z * v[j];
        //     cur_z *= z;
        // }
        // assert(t_0 == check_t);

    Fr tao_1, tao_2;
    tao_1.setRand();
    tao_2.setRand();
    G1 T_1, T_2;
    T_1 = g * t_1 + h * tao_1;
    T_2 = g * t_2 + h * tao_2;
    end = std::chrono::high_resolution_clock::now();
    prover_time += end - start;

    proof_size += 2 * G1::getSerializedByteSize(); //T_1, T_2
    //T_1, T_2 sent
    Fr x;
    x.setRand();

    start = std::chrono::high_resolution_clock::now();
    std::vector<Fr> l_vec = l_poly[0];
    for (int i = 0; i < n * l; i++) l_vec[i] += l_poly[1][i] * x;
    std::vector<Fr> r_vec = r_poly[0];
    for (int i = 0; i < n * l; i++) r_vec[i] += r_poly[1][i] * x;
    Fr t_hat = vec_product(l_vec, r_vec);
    Fr tao_x = tao_2 * x * x + tao_1 * x;
    cur_z = z * z;
    for (int j = 0; j < l; j++) {
        tao_x += cur_z * gamma[j];
        cur_z *= z;
    }
    Fr miu = alpha + rou * x;
    end = std::chrono::high_resolution_clock::now();
    prover_time += end - start;

    proof_size += sizeof(tao_x) + sizeof(miu) + sizeof(t_hat) + sizeof(l_vec) + sizeof(r_vec);

    start = std::chrono::high_resolution_clock::now();
    std::vector<G1> h_vec_prime(n * l);
    Fr cur_y = 1;
    for (int i = 0; i < n * l; i++) {
        h_vec_prime[i] = h_vec[i] * cur_y;
        cur_y *= (1 / y);
    }

    G1 lhs, rhs;
    lhs = g * t_hat + h * tao_x;
    rhs = g * delta + T_1 * x + T_2 * (x * x);
    cur_z = 1;
    for (int j = 0; j < l; j++) {
        rhs += V[j] * (z * z * cur_z);
        cur_z *= z;
    }
    assert(lhs == rhs);

    G1 P = A + S * x;
    cur_y = 1;
    for (int i = 0; i < n * l; i++) {
        P += g_vec[i] * (-z) + h_vec_prime[i] * (z * cur_y);
        cur_y *= y;
    }
    cur_z = z * z;
    for (int j = 0; j < l; j++) {
        cur = 1;
        for (int i = 0; i < n; i++) {
            P += h_vec_prime[j * n + i] * cur_z * cur;
            cur *= 2;
        }
        cur_z *= z;
    }

    rhs = h * miu;
    for (int i = 0; i < n * l; i++) {
        rhs += g_vec[i] * l_vec[i] + h_vec_prime[i] * r_vec[i];
    }
    assert(P == rhs);
    end = std::chrono::high_resolution_clock::now();
    verifier_time += end - start;

    seed.setRand();
    G1 u = g * seed;
    Fr xx;
    xx.setRand();
    G1 P_prime = P + h * (-miu) + u * (x * t_hat);
    assert(Inner_Product_Argument(g_vec, h_vec_prime, u * x, P_prime, l_vec, r_vec));

    std::cout << "1" << std::endl;

    std::cout << "Prover Time:" << prover_time.count() << std::endl;
    std::cout << "Verifier Time:" << verifier_time.count() << std::endl;
    std::cout << "Proof Size:" << ((double) proof_size) / 1024.0 << " KB" << std::endl;
    return 0;
}
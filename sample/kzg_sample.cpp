#include <iostream>
#include <vector>
#include <mcl/bls12_381.hpp>
#include <random>
#include "kzg.h"
#include "subgroups.h"

using namespace mcl::bls12;

int generateRandomNumber() {
    // 创建一个随机设备和一个随机数生成器
    std::random_device rd;               // 用于生成种子
    std::mt19937 gen(rd());              // 随机数生成器，使用Mersenne Twister算法
    std::uniform_int_distribution<> dis(0, 3); // 生成范围为0到3的均匀分布

    return dis(gen);
}

std::vector<Fr> generate_polynomial(const int &l) {
    std::vector<Fr> coeffs;
    coeffs.reserve(l + 1); 
    cybozu::RandomGenerator rg;
    Fr r;
    for (int i = 0; i <= l; i++) {
        r.setRand(rg);
        coeffs.push_back(r);
    }
    return coeffs;
}

std::vector<std::vector<Fr>> generate_polynomial_d2(const __int128 &l, const __int128 &k) {
    std::vector<std::vector<Fr>> coeffs;
    for (int i = 0; i <= l; i++) {
        std::vector<Fr> tmp;
        coeffs.push_back(tmp);
    }
    cybozu::RandomGenerator rg;
    Fr r;
    for (__int128 i = 0; i <= l; i++) {
        for (__int128 j = 0; j <= k; j++) {
            // r = generateRandomNumber();
            r.setRand(rg);
            coeffs[i].push_back(r);
        }
    }
    return coeffs;
}

int main() {
    initPairing(mcl::BLS12_381);
    G1 g1;
    hashAndMapToG1(g1, "ab", 2);
    G2 g2;
    hashAndMapToG2(g2, "abc", 3);

    int l = 8, k = 8; //degree of the polynomial
    std::vector<std::vector<Fr>> coeffs = generate_polynomial_d2(l, k);

    //kzg setup
    CRS crs = kzg_setup(g2);

    //P generate ploynomial commitment
    G1 commitment = kzg_commit_d2(coeffs, g1, crs.SK, crs.SK2);

    //V generate random challenging i
    Fr i, j;
    i.setRand();
    j.setRand();

    //P create witness
    std::pair<Fr, std::pair<G1, G1>> witness_pair = kzg_createWitness_d2(coeffs, g1, crs.SK, crs.SK2, i, j);

    //verify the eval
    bool flag = kzg_verifyEval_d2(commitment, g1, g2, crs.PK, crs.PK2, i, j, witness_pair.first, witness_pair.second);

    std::cout << flag << std::endl;
    return 0;
}
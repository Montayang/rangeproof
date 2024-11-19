#include "kzg.h"

Fr evaluate_polynomial_d2(const std::vector<std::vector<Fr>>& coeffs, const Fr& x, const Fr& y) {
    __int128 n = coeffs.size();       
    __int128 m = coeffs[0].size();
    Fr result = 0;
    for (__int128 i = 0; i < n; ++i) {
        Fr x_term = 1;
        for (__int128 xi = 0; xi < i; ++xi) {
            x_term = x_term * x;
        }
        for (__int128 j = 0; j < m; ++j) {
            Fr y_term = 1;
            for (__int128 yj = 0; yj < j; ++yj) {
                y_term = y_term * y;
            }
            result = (result + coeffs[i][j] * x_term * y_term);
        }
    }
    return result;
}

Fr evaluate_polynomial(const std::vector<Fr> &coeffs, const Fr &x) {
    Fr result = 0;
    Fr power = 1;
    for (const auto &coeff : coeffs) {
        result += coeff * power;
        power *= x;
    }
    return result;
}

std::pair<std::vector<Fr>, std::vector<Fr>> polynomial_division(const std::vector<Fr> &f, const std::vector<Fr> &d) {
    __int128 m = f.size() - 1;
    __int128 n = d.size() - 1;
    std::vector<Fr> quotient(f.size(), 0);
    std::vector<Fr> remainder = f;
    for (__int128 i = m; i >= n; --i) {
        Fr coefficient = remainder[i] / d[n];
        quotient[i - n] = coefficient;

        for (__int128 j = 0; j <= n; ++j) {
            remainder[i - j] -= coefficient * d[n - j];
        }
    }
    return make_pair(quotient, remainder);
}

std::pair<std::vector<std::vector<Fr>>, std::vector<std::vector<Fr>>> polynomial_division_d2(const std::vector<std::vector<Fr>> &f, const std::vector<Fr> &d1, const std::vector<Fr> &d2) {
    std::vector<std::vector<Fr>> quotient_x = f;
    std::vector<std::vector<Fr>> quotient_y = f;
    std::vector<std::vector<Fr>> tmp_remainder = f;
    for (__int128 i = 0; i < f.size(); ++i) {
        std::pair<std::vector<Fr>, std::vector<Fr>> tmp = polynomial_division(f[i], d2);
        quotient_y[i] = tmp.first;
        tmp_remainder[i] = tmp.second;
    }

    for (__int128 j = 0; j < tmp_remainder[0].size(); ++j) {
        std::vector<Fr> tmp_dividend;
        for (__int128 i = 0; i < tmp_remainder.size(); ++i) {
            tmp_dividend.push_back(tmp_remainder[i][j]);
        }
        std::pair<std::vector<Fr>, std::vector<Fr>> tmp = polynomial_division(tmp_dividend, d1);
        std::vector<Fr> tmpp = tmp.first;
        std::vector<Fr> tmppp = tmp.second;
        //for (auto i : tmppp) assert(i == 0);
        for (__int128 i = 0; i < tmp_remainder.size(); ++i) {
            if (i < tmpp.size()) quotient_x[i][j] = tmpp[i];
            else quotient_x[i][j] = 0;
        }  
    }
    // for (__int128 i = 0; i < 10; i++) {
    //     Fr ff1, ff2;
    //     ff1.setRand();
    //     ff2.setRand();
    //     assert(evaluate_polynomial_d2(f,ff1,ff2) == evaluate_polynomial_d2(quotient_x,ff1,ff2)*evaluate_polynomial(d1,ff1) +  evaluate_polynomial_d2(quotient_y,ff1,ff2)*evaluate_polynomial(d2,ff2));
    // }
    return make_pair(quotient_x, quotient_y);
}

// CRS kzg_setup(const G1 &g1, const G2 &g2, const __int128 &t1, const __int128 &t2) {
//     CRS crs;
//     crs.SK1.setRand();
//     crs.SK2.setRand();
//     crs.g1 = g1;
//     crs.g2 = g2;
//     std::vector<G1> tmp;
//     Fr cur = 1;
//     for (__int128 i = 0; i <= t1; i++) {
//         tmp.push_back(g1 * cur);
//         cur *= crs.SK1;
//     }
//     crs.PK1 = make_pair(tmp, g2 * crs.SK1);
//     std::vector<std::vector<G1>> tmp2;
//     cur = 1;
//     for (__int128 i = 0; i <= t1; i++) {
//         std::vector<G1> vec;
//         tmp2.push_back(vec);
//         Fr cur2 = 1;
//         for (__int128 j = 0; j <= t2; j++) {
//             tmp2[i].push_back(g1 * cur * cur2);
//             cur2 *= crs.SK2;
//         }
//         cur *= crs.SK1;
//     }
//     crs.PK2 = make_pair(tmp2, g2 * crs.SK2);
//     return crs;
// }

CRS kzg_setup(const G1 &g1, const G2 &g2) {
    CRS crs;
    crs.SK1.setRand();
    crs.SK2.setRand();
    crs.g1 = g1;
    crs.g2 = g2;
    crs.PK1 = g2 * crs.SK1;
    crs.PK2 = g2 * crs.SK2;
    return crs;
}

G1 kzg_commit(const std::vector<Fr> &coeffs, const CRS &crs) {
    G1 commitment;
    commitment = crs.g1 * evaluate_polynomial(coeffs, crs.SK1);
    return commitment;
}

G1 kzg_commit_d2(const std::vector<std::vector<Fr>> &coeffs, const CRS &crs) {
    G1 commitment;
    commitment = crs.g1 * evaluate_polynomial_d2(coeffs, crs.SK1, crs.SK2);
    return commitment;
}

// bool verify_kzg_open(const G1 &commitment, const std::vector<Fr> &coeffs, const CRS &crs) {
//     G1 eval;
//     eval -= eval; //set eval = 0
//     for (__int128 i = 0; i < coeffs.size(); ++i) {
//         eval += crs.PK1.first[i] * coeffs[i];
//     }
//     return commitment == eval;
// }

// bool verify_kzg_open_d2(const G1 &commitment, const std::vector<std::vector<Fr>> &coeffs, const CRS &crs) {
//     G1 eval;
//     eval -= eval; //set eval = 0
//     for (__int128 i = 0; i < coeffs.size(); ++i) {
//         for (__int128 j = 0; j < coeffs[i].size(); ++j) {
//             eval += crs.PK2.first[i][j] * coeffs[i][j];
//         }
//     }
//     return commitment == eval;
// }

std::pair<Fr, G1> kzg_createWitness(const std::vector<Fr> &coeffs, const CRS &crs, const Fr &i) {
    Fr Phi_i;
    G1 witness;
    Phi_i = evaluate_polynomial(coeffs, i);
    Fr Psi_alpha = (evaluate_polynomial(coeffs, crs.SK1) - Phi_i) / (crs.SK1 - i);
    witness = crs.g1 * Psi_alpha;
    return std::make_pair(Phi_i, witness);
}

std::pair<Fr, std::pair<G1, G1>> kzg_createWitness_d2(const std::vector<std::vector<Fr>> &coeffs, const CRS &crs, const Fr &i, const Fr &j) {
    Fr Phi_ij;
    G1 witness1, witness2;
    Phi_ij = evaluate_polynomial_d2(coeffs, i, j);
    std::vector<std::vector<Fr>> dividend = coeffs;
    dividend[0][0] -= Phi_ij;
    std::vector<Fr> d1 = {-i, 1}, d2 = {-j, 1};
    std::pair<std::vector<std::vector<Fr>>, std::vector<std::vector<Fr>>> Psi = polynomial_division_d2(dividend, d1, d2);
    // std::cout << "i:" << i << std::endl;
    // std::cout << "j:" << j << std::endl;
    // std::cout << "Phi:" << Phi_ij << std::endl;
    // std::cout << "coeffs:" <<std::endl;
    // for (auto i : coeffs) {
    //     for (auto j : i) std::cout << j <<"  ";
    //     std::cout << std::endl;
    // }
    // std::cout << "Psi1:" <<std::endl;
    // for (auto i : Psi.first) {
    //     for (auto j : i) std::cout << j <<"  ";
    //     std::cout << std::endl;
    // }
    // std::cout << "Psi2:" <<std::endl;
    // for (auto i : Psi.second) {
    //     for (auto j : i) std::cout << j <<"  ";
    //     std::cout << std::endl;
    // }
    witness1 = crs.g1 * evaluate_polynomial_d2(Psi.first, crs.SK1, crs.SK2);
    witness2 = crs.g1 * evaluate_polynomial_d2(Psi.second, crs.SK1, crs.SK2);
    return std::make_pair(Phi_ij, std::make_pair(witness1, witness2));
}

bool kzg_verifyEval(const G1 &commitment, const CRS &crs, const Fr &i, const Fr &Phi_i, const G1 &witness) {
    GT lhs, rhs1, rhs2;
    pairing(lhs, commitment, crs.g2);
    pairing(rhs1, witness, crs.PK1 - crs.g2 * i);
    pairing(rhs2, crs.g1, crs.g2);
    GT::pow(rhs2, rhs2, Phi_i);
    return lhs == rhs1 * rhs2;
}

bool kzg_verifyEval_d2(const G1 &commitment, const CRS &crs, const Fr &i, const Fr &j, const Fr &Phi_ij, const std::pair<G1, G1> &witness) {
    GT lhs, rhs1, rhs2, rhs3;
    pairing(lhs, commitment, crs.g2);
    pairing(rhs1, witness.first, crs.PK1 - crs.g2 * i);
    pairing(rhs2, witness.second, crs.PK2 - crs.g2 * j);
    pairing(rhs3, crs.g1, crs.g2);
    GT::pow(rhs3, rhs3, Phi_ij);
    return lhs == rhs1 * rhs2 * rhs3;
}
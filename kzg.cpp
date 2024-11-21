#include "kzg.h"
#include <unordered_map>

// G1 pippenger(const std::vector<Fr>& scalars, const std::vector<G1>& points, int window_size = 4) {
//     size_t n = scalars.size();

//     // 计算窗口的最大值和窗口数量
//     int max_window_value = (1 << window_size) - 1; // 最大窗口值 2^w - 1
//     int num_windows = (256 + window_size - 1) / window_size; // 标量为256位，计算窗口数（向上取整）

//     // 初始化结果点
//     G1 result;
//     result -= result;

//     // 遍历每个窗口
//     for (int window_index = 0; window_index < num_windows; ++window_index) {
//         // 预计算每个窗口值的点和
//         std::unordered_map<int, G1> precomputed_sums;
//         for (int i = 1; i <= max_window_value; ++i) {
//             precomputed_sums[i] -= precomputed_sums[i]; // 初始化每个窗口值的点和为零点
//         }

//         // 遍历所有标量和点
//         for (size_t i = 0; i < n; ++i) {
//             // 提取当前窗口的值
//             int shift = window_index * window_size;
//             Fr mask = Fr(max_window_value);
//             Fr window_value = (scalars[i] >> shift) & mask;

//             // 将点添加到对应的窗口值
//             if (window_value != 0) { // 忽略零窗口值
//                 precomputed_sums[window_value] = precomputed_sums[window_value] + points[i];
//             }
//         }

//         // 从最大到最小遍历窗口值，累加部分和
//         G1 window_sum;
//         window_sum -= window_sum;
//         for (int k = max_window_value; k > 0; --k) {
//             window_sum = window_sum + precomputed_sums[k];
//         }

//         // 将窗口和按位移权重累加到结果中
//         for (int j = 0; j < window_size; ++j) {
//             result = result + result; // 等价于乘以 2
//         }
//         result = result + window_sum;
//     }

//     return result;
// }

Fr evaluate_polynomial_d2(const std::vector<std::vector<Fr>>& coeffs, const Fr& x, const Fr& y) {
    int n = coeffs.size();       
    int m = coeffs[0].size();
    Fr result = 0;
    for (int i = 0; i < n; ++i) {
        Fr x_term = 1;
        for (int xi = 0; xi < i; ++xi) {
            x_term = x_term * x;
        }
        for (int j = 0; j < m; ++j) {
            Fr y_term = 1;
            for (int yj = 0; yj < j; ++yj) {
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
    int m = f.size() - 1;
    int n = d.size() - 1;
    std::vector<Fr> quotient(f.size(), 0);
    std::vector<Fr> remainder = f;
    for (int i = m; i >= n; --i) {
        Fr coefficient = remainder[i] / d[n];
        quotient[i - n] = coefficient;

        for (int j = 0; j <= n; ++j) {
            remainder[i - j] -= coefficient * d[n - j];
        }
    }
    return make_pair(quotient, remainder);
}

std::pair<std::vector<std::vector<Fr>>, std::vector<std::vector<Fr>>> polynomial_division_d2(const std::vector<std::vector<Fr>> &f, const std::vector<Fr> &d1, const std::vector<Fr> &d2) {
    std::vector<std::vector<Fr>> quotient_x = f;
    std::vector<std::vector<Fr>> quotient_y = f;
    std::vector<std::vector<Fr>> tmp_remainder = f;
    for (int i = 0; i < f.size(); ++i) {
        std::pair<std::vector<Fr>, std::vector<Fr>> tmp = polynomial_division(f[i], d2);
        quotient_y[i] = tmp.first;
        tmp_remainder[i] = tmp.second;
    }

    for (int j = 0; j < tmp_remainder[0].size(); ++j) {
        std::vector<Fr> tmp_dividend;
        for (int i = 0; i < tmp_remainder.size(); ++i) {
            tmp_dividend.push_back(tmp_remainder[i][j]);
        }
        std::pair<std::vector<Fr>, std::vector<Fr>> tmp = polynomial_division(tmp_dividend, d1);
        std::vector<Fr> tmpp = tmp.first;
        std::vector<Fr> tmppp = tmp.second;
        //for (auto i : tmppp) assert(i == 0);
        for (int i = 0; i < tmp_remainder.size(); ++i) {
            if (i < tmpp.size()) quotient_x[i][j] = tmpp[i];
            else quotient_x[i][j] = 0;
        }  
    }
    // for (int i = 0; i < 10; i++) {
    //     Fr ff1, ff2;
    //     ff1.setRand();
    //     ff2.setRand();
    //     assert(evaluate_polynomial_d2(f,ff1,ff2) == evaluate_polynomial_d2(quotient_x,ff1,ff2)*evaluate_polynomial(d1,ff1) +  evaluate_polynomial_d2(quotient_y,ff1,ff2)*evaluate_polynomial(d2,ff2));
    // }
    return make_pair(quotient_x, quotient_y);
}

CRS kzg_setup(const G1 &g1, const G2 &g2, const int &t1, const int &t2) {
    CRS crs;
    Fr SK1, SK2;
    SK1.setRand();
    SK2.setRand();
    crs.g1 = g1;
    crs.g2 = g2;
    std::vector<G1> tmp;
    Fr cur = 1;
    for (int i = 0; i <= t1; i++) {
        tmp.push_back(g1 * cur);
        cur *= SK1;
    }
    crs.PK1 = make_pair(tmp, g2 * SK1);
    std::vector<std::vector<G1>> tmp2;
    cur = 1;
    for (int i = 0; i <= t1; i++) {
        std::vector<G1> vec;
        tmp2.push_back(vec);
        Fr cur2 = 1;
        for (int j = 0; j <= t2; j++) {
            tmp2[i].push_back(g1 * cur * cur2);
            cur2 *= SK2;
        }
        cur *= SK1;
    }
    crs.PK2 = make_pair(tmp2, g2 * SK2);
    return crs;
}

// CRS kzg_setup(const G1 &g1, const G2 &g2) {
//     CRS crs;
//     crs.SK1.setRand();
//     crs.SK2.setRand();
//     crs.g1 = g1;
//     crs.g2 = g2;
//     crs.PK1 = g2 * crs.SK1;
//     crs.PK2 = g2 * crs.SK2;
//     return crs;
// }

G1 kzg_commit(const std::vector<Fr> &coeffs, const CRS &crs) {
    G1 eval;
    int n = coeffs.size();
    G1* PK1Array = new G1[n];
    Fr* coeffsArray = new Fr[n];
    for (int i = 0; i < n; i++) {
        PK1Array[i] = crs.PK1.first[i];
        coeffsArray[i] = coeffs[i];
    }
    G1::mulVec(eval, PK1Array, coeffsArray, n);

    delete[] PK1Array;
    delete[] coeffsArray;
    return eval;
}

G1 kzg_commit_d2(const std::vector<std::vector<Fr>> &coeffs, const CRS &crs) {
    G1 eval;
    int n = 0;
    for (int i = 0; i < coeffs.size(); ++i) {
        n += coeffs[i].size();
    }
    G1* PK2Array = new G1[n];
    Fr* coeffsArray = new Fr[n];
    int cnt = 0;
    for (int i = 0; i < coeffs.size(); ++i) {
        for (int j = 0; j < coeffs[i].size(); ++j) {
            PK2Array[cnt] = crs.PK2.first[i][j];
            coeffsArray[cnt] = coeffs[i][j];
            cnt++;
        }
    }
    G1::mulVec(eval, PK2Array, coeffsArray, n);

    delete[] PK2Array;
    delete[] coeffsArray;
    return eval;
}

// bool verify_kzg_open(const G1 &commitment, const std::vector<Fr> &coeffs, const CRS &crs) {
//     G1 eval;
//     eval -= eval; //set eval = 0
//     for (int i = 0; i < coeffs.size(); ++i) {
//         eval += crs.PK1.first[i] * coeffs[i];
//     }
//     return commitment == eval;
// }

// bool verify_kzg_open_d2(const G1 &commitment, const std::vector<std::vector<Fr>> &coeffs, const CRS &crs) {
//     G1 eval;
//     eval -= eval; //set eval = 0
//     for (int i = 0; i < coeffs.size(); ++i) {
//         for (int j = 0; j < coeffs[i].size(); ++j) {
//             eval += crs.PK2.first[i][j] * coeffs[i][j];
//         }
//     }
//     return commitment == eval;
// }

std::pair<Fr, G1> kzg_createWitness(const std::vector<Fr> &coeffs, const CRS &crs, const Fr &i) {
    Fr Phi_i;
    Phi_i = evaluate_polynomial(coeffs, i);
    std::vector<Fr> Psi = coeffs;
    Psi[0] -= Phi_i;
    Psi = polynomial_division(Psi, {-i, 1}).first;
    G1 witness;
    int n = Psi.size();
    G1* PK1Array = new G1[n];
    Fr* PsiArray = new Fr[n];
    for (int t = 0; t < n; t++) {
        PK1Array[t] = crs.PK1.first[t];
        PsiArray[t] = Psi[t];
    }
    G1::mulVec(witness, PK1Array, PsiArray, n);

    delete[] PK1Array;
    delete[] PsiArray;
    return std::make_pair(Phi_i, witness);
}

std::pair<Fr, std::pair<G1, G1>> kzg_createWitness_d2(const std::vector<std::vector<Fr>> &coeffs, const CRS &crs, const Fr &i, const Fr &j) {
    Fr Phi_ij;
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
    G1 witness1, witness2;
    int n1 = 0;
    for (int t = 0; t < Psi.first.size(); ++t) {
        n1 += Psi.first[t].size();
    }
    G1* PK2Array = new G1[n1];
    Fr* PsiArray = new Fr[n1];
    int cnt = 0;
    for (int t = 0; t < Psi.first.size(); ++t) {
        for (int k = 0; k < Psi.first[t].size(); ++k) {
            PK2Array[cnt] = crs.PK2.first[t][k];
            PsiArray[cnt] = Psi.first[t][k];
            cnt++;
        }
    }
    G1::mulVec(witness1, PK2Array, PsiArray, n1);

    delete[] PK2Array;
    delete[] PsiArray;
    
    int n2 = 0;
    for (int t = 0; t < Psi.second.size(); ++t) {
        n2 += Psi.second[t].size();
    }
    G1* PK2Array_ = new G1[n2];
    Fr* PsiArray_ = new Fr[n2];
    cnt = 0;
    for (int t = 0; t < Psi.second.size(); ++t) {
        for (int k = 0; k < Psi.second[t].size(); ++k) {
            PK2Array_[cnt] = crs.PK2.first[t][k];
            PsiArray_[cnt] = Psi.second[t][k];
            cnt++;
        }
    }
    G1::mulVec(witness2, PK2Array_, PsiArray_, n2);

    delete[] PK2Array_;
    delete[] PsiArray_;
    
    return std::make_pair(Phi_ij, std::make_pair(witness1, witness2));
}

bool kzg_verifyEval(const G1 &commitment, const CRS &crs, const Fr &i, const Fr &Phi_i, const G1 &witness) {
    GT lhs, rhs1, rhs2;
    pairing(lhs, commitment, crs.g2);
    pairing(rhs1, witness, crs.PK1.second - crs.g2 * i);
    pairing(rhs2, crs.g1, crs.g2);
    GT::pow(rhs2, rhs2, Phi_i);
    return lhs == rhs1 * rhs2;
}

bool kzg_verifyEval_d2(const G1 &commitment, const CRS &crs, const Fr &i, const Fr &j, const Fr &Phi_ij, const std::pair<G1, G1> &witness) {
    GT lhs, rhs1, rhs2, rhs3;
    pairing(lhs, commitment, crs.g2);
    pairing(rhs1, witness.first, crs.PK1.second - crs.g2 * i);
    pairing(rhs2, witness.second, crs.PK2.second - crs.g2 * j);
    pairing(rhs3, crs.g1, crs.g2);
    GT::pow(rhs3, rhs3, Phi_ij);
    return lhs == rhs1 * rhs2 * rhs3;
}
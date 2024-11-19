#include <iostream>
#include <vector>
#include <mcl/bls12_381.hpp>
#include <random>
#include <cstring>
#include <cmath>
#include "kzg.h"
#include "subgroups.h"

const double PI = acos(-1);

using std::vector;

// 快速傅里叶变换（FFT）函数
void fft(vector<complex<double>>& a, bool invert) {
    int n = a.size();
    
    // 重新排列数据
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        if (i < j) {
            std::swap(a[i], a[j]);
        }
    }

    // 计算 FFT
    for (int len = 2; len <= n; len <<= 1) {
        double angle = 2 * PI / len * (invert ? -1 : 1);
        complex<double> wlen(cos(angle), sin(angle));
        for (int i = 0; i < n; i += len) {
            complex<double> w(1);
            for (int j = 0; j < len / 2; j++) {
                complex<double> u = a[i + j];
                complex<double> v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }

    // 反转 FFT 时，将所有值除以 n
    if (invert) {
        for (complex<double>& x : a) {
            x /= n;
        }
    }
}

// 使用 FFT 计算两个多项式乘积的系数
vector<Fr> multiply_polynomials(const vector<Fr>& poly1, const vector<Fr>& poly2) {
    int n = 1;
    while (n < poly1.size() + poly2.size()) {
        n <<= 1;
    }

    // 扩展多项式，使其长度为 n 的最小 2 的幂
    vector<Fr> f1(poly1.begin(), poly1.end());
    vector<Fr> f2(poly2.begin(), poly2.end());
    f1.resize(n);
    f2.resize(n);

    // 执行 FFT
    fft(f1, false);
    fft(f2, false);

    // 点乘结果
    for (int i = 0; i < n; i++) {
        f1[i] *= f2[i];
    }

    // 执行逆 FFT
    fft(f1, true);

    // 将复数部分舍弃，返回实数部分
    vector<double> result(n);
    for (int i = 0; i < n; i++) {
        result[i] = round(f1[i].real());  // 取实部，并四舍五入
    }

    return result;
}

int main() {
    // 测试用例
    vector<double> poly1 = {1, 2, 3};    // 多项式 A(x) = 1 + 2x + 3x^2
    vector<double> poly2 = {4, 5, 6};    // 多项式 B(x) = 4 + 5x + 6x^2

    vector<double> result = multiply_polynomials(poly1, poly2);

    // 输出乘积多项式的系数
    std::cout << "Product coefficients: ";
    for (double coeff : result) {
        std::cout << coeff << " ";
    }
    std::cout << std::endl;

    return 0;
}

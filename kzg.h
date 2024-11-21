#ifndef KZG_H
#define KZG_H

#include <iostream>
#include <vector>
#include <mcl/bls12_381.hpp>
#include <random>

using namespace mcl::bls12;

Fr evaluate_polynomial_d2(const std::vector<std::vector<Fr>>& coeffs, const Fr& x, const Fr& y);

Fr evaluate_polynomial(const std::vector<Fr> &coeffs, const Fr &x);

std::pair<std::vector<Fr>, std::vector<Fr>> polynomial_division(const std::vector<Fr> &f, const std::vector<Fr> &d);

//only for exact division
std::pair<std::vector<std::vector<Fr>>, std::vector<std::vector<Fr>>> polynomial_division_d2(const std::vector<std::vector<Fr>> &f, const std::vector<Fr> &d1, const std::vector<Fr> &d2);

struct CRS {
    G1 g1;
    G2 g2;
    std::pair<std::vector<G1>, G2> PK1; //first is i * g1, second is sk * g2
    std::pair<std::vector<std::vector<G1>>, G2> PK2; //for kzg_d2
};

CRS kzg_setup(const G1 &g1, const G2 &g2, const int &t1, const int &t2);

G1 kzg_commit(const std::vector<Fr> &coeffs, const CRS &crs);

G1 kzg_commit_d2(const std::vector<std::vector<Fr>> &coeffs, const CRS &crs);

//bool verify_kzg_open(const G1 &commitment, const std::vector<Fr> &coeffs, const CRS &crs);

//bool verify_kzg_open_d2(const G1 &commitment, const std::vector<std::vector<Fr>> &coeffs, const CRS &crs);

std::pair<Fr, G1> kzg_createWitness(const std::vector<Fr> &coeffs, const CRS &crs, const Fr &i);

std::pair<Fr, std::pair<G1, G1>> kzg_createWitness_d2(const std::vector<std::vector<Fr>> &coeffs, const CRS &crs ,const Fr &i, const Fr &j);

bool kzg_verifyEval(const G1 &commitment, const CRS &crs, const Fr &i, const Fr &Phi_i, const G1 &witness);

bool kzg_verifyEval_d2(const G1 &commitment, const CRS &crs, const Fr &i, const Fr &j, const Fr &Phi_ij, const std::pair<G1, G1> &witness);

#endif
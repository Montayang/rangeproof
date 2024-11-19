#ifndef SUBGROUPS_H
#define SUBGROUPS_H

#include <iostream>
#include <vector>
#include <mcl/bls12_381.hpp>

using namespace mcl::bls12;

class subgroup {
public:
    __int128 order;
    Fr generator;
    //std::vector<Fr> elements;

public:
    subgroup(__int128 size);

    const __int128& getOrder() const;

    const Fr& getGenerator() const;

    const Fr& operator[](int index) const;
};

#endif
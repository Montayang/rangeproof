#ifndef SUBGROUPS_H
#define SUBGROUPS_H

#include <iostream>
#include <vector>
#include <mcl/bn.hpp>

using namespace mcl::bn;

class subgroup {
public:
    Fr generator;
    //std::vector<Fr> elements;

public:
    subgroup(int size);

    // const int& getOrder() const;

    const Fr& getGenerator() const;

    //const Fr& operator[](int index) const;
};

#endif
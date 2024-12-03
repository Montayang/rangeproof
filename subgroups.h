<<<<<<< HEAD
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

=======
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

>>>>>>> ea97283b7d0d33a96163f269cae4f9e92a8490e0
#endif
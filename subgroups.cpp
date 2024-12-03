#include "subgroups.h"

subgroup::subgroup(int size) {
    if (size == 0) return;
    if (size == 1) {
        generator = 1;
        //elements.push_back(1);
        return;
    }
    size = __builtin_ctz(size);
    generator = -1;
    while (--size) {
        Fr::squareRoot(generator, generator);
    }
    // Fr minus_1 = -1, root = 3;
    // Fr tmp = size;
    // Fr::div(tmp, minus_1, tmp);
    // Fr::pow(generator, root, tmp);
}

// const int& subgroup::getOrder() const { 
//     return order;
// }

const Fr& subgroup::getGenerator() const { 
    return generator;
}

// const Fr& subgroup::operator[](int index) const {
//     return elements.at(index);
// }
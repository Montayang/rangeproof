#include "subgroups.h"

subgroup::subgroup(__int128 size) {
    order = size;
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
    // Fr data;
    // data = 1;
    // elements.push_back(data);
    // for (__int128 i = 1; i < order; i++) {
    //     data = data * generator;
    //     elements.push_back(data);
    // }
}

const __int128& subgroup::getOrder() const { 
    return order;
}

const Fr& subgroup::getGenerator() const { 
    return generator;
}

// const Fr& subgroup::operator[](int index) const {
//     return elements.at(index);
// }
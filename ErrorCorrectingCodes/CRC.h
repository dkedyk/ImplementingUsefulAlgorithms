#ifndef IGMDK_CRC_H
#define IGMDK_CRC_H
#include <limits>

namespace igmdk{

using namespace std;

class CRC32
{
    uint32_t polynomial, constant[256];
public:
    CRC32(uint32_t thePolynomial = 0xFA567D89u):polynomial(thePolynomial)
    {
        for(int i = 0; i < 256; ++i)
        {
            constant[i] = i << 24;//make extended c
            for(int j = 0; j < 8; ++j) constant[i] =
                (constant[i] << 1) ^ (constant[i] >> 31 ? polynomial : 0);
        }
    }
    uint32_t hash(unsigned char* array, int size, uint32_t crc = 0)
    {
        assert(numeric_limits<unsigned char>::digits == 8);
        for(int i = 0; i < size; ++i)
            crc = (crc << 8) ^ constant[(crc >> 24) ^ array[i]];
        return crc;
    }
};

}//end namespace
#endif

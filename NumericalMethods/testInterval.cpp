#include <cassert>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <limits>
#include <memory> //for shared ptr
#include "IntervalNumber.h"
using namespace std;
using namespace igmdk;

int main()
{
    IntervalNumber<> iv1(2.0), iv2(3.14);
    DEBUG("+");
    (iv1 + iv2).debug();
    DEBUG("-");
    (iv1 - iv2).debug();
    DEBUG("*");
    (iv1 * iv2).debug();
    DEBUG("/");
    (iv1 / iv2).debug();
    DEBUG(iv1.isfinite());
    DEBUG(iv1.isnan());

    return 0;
}

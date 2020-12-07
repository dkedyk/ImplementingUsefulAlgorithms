#include "Sobol.h"
using namespace igmdk;

void testSobol()
{
    DEBUG(Sobol::maxD());
    Sobol so(2);
    //estimate Pi
    double n = 0, nTotal = pow(10, 8);
    for(int i = 0; i < nTotal; ++i)
    {
        double x = so.getU01Value(0), y = so.getU01Value(1);
        if(x * x + y * y <= 1) ++n;
        so.next();
    }
    double pi = 4 * n/nTotal;
    DEBUG(pi);
    double error = 4 * 2 * PI() * log(nTotal) * log(nTotal)/nTotal;
    DEBUG(error);
}

int main(int argc, char *argv[])
{
    testSobol();
    return 0;
}

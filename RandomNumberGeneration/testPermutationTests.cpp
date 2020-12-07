#include "PermutationTests.h"
#include "Statistics.h"
#include "../NumericalMethods/NumericalMethods.h"
#include "../NumericalMethods/Matrix.h"
#include "../Utils/Debug.h"
#include "../ExternalMemoryAlgorithms/CSV.h"
using namespace igmdk;

struct meaner
{
    double operator()(Vector<double> const& observations)const
    {
        double result = observations[0];
        for(int i = 1; i < observations.getSize(); ++i)
        {
            result += observations[i];
        }
        return result / observations.getSize();
    }
};
//duplicated
double permutationPairedTest(Vector<pair<double, double> > const& data,
    int b = 10000)
{
    PairedTestLocationPermuter<meaner> p = {data};
    return permutationTest(p, b);
}

void testPermPair()
{
    Vector<pair<double, double> > a;
    a.append(pair<double, double>(51.2, 45.8));
    a.append(pair<double, double>(46.5, 41.3));
    a.append(pair<double, double>(24.1, 15.8));
    a.append(pair<double, double>(10.2, 11.1));
    a.append(pair<double, double>(65.3, 58.5));
    a.append(pair<double, double>(92.1, 70.3));
    a.append(pair<double, double>(30.3, 31.6));
    a.append(pair<double, double>(49.2, 35.4));
    DEBUG(permutationPairedTest(a, 10000000));//seems to converge to 0.0236, t-test gives 0.0283
}

template<typename LOCATION_F> pair<double, double> permutationPairedConf(
    Vector<pair<double, double> > const& data, double a = 0.05, int b = 10000)
{
    PairedTestLocationPermuter<LOCATION_F> p = {data};
    return permutationConf(p, a, b);
}

void testPermConf()
{
    Vector<pair<double, double> > a;
    a.append(pair<double, double>(51.2, 45.8));
    a.append(pair<double, double>(46.5, 41.3));
    a.append(pair<double, double>(24.1, 15.8));
    a.append(pair<double, double>(10.2, 11.1));
    a.append(pair<double, double>(65.3, 58.5));
    a.append(pair<double, double>(92.1, 70.3));
    a.append(pair<double, double>(30.3, 31.6));
    a.append(pair<double, double>(49.2, 35.4));
    pair<double, double> conf = permutationPairedConf<meaner>(a);
    DEBUG(conf.first);
    DEBUG(conf.second);
    //DEBUG(permutationPairedTest(a, 10000000));//seems to converge to 0.0236
}

int main(int argc, char *argv[])
{
    testPermPair();
    return 0;
    testPermConf();
    return 0;

    return 0;
}

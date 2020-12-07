#ifndef IGMDK_FRIEDMAN_TEST_H
#define IGMDK_FRIEDMAN_TEST_H
#include "../RandomNumberGeneration/Statistics.h"
#include "../RandomNumberGeneration/DistributionTests.h"
#include "../Utils/Vector.h"
namespace igmdk{

pair<double,Vector<double> > FriedmanPValue(Vector<Vector<double> > const& a)
{//a[i] is vector of responses on domain i
    assert(a.getSize() > 0 && a[0].getSize() > 1);
    int n = a.getSize(), k = a[0].getSize();
    double aveRank = (k + 1)/2.0, SSAlternative = 0, SSTotal = 0;
    Vector<double> alternativeRankSums(k);
    for(int i = 0; i < n; ++i)
    {
        assert(a[i].getSize() == k);
        Vector<double> ri = convertToRanks(a[i]);
        for(int j = 0; j < k; ++j)
        {
            alternativeRankSums[j] += ri[j];
            SSTotal += (ri[j] - aveRank) * (ri[j] - aveRank);
        }
    }
    for(int j = 0; j < k; ++j)
    {
        double temp = alternativeRankSums[j] - n * aveRank;
        SSAlternative += temp * temp;
    }
    double p =
        1 - evaluateChiSquaredCdf(SSAlternative * (k - 1)/SSTotal, k - 1);
    return make_pair(p, alternativeRankSums);
}

}
#endif

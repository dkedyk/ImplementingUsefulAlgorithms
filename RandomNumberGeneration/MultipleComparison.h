#ifndef IGMDK_MULTIPLE_COMPARISON_H
#define IGMDK_MULTIPLE_COMPARISON_H
#include "Statistics.h"
#include "../Sorting/Sort.h"
#include "../NumericalMethods/Matrix.h"
#include <cmath>
namespace igmdk{

double find2SidedConfZBonf(int k, double conf = 0.95)
{
    assert(k > 0);
    return find2SidedConfZ(1 - (1 - conf)/k);
}

double findNemenyiSignificantAveRankDiff(int k, int n, bool forControl = false,
    double conf = 0.95)
{//invert rank sum formula
    int nPairs = k * (k + 1)/2;
    double q = sqrt(nPairs/(3.0 * n)),
        z = find2SidedConfZBonf(forControl ? k : nPairs, conf);
    return q * z/n;//for rank average, not sum
}

void HolmAdjust(Vector<double>& pValues)
{
    int k = pValues.getSize();
    Vector<int> indices(k);
    for(int i = 0; i < k; ++i) indices[i] = i;
    IndexComparator<double> c(pValues.getArray());
    quickSort(indices.getArray(), 0, k - 1, c);
    for(int i = 0; i < k; ++i) pValues[indices[i]] =  min(1.0, max(i > 0 ?
        pValues[indices[i - 1]] : 0, (k - i) * pValues[indices[i]]));
}

void FDRAdjust(Vector<double>& pValues)
{
    int k = pValues.getSize();
    Vector<int> indices(k);
    for(int i = 0; i < k; ++i) indices[i] = i;
    IndexComparator<double> c(pValues.getArray());
    quickSort(indices.getArray(), 0, k - 1, c);
    for(int i = k - 1; i >= 0; --i) pValues[indices[i]] = min(i < k - 1 ?
        pValues[indices[i + 1]] : 1, pValues[indices[i]] * k/(i + 1));
}

Vector<double> FriedmanRankSums(Vector<Vector<double> > const& a)
{//a[i] is vector of responses on domain i
    assert(a.getSize() > 0 && a[0].getSize() > 1);
    int n = a.getSize(), k = a[0].getSize();
    Vector<double> alternativeRankSums(k);
    for(int i = 0; i < n; ++i)
    {
        assert(a[i].getSize() == k);
        Vector<double> ri = convertToRanks(a[i]);
        for(int j = 0; j < k; ++j) alternativeRankSums[j] += ri[j];
    }
    return alternativeRankSums;
}

double NemenyiAllPairsPValueUnadjusted(double r1, double r2, int n, int k)
    {return 1 - approxNormalCDF(abs(r1 - r2)/sqrt(n * k * (k + 1)/6.0));}
Matrix<double> RankTestAllPairs(Vector<Vector<double> > const& a,
    double NemenyiALevel = 0.05, bool useFDR = false)
{
    Vector<double> rankSums = FriedmanRankSums(a);
    int n = a.getSize(), k = rankSums.getSize();
    Vector<double> temp(k * (k - 1)/2);
    for(int i = 1, index = 0; i < k; ++i) for(int j = 0; j < i; ++j)
        temp[index++] = NemenyiAllPairsPValueUnadjusted(rankSums[i],
            rankSums[j], n, k);
    if(useFDR) FDRAdjust(temp);
    else HolmAdjust(temp);
    Matrix<double> result = Matrix<double>::identity(k);
    for(int i = 1, index = 0; i < k; ++i) for(int j = 0; j < i; ++j)
        result(i, j) = result(j, i) = temp[index++];
    return result;
}

}//end namespace
#endif

#include "SensitivityAnalysis.h"
#include "Statistics.h"
using namespace igmdk;

struct SobolIndexFunctor
{
    double operator()(Vector<double> const& x)const
    {
        double sum = 0;
        for(int j = 0; j < x.getSize(); ++j) sum += x[j];
        return sum;
    }
};
void testSobolIndex()
{
    int n = 10000, D = 5;
    Vector<pair<Vector<double>, double> > data(n);
    SobolIndexFunctor f;
    for(int i = 0; i < n; ++i)
    {
        Vector<double> x(D);
        double y = 0;
        for(int j = 0; j < D; ++j) x[j] = GlobalRNG().normal01() * j;
        data[i] = make_pair(x, f(x));
    }
    pair<Vector<double>, Vector<pair<double, double>>> result = findSobolIndicesSaltelli(data, f);
    Vector<double> indices = result.first;
    for(int j = 0; j < D; ++j) DEBUG(indices[j]);
    for(int j = 0; j < D; ++j)
    {
        DEBUG(result.second[j].first);
        DEBUG(result.second[j].second);
    }

    Vector<double> correct(D);
    for(int j = 0; j < D; ++j) correct[j] = 1.0 * j * j;
    normalizeProbs(correct);
    double mse = 0;
    for(int j = 0; j < D; ++j) mse += pow(correct[j] - indices[j], 2);
    DEBUG(sqrt(mse));
}

int main(int argc, char *argv[])
{
    testSobolIndex();
    return 0;
}

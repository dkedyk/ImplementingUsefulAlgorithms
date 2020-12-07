#include "Correlation.h"
#include "Bootstrap.h"
#include "../NumericalMethods/NumericalMethods.h"
#include "../NumericalMethods/Matrix.h"
#include "../Utils/Debug.h"
#include "../ExternalMemoryAlgorithms/CSV.h"
#include "testCommon.h"
using namespace igmdk;

template<typename SAMPLER>
void testPearsonConfs(Vector<Vector<string> >& matrix)
{
    int n = 30;
    typedef pair<double, double> X;
    typedef correr F;
    F f;
    pair<double, pair<double, double> > target = testBootstrapHelper<SAMPLER, X>(f, n);
    DEBUG("Exact Simulation");
    DEBUG(target.first);
    DEBUG(target.second.first);
    DEBUG(target.second.second);
    SAMPLER s;
    BootstrapSimulationResult p("Normal"),
        bt("Bootstrap-t-capped");
    for(int i = 0; i < 10000; ++i)
    {
        Vector<X> samples(n);
        for(int j = 0; j < n; ++j) samples[j] = s();
        double stat = f(samples);
        double left = target.first + stat - target.second.second;
        double right = target.first + stat - target.second.first;
        pair<double, double> exact(left, right);
        BasicBooter<F, X> booter(samples);
        p.addValue(PearsonCorrelationConf(PearsonCorrelation(samples), samples.getSize()), target.first, exact, stat);
        //p.addValue(PearsonCorrelationConf(PearsonCorrelation(samples), samples.getSize()), target.first, exact, stat);
        booter.resample = booter.data;//reset
        bt.addValue(bootstrapTIntervalCapped(booter), target.first, exact, stat);
    }
    p.print(matrix);
    //bc.print(matrix);
    //bmr.print(matrix);
    //bm.print(matrix);
    bt.print(matrix);
    //br.print(matrix);
}

void testPearsonConfsDriver()
{
    Vector<Vector<string> > matrix;
    DEBUG("Normal Error Spearman");
    matrix.append(Vector<string>());
    matrix.lastItem().append("Normal Error Almost Line");
    testPearsonConfs<NErrorPairSampler>(matrix);
    matrix.append(Vector<string>());
    matrix.lastItem().append("Normal Error Cubic");
    testPearsonConfs<CubicPairSampler>(matrix);

    makeConfsReport("Pearson", matrix);
}

int main(int argc, char *argv[])
{
    testPearsonConfsDriver();
    return 0;
}

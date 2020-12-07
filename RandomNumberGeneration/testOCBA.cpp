#include "OCBA.h"
#include "Statistics.h"
#include "../NumericalMethods/NumericalMethods.h"
#include "../NumericalMethods/Matrix.h"
#include "../Utils/Debug.h"
#include "../ExternalMemoryAlgorithms/CSV.h"
using namespace igmdk;

void testNormalSummary()
{

	Vector<NormalSummary> hha;
	NormalSummary n1(10, 100.0/8);
	NormalSummary n2(20, 81.0/8);
	NormalSummary n3(22, 144.0/8);
	hha.append(n1);
	hha.append(n2);
	hha.append(n3);
	DEBUG(isNormal0BestBonf(hha, 0.05));//wont match Chen book - k adjustment

	Vector<NormalSummary> hha2;
	NormalSummary n4(1, 0.1);
	NormalSummary n5(2, 0.2);
	NormalSummary n6(3, 0.3);
	hha2.append(n4);
	hha2.append(n5);
	hha2.append(n6);
	DEBUG(isNormal0BestBonf(hha2, 0.05));
}

struct OCBATest
{
    mutable int nEvals;
    OCBATest(): nEvals(0){}
    int getSize()const{return 6;}
    double operator()(int i)const
    {
        ++nEvals;
        if(i == 0) return GlobalRNG().normal(1, 9);
        if(i == 1) return GlobalRNG().normal(2, 8);
        if(i == 2) return GlobalRNG().normal(3, 7);
        if(i == 3) return GlobalRNG().normal(4, 6);
        if(i == 4) return GlobalRNG().normal(5, 5);
        else return GlobalRNG().normal(6, 4);
    }
};

template<typename MULTI_FUNCTION> int
    simulateSelectBest(MULTI_FUNCTION& f, int n0 = 30, int T = 100000,
    double meanPrecision = 0, double aLevel = 0.05)
{
    int D = f.getSize(), winner = -1;
    assert(D > 1 && n0 > 1 && T > n0 * D);
    Vector<IncrementalStatistics> data(D);
    for(int i = 0; i < D; ++i)
        for(int j = 0; j < n0; ++j) data[i].addValue(f(i));
    int k = n0 * D;
    for(; k < T;)
    {
        Vector<NormalSummary> s;
        for(int i = 0; i < D; ++i) s.append(data[i].getStandardErrorSummary());
        int bestIndex = 0;
        double bestMean = s[0].mean;
        for(int i = 1; i < D; ++i)
            if(s[i].mean < bestMean) bestMean = s[bestIndex = i].mean;
        swap(s[0], s[bestIndex]);
        for(int i = 0; i < D; ++i)
            if(isPowerOfTwo(++k) && isNormal0BestBonf(s, aLevel/lgCeiling(T)))
            {
                return lgCeiling(T);
            }
        for(int i = 0; i < D; ++i) data[i].addValue(f(i));
    }
    return lgCeiling(T);
}
void testOCBA()
{
    IncrementalStatistics sO, sN;
    for(int i = 0; i < 100; ++i)
    {
        OCBATest to;
        OCBA<OCBATest> o(to);
        int nTests = o.simulateTillBest();
        sO.addValue(to.nEvals);
        pair<Vector<NormalSummary>, int> best = o.findBest();
        DEBUG(isNormal0BestBonf(best.first, 0.05/nTests));
        DEBUG(best.second);
        OCBATest t;
        simulateSelectBest(t);
        sN.addValue(t.nEvals);
    }
    DEBUG(sO.getMean());
    DEBUG(sO.getStandardErrorSummary().error95());
    DEBUG(sN.getMean());
    DEBUG(sN.getStandardErrorSummary().error95());
}

int main(int argc, char *argv[])
{
    testOCBA();
    return 0;
}

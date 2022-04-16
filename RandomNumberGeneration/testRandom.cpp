#include "Random.h"
#include "../NumericalMethods/Matrix.h"

using namespace igmdk;

void testSumHeap()
{
    int N = 1000000;
    SumHeap<double> sh;
    sh.add(0);
	sh.add(1.0/36);
	sh.add(2.0/36);
	sh.add(3.0/36);
	sh.add(4.0/36);
	sh.add(5.0/36);
	sh.add(6.0/36);
	sh.add(5.0/36);
	sh.add(4.0/36);
	sh.add(3.0/36);
	sh.add(2.0/36);
	sh.add(1.0/36);
    int sum = 0;
    for(int i = 0 ; i < N; ++i) sum += sh.next();
    DEBUG(sum*1.0/N);
}

void DDDAlias()
{
    Vector<double> probabilities;
    for(int i = 0; i < 5; ++i) probabilities.append((i-2)*(i-2)+1);
    normalizeProbs(probabilities);
    AliasMethod alias(probabilities);
    cout << "breakpoint" << endl;
}

void DDDSumHeap()
{
    Vector<double> probabilities;
    for(int i = 0; i < 5; ++i) probabilities.append((i-2)*(i-2)+1);
    normalizeProbs(probabilities);
    SumHeap<double> sumHeap;
    for(int i = 0; i < 5; ++i) sumHeap.add(probabilities[i]);
    cout << "breakpoint" << endl;
}

void testGenerators()
{
    //ARC4 x;
    MRG32k3a x;
    x.jumpAhead();
    unsigned long long N = 1 << 3;
    unsigned long long dummy = 0;
    while(N--) dummy += x.next();
    DEBUG(dummy);

    SumHeap<double> st;
    st.add(0.2);
    st.add(0.2);
    st.add(0.2);
    st.add(0.2);
    st.add(0.2);
    DEBUG(st.total());
    DEBUG(st.find(0.1));
    DEBUG(st.find(0.3));
    DEBUG(st.find(0.6));
    DEBUG(st.find(0.9));

    DEBUG(st.find(0));
    DEBUG(st.find(0));
    DEBUG(st.find(0.5));
    DEBUG(st.find(1));

    DEBUG(st.cumulative(0));
    DEBUG(st.cumulative(1));
    DEBUG(st.cumulative(2));
    DEBUG(st.cumulative(3));
    DEBUG(st.cumulative(4));
    DEBUG(st.find(st.cumulative(0)));
    DEBUG(st.find(st.cumulative(1)));
    DEBUG(st.find(st.cumulative(2)));
    DEBUG(st.find(st.cumulative(3)));
    DEBUG(st.find(st.cumulative(4)));
    for(int i = 0; i < 100; ++i)
    {
        DEBUG(st.next());
    }

    clock_t start = clock();
    MersenneTwister64 random;
    //Xorshift64 random;
    //QualityXorshift64 random;
    unsigned long long sum = 0;
    for(int i = 0; i < 1000000; ++i)
    {
		sum += random.next();
	}
	DEBUG(sum);
	clock_t end = clock();
	int time = (end - start);
    cout << "IX: " << time << endl;

    if(true)
    {
        DEBUG(GlobalRNG().uniform01());
        DEBUG(GlobalRNG().uniform(10, 20));
        DEBUG(GlobalRNG().normal01());
        DEBUG(GlobalRNG().normal(10, 20));
        DEBUG(GlobalRNG().exponential(1));
        DEBUG(GlobalRNG().gamma1(0.5));
        DEBUG(GlobalRNG().gamma1(1.5));
        DEBUG(GlobalRNG().weibull1(20));
        DEBUG(GlobalRNG().erlang(10, 2));
        DEBUG(GlobalRNG().chiSquared(10));
        DEBUG(GlobalRNG().t(10));
        DEBUG(GlobalRNG().logNormal(10, 20));
        DEBUG(GlobalRNG().beta(0.5, 0.5));
        DEBUG(GlobalRNG().F(10 ,20));
        DEBUG(GlobalRNG().cauchy(0, 1));
        DEBUG(GlobalRNG().Levy());
        DEBUG(GlobalRNG().symmetrizedLevy());
        DEBUG(GlobalRNG().binomial(0.7, 20));
        DEBUG(GlobalRNG().geometric(0.7));
        DEBUG(GlobalRNG().poisson(0.7));
        DEBUG(GlobalRNG().triangular01(0.7));
		//system("PAUSE");
	}
	int M = 100000;
	double average = 0;
	for(int i = 0; i < M; ++i)
	{
        average += GlobalRNG().beta(0.5, 0.5);
	}
	DEBUG(average/M);
}

void testMultivarNormal()
{//code in Matrix.h
    Vector<double> m(2, 0);
    Matrix<double> var = Matrix<double>::identity(2);
    MultivariateNormal mn(m, var);
    Vector<double> sample = mn.next();
    sample.debug();
}

struct XYZ
{
    double operator()()const
    {
        //return GlobalRNG().bernoulli(0.95);
        return GlobalRNG().uniform01();
    }
};
void testIncremental()
{
    XYZ xyz;
    IncrementalStatistics si = MonteCarloSimulate(xyz, 1000);
    NormalSummary s = si.getStandardErrorSummary();
    DEBUG(s.mean);
    DEBUG(s.error95());

}

int main(int argc, char *argv[])
{
    testMultivarNormal();
    return 0;
    testGenerators();
    return 0;
    DDDAlias();
    return 0;
    DDDSumHeap();
    return 0;
    testSumHeap();
    return 0;
    testIncremental();
    return 0;
}

#ifndef IGMDK_KERNEL_SVM_H
#define IGMDK_KERNEL_SVM_H
#include "ClassificationCommon.h"
#include "../Utils/Vector.h"
#include "../HashTable/LinearProbingHashTable.h"
#include "../NumericalOptimization/DiscreteGlobalOptimization.h"
#include "../RandomNumberGeneration/Statistics.h"
#include <cmath>
namespace igmdk{

struct GaussianKernel
{
    double a;
    GaussianKernel(double theA): a(theA) {}
    double operator()(NUMERIC_X const& x, NUMERIC_X const& y)const
    {//beware - if x - y is too large, result = 0
        NUMERIC_X temp = x - y;
        return exp(-a * dotProduct(temp, temp));
    }
};
template<typename KERNEL = GaussianKernel, typename X = NUMERIC_X> struct SVM
{
    Vector<X> supportVectors;
    Vector<double> supportCoefficients;
    double bias;
    KERNEL K;
    template<typename DATA> double evalK(LinearProbingHashTable<long long,
        double>& cache, long long i, long long j, DATA const& data)
    {
        long long key = i * data.getSize() + j;
        double* result = cache.find(key);
        if(result) return *result;
        else
        {
            double value = K(data.getX(i), data.getX(j));
            cache.insert(key, value);
            return value;
        }
    }
    int makeY(bool label){return label * 2 - 1;}
    double lowDiff(bool label, double C, double d){return label ? d : d + C;}
    double highDiff(bool label, double C, double d){return label ? C - d: d;}
public:
    template<typename DATA> SVM(DATA const& data, pair<KERNEL, double> const&
        params, int maxRepeats = 10, int maxConst = 10000): K(params.first)
    {
        double C = params.second;
        assert(data.getSize() > 0 && C > 0);
        bias = makeY(data.getY(0));//just in case have 1 class only
        LinearProbingHashTable<long long, double> cache;
        int n = data.getSize(), maxIters = max(maxConst, n * maxRepeats);
        Vector<double> d(n, 0), g(n);
        for(int k = 0; k < n; ++k) g[k] = makeY(data.getY(k));
        while(maxIters--)
        {//select directions using max violating pair
            int i = -1, j = -1;//i can increase, j can decrease
            for(int k = 0; k < n; ++k)
            {//find max gi and min gj
                if(highDiff(data.getY(k), C, d[k]) > 0 && (i == -1 ||
                    g[k] > g[i])) i = k;
                if(lowDiff(data.getY(k), C, d[k]) > 0 && (j == -1 ||
                    g[k] < g[j])) j = k;
            }
            if(i == -1 || j == -1) break;
            bias = (g[i] + g[j])/2;//ave for stability
            //check optimality condition
            double optGap = g[i] - g[j];
            if(optGap < 0.001) break;
            //compute direction-based minimum and box bounds
            double denom = evalK(cache, i, i, data) -
                2 * evalK(cache, i, j, data) + evalK(cache, j, j, data),
                step = min(highDiff(data.getY(i), C, d[i]),
                lowDiff(data.getY(j), C, d[j]));
            //shorten step to box bounds if needed, check for numerical
            //error in kernel calculation or duplicate data, if error
            //move points to box bounds
            if(denom > 0) step = min(step, optGap/denom);
            //update support vector coefficients and gradient
            d[i] += step;
            d[j] -= step;
            for(int k = 0; k < n; ++k) g[k] += step *
                (evalK(cache, j, k, data) - evalK(cache, i, k, data));
        }//determine support vectors
        for(int k = 0; k < n; ++k) if(abs(d[k]) > defaultPrecEps)
            {
                supportCoefficients.append(d[k]);
                supportVectors.append(data.getX(k));
            }
    }
    int predict(X const& x)const
    {
        double sum = bias;
        for(int i = 0; i < supportVectors.getSize(); ++i)
            sum += supportCoefficients[i] * K(supportVectors[i], x);
        return sum >= 0;
    }
};

template<typename KERNEL = GaussianKernel, typename X = NUMERIC_X>
class MulticlassSVM
{//need buffer for speed
    typedef pair<KERNEL, double> P;
    MulticlassLearner<BufferLearner<SVM<KERNEL, X>, InMemoryData<X, int>, P>,
        P> mcl;
public:
    template<typename DATA> MulticlassSVM(DATA const& data,
        pair<KERNEL, double> const& params): mcl(data, params) {}
    int predict(X const& x)const{return mcl.predict(x);}
};

struct NoParamsSVM
{
    MulticlassSVM<> model;
    struct CVSVMFunctor
    {
        typedef Vector<double> PARAMS;
        MulticlassSVM<> model;
        template<typename DATA> CVSVMFunctor(DATA const& data,
            PARAMS const& p):
            model(data, make_pair(GaussianKernel(p[0]), p[1])) {}
        int predict(NUMERIC_X const& x)const{return model.predict(x);}
    };
    template<typename DATA> static pair<GaussianKernel, double>
        gaussianMultiClassSVM(DATA const& data, int CLow = -5,
        int CHigh = 15, int yLow = -15, int yHi = 3)
    {
        Vector<Vector<double> > sets(2);
        for(int i = yLow; i <= yHi; i += 2) sets[0].append(pow(2, i));
        for(int i = CLow; i <= CHigh; i += 2) sets[1].append(pow(2, i));
        Vector<double> best = compassDiscreteMinimize(sets,
            SCVRiskFunctor<CVSVMFunctor, Vector<double>, DATA>(data),
            10);
        //Vector<double> best = gridMinimize(sets,
        //    SCVRiskFunctor<CVSVMFunctor, Vector<double> >(data));
        return make_pair(GaussianKernel(best[0]), best[1]);
    }
    template<typename DATA> NoParamsSVM(DATA const& data): model(data,
        gaussianMultiClassSVM(data)) {}
    int predict(NUMERIC_X const& x)const{return model.predict(x);}
};
typedef ScaledLearner<NoParamsLearner<NoParamsSVM, int>, int> SSVM;

}//end namespace
#endif


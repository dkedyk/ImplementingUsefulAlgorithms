#ifndef IGMDK_LASSO_H
#define IGMDK_LASSO_H
#include "../NumericalOptimization/NumericalOptimization.h"
#include "../NumericalOptimization/DiscreteGlobalOptimization.h"
#include "RegressionCommon.h"
#include <cmath>

namespace igmdk{

class L1LinearReg
{
    Vector<double> w;
    double b, l, r;
    int learnedCount;//this and r are only for online learning with SGD
    double f(NUMERIC_X const& x)const{return dotProduct(w, x) + b;}
    template<typename DATA> void coordinateDescent(DATA const& data,
        int maxIterations, double eps)
    {
        assert(data.getSize() > 0);
        int D = getD(data);
        Vector<double> sums(data.getSize());
        for(int i = 0; i < data.getSize(); ++i) sums[i] = data.getY(i);
        bool done = false;
        while(!done && maxIterations--)
        {
            done = true;
            for(int j = -1; j < D; ++j)
            {
                double oldVar = j == -1 ? b : w[j];
                //remove current var from sum
                for(int i = 0; i < data.getSize(); ++i)
                    sums[i] += j == -1 ? b : w[j] * data.getX(i, j);
                //solve for opt current var
                if(j == -1)
                {//update bias
                    IncrementalStatistics s;
                    for(int i = 0; i < data.getSize(); ++i)
                        s.addValue(sums[i]);
                    b = s.getMean();
                }
                else
                {//update weight
                    double a = 0, c = 0;
                    for(int i = 0; i < data.getSize(); ++i)
                    {
                        double xij = data.getX(i, j);
                        a += sums[i] * xij;
                        c += l * xij * xij;
                    }
                    if(a < -l) w[j] = (a - l)/c;
                    else if(a > l) w[j] = (a + l)/c;
                    else w[j] = 0;
                }
                //add back current var to up
                for(int i = 0; i < data.getSize(); ++i)
                    sums[i] -= j == -1 ? b : w[j] * data.getX(i, j);
                if(abs((j == -1 ? b : w[j]) - oldVar) > eps) done = false;
            }
        }
    }
public:
    template<typename DATA> L1LinearReg(DATA const& data, double theL,
        int nCoord = 1000): l(theL/2), b(0), w(getD(data)), learnedCount(-1)
        {coordinateDescent(data, nCoord, pow(10, -6));}
    typedef pair<int, pair<double, double> > PARAM;//D/l/r
    L1LinearReg(PARAM const& p): l(p.second.first/2), r(p.second.second),
        b(0), w(p.first), learnedCount(0) {}
    void learn(NUMERIC_X const& x, double y)
    {
        assert(learnedCount != -1);//can't mix batch and offline
        double rate = r * RMRate(learnedCount++), err = y - f(x);
        //l/n*|w| + (y - wx + b)2
        //dw = l/n*sign(w) - x(y - (wx + b));
        //db = - (y - (wx + b))
        for(int i = 0; i < w.getSize(); ++i) w[i] +=
            rate * (x[i] * err - (w[i] > 0 ? 1 : -1) * l/learnedCount);
        b += rate * err;
    }
    double predict(NUMERIC_X const& x)const{return f(x);}
    template<typename MODEL, typename DATA>
    static double findL(DATA const& data)
    {
        int lLow = -15, lHigh = 5;
        Vector<double> regs;
        for(double j = lHigh; j > lLow; j -= 2) regs.append(pow(2, j));
        return valMinFunc(regs.getArray(), regs.getSize(),
            RRiskFunctor<MODEL, double, DATA>(data));
    }
};
struct NoParamsL1LinearReg
{
    L1LinearReg model;
    template<typename DATA> NoParamsL1LinearReg(DATA const& data):
        model(data, L1LinearReg::findL<L1LinearReg>(data)) {}
    double predict(NUMERIC_X const& x)const{return model.predict(x);}
};
typedef ScaledLearner<NoParamsLearner<NoParamsL1LinearReg, double>, double>
    SLasso;

class SRaceLasso
{
    ScalerMQ s;
    RaceLearner<L1LinearReg, L1LinearReg::PARAM> model;
    static Vector<L1LinearReg::PARAM> makeParams(int D)
    {
        Vector<L1LinearReg::PARAM> result;
        int lLow = -15, lHigh = 5, rLow = -15, rHigh = 5;
        for(int j = lHigh; j > lLow; j -= 2)
        {
            double l = pow(2, j);
            for(int i = rHigh; i > rLow; i -= 2) result.append(
                L1LinearReg::PARAM(D, pair<double, double>(l, pow(2, i))));
        }
        return result;
    }
public:
    template<typename DATA> SRaceLasso(DATA const& data):
        model(makeParams(getD(data))), s(getD(data))
    {
        for(int j = 0; j < 1000000; ++j)
        {
            int i = GlobalRNG().mod(data.getSize());
            learn(data.getX(i), data.getY(i));
        }
    }
    SRaceLasso(int D): model(makeParams(D)), s(D) {}
    void learn(NUMERIC_X const& x, double y)
    {
        s.addSample(x);
        model.learn(s.scale(x), y);
    }
    double predict(NUMERIC_X const& x)const
        {return model.predict(s.scale(x));}
};

}//end namespace
#endif


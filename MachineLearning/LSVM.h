#ifndef IGMDK_LSVM_H
#define IGMDK_LSVM_H
#include "ClassificationCommon.h"
#include "../Utils/Vector.h"
#include "../NumericalOptimization/NumericalOptimization.h"
#include <cmath>
namespace igmdk{

class BinaryLSVM
{
    Vector<double> w;
    double b, l;
    int learnedCount;
    static int y(bool label){return label * 2 - 1;}
    static double loss(double fxi, double yi){return max(0.0, 1 - fxi * yi);}
    double f(NUMERIC_X const& x)const{return dotProduct(w, x) + b;}
    template<typename DATA> class GSL1Functor
    {
        DATA const& data;
        mutable Vector<double> sums;
        Vector<double> &w;
        double &b, l;
        int j, D;
        mutable int evalCount;
        double getSumI(double wj, int i)const
        {
            return sums[i] + (j == D ? wj - b : (wj - w[j]) * data.getX(i, j));
        }
    public:
        GSL1Functor(DATA const& theData, double& theB, Vector<double>& theW,
            double theL): data(theData), sums(theData.getSize(), theB),
            w(theW), b(theB), l(theL), j(0), D(getD(data)), evalCount(0)
        {
            for(int i = 0; i < data.getSize(); ++i)
                sums[i] += dotProduct(w, data.getX(i)) + b;
        }
        void setCurrentDimension(int theJ)
        {
            assert(theJ >= 0 && theJ < D + 1);
            j = theJ;
        }
        int getEvalCount()const{return evalCount;}
        int getSize()const{return D + 1;}
        Vector<double> getX()const
        {
            Vector<double> x = w;
            x.append(b);
            return x;
        }
        double getXi()const{return j == D ? b : w[j];}
        double operator()(double wj)const
        {
            ++evalCount;
            double result = j == D ? 0 : l * abs(wj);
            for(int i = 0; i < data.getSize(); ++i)
                result += loss(getSumI(wj, i), y(data.getY(i)));
            return result/data.getSize();
        }
        void bind(double wjNew)
        {//first update sums
            for(int i = 0; i < data.getSize(); ++i)
                sums[i] = getSumI(wjNew, i);
            (j == D ? b : w[j]) = wjNew;
        }
    };
public:
    BinaryLSVM(pair<int, double> const& p): w(p.first), l(p.second), b(0),
        learnedCount(0) {}
    template<typename DATA> BinaryLSVM(DATA const& data, double theL,
        int nGoal = 100000, int nEvals = 100000): l(theL), b(0),
        w(getD(data)), learnedCount(0)
    {//first SGD
        for(int j = 0; j < ceiling(nGoal, data.getSize()); ++j)
            for(int i = 0; i < data.getSize(); ++i)
                learn(data.getX(i), data.getY(i), data.getSize());
        //then coordinate descent
        GSL1Functor<DATA> f(data, b, w, l);
        unimodalCoordinateDescent(f, nEvals);
    }
    int getLearnedCount(){return learnedCount;}
    void learn(NUMERIC_X const& x, int label, int n = -1)
    {//online mode uses SGD only
        if(n == -1) n = learnedCount + 1;
        double rate = RMRate(learnedCount++), yl = y(label);
        for(int i = 0; i < w.getSize(); ++i)
            w[i] -= rate * (w[i] > 0 ? 1 : -1) * l/n;
        if(yl * f(x) < 1)
        {
            w -= x * (-yl * rate);
            b += rate * yl;
        }
    }
    int predict(NUMERIC_X const& x)const{return f(x) >= 0;}
    template<typename MODEL, typename DATA>
    static double findL(DATA const& data)
    {//used for regression as well
        int lLow = -15, lHigh = 5;
        Vector<double> regs;
        for(double j = lHigh; j > lLow; j -= 2) regs.append(pow(2, j));
        return valMinFunc(regs.getArray(), regs.getSize(),
            SCVRiskFunctor<MODEL, double, DATA>(data));
    }
};

struct NoParamsLSVM
{
    typedef MulticlassLearner<BinaryLSVM, double> MODEL;
    MODEL model;
    template<typename DATA> NoParamsLSVM(DATA const& data): model(data,
        BinaryLSVM::findL<MODEL, DATA>(data)) {}
    int predict(NUMERIC_X const& x)const{return model.predict(x);}
};
typedef ScaledLearner<NoParamsLearner<NoParamsLSVM, int>, int> SLSVM;

class SRaceLSVM
{
    ScalerMQ s;
    typedef pair<int, double> P;
    RaceLearner<OnlineMulticlassLearner<BinaryLSVM, P>, P> model;
    static Vector<P> makeParams(int D)
    {
        Vector<P> result;
        int lLow = -15, lHigh = 5;
        for(int j = lHigh; j > lLow; j -= 2)
        {
            double l = pow(2, j);
            result.append(P(D, l));
        }
        return result;
    }
public:
    template<typename DATA> SRaceLSVM(DATA const& data):
        model(makeParams(getD(data))), s(getD(data))
    {
        for(int j = 0; j < 1000000; ++j)
        {
            int i = GlobalRNG().mod(data.getSize());
            learn(data.getX(i), data.getY(i));
        }
    }
    SRaceLSVM(int D): model(makeParams(D)), s(D){}
    void learn(NUMERIC_X const& x, int label)
    {
        s.addSample(x);
        model.learn(s.scale(x), label);
    }
    int predict(NUMERIC_X const& x)const{return model.predict(s.scale(x));}
};

}//end namespace
#endif


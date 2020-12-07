#include <cassert>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <limits>
#include <memory> //for shared ptr
#include "GlobalNumericalOptimization.h"
#include "testOptCommon.h"
using namespace std;
using namespace igmdk;

struct TestFunctionsGlobalBoxMin
{
    struct BaseF
    {
        virtual double operator()(Vector<double> const& x)const = 0;
        virtual string name()const = 0;
        virtual Vector<pair<double, double> > getBox()const = 0;
        virtual Vector<double> getAnswer()const = 0;
    };
    //separability annotations from Jamil & Yang
    //many visualizations at
    //1. http://infinity77.net/global_optimization/test_functions.html
    //2. https://www.sfu.ca/~ssurjano/optimization.html
    //3. http://al-roomi.org/benchmarks/unconstrained
    struct Ackley: public BaseF
    {//Quadratic with noise; from Simon; non-separable
        int n;
        Ackley(int theN = 2): n(theN) {assert(theN % 2 == 0);}
        double operator()(Vector<double> const& x)const
        {
            double temp = 0;
            for(int i = 0; i < n; ++i)
                temp += cos(2 * PI() * x[i]);
            return 20 + exp(1) -20 * exp(-0.2 * norm(x)) - exp(temp/n);
        }
        string name()const{return "Ackley" + to_string(n);}
        Vector<pair<double, double> > getBox()const
            {return Vector<pair<double, double> >(n, make_pair(-30, 30));}
        Vector<double> getAnswer()const{return Vector<double>(n, 0);}
    };
    struct FletcherPowell: public BaseF
    {//Unpredictable landscape; from Simon; non-separable (self-concluded)
        int n;
        Vector<double> a;
        Vector<Vector<double> > aij, bij;
        FletcherPowell(int theN = 2): n(theN), a(n), aij(n, Vector<double>(n)),
            bij(aij)
        {
            assert(n > 1);
            for(int i = 0; i < n; ++i)
            {
                a[i] = GlobalRNG().uniform(-PI(), PI());
                for(int j = 0; j < n; ++j)
                {
                    aij[i][j] = GlobalRNG().uniform(-100, 100);
                    bij[i][j] = GlobalRNG().uniform(-100, 100);
                }
            }
        }
        double operator()(Vector<double> const& x)const
        {
            double sum = 0;
            for(int i = 0; i < n; ++i)
            {
                double Ai = 0, Bi = 0;
                for(int j = 0; j < n; ++j)
                {
                    Ai += aij[i][j] * sin(a[j]) + bij[i][j] * cos(a[j]);
                    Bi += aij[i][j] * sin(x[j]) + bij[i][j] * cos(x[j]);
                }
                sum += (Ai - Bi) * (Ai - Bi);
            }
            return sum;
        }
        string name()const{return "FletcherPowell" + to_string(n);}
        Vector<pair<double, double> > getBox()const
        {
            return Vector<pair<double, double> >(n, make_pair(-PI(), PI()));
        }
        Vector<double> getAnswer()const{return a;}
    };
    struct Griewank: public BaseF
    {//Very noisy non-separable; from Simon; non-separable
        int n;
        Griewank(int theN = 2): n(theN) {assert(theN % 2 == 0);}
        double operator()(Vector<double> const& x)const
        {
            double temp = 1 + dotProduct(x, x)/4000, prod = 1;
            for(int i = 0; i < n; ++i) prod *= cos(x[i]/sqrt(i + 1));
            return temp - prod;
        }
        string name()const{return "Griewank" + to_string(n);}
        Vector<pair<double, double> > getBox()const
            {return Vector<pair<double, double> >(n, make_pair(-600, 600));}
        Vector<double> getAnswer()const{return Vector<double>(n, 0);}
    };
    struct Rastrigin: public BaseF
    {//Exponentially many in n local minima; from Simon; separable
        int n;
        Rastrigin(int theN = 2): n(theN) {assert(theN % 2 == 0);}
        double operator()(Vector<double> const& x)const
        {
            double temp = 0;
            for(int i = 0; i < n; ++i)
                temp += x[i] * x[i] - 10 * cos(2 * PI() * x[i]);
            return 10 * n + temp;
        }
        string name()const{return "Rastrigin" + to_string(n);}
        Vector<pair<double, double> > getBox()const
            {return Vector<pair<double, double> >(n, make_pair(-5.12, 5.12));}
        Vector<double> getAnswer()const{return Vector<double>(n, 0);}
    };
    struct SchwefelDoubleSum: public BaseF
    {//High condition number; from Simon; non-separable
        int n;
        SchwefelDoubleSum(int theN = 2): n(theN) {assert(theN % 2 == 0);}
        double operator()(Vector<double> const& x)const
        {
            double sum = 0;
            for(int i = 0; i < n; ++i)
            {
                double sum2 = 0;
                for(int j = 0; j < i; ++j) sum2 += x[j];
                sum += sum2 * sum2;
            }
            return sum;
        }
        string name()const{return "SchwefelDoubleSum" + to_string(n);}
        Vector<pair<double, double> > getBox()const
            {return Vector<pair<double, double> >(n, make_pair(-65.356, 65.356));}
        Vector<double> getAnswer()const{return Vector<double>(n, 0);}
    };
    struct SchwefelMax: public BaseF
    {//Only largest value matters; from Simon; separable
        int n;
        SchwefelMax(int theN = 2): n(theN) {assert(theN % 2 == 0);}
        double operator()(Vector<double> const& x)const{return normInf(x);}
        string name()const{return "SchwefelMax" + to_string(n);}
        Vector<pair<double, double> > getBox()const
            {return Vector<pair<double, double> >(n, make_pair(-100, 100));}
        Vector<double> getAnswer()const{return Vector<double>(n, 0);}
    };
    struct StepFunction: public BaseF
    {//Quadratic with plateaus; from Simon; separable
        int n;
        StepFunction(int theN = 2): n(theN) {assert(theN % 2 == 0);}
        double operator()(Vector<double> const& x)const
        {
            double sum = 0;
            for(int i = 0; i < n; ++i)
            {
                double temp = int(x[i] + 0.5);
                sum += temp * temp;
            }
            return sum;
        }
        string name()const{return "StepFunction" + to_string(n);}
        Vector<pair<double, double> > getBox()const
            {return Vector<pair<double, double> >(n, make_pair(-100, 100));}
        Vector<double> getAnswer()const{return Vector<double>(n, 0);}
    };
    struct Weierstrass: public BaseF
    {//Nowhere differentiable quadratic; from Simon; separable (!? maybe typo in paper)
        int n;
        Weierstrass(int theN = 2): n(theN) {assert(theN % 2 == 0);}
        double operator()(Vector<double> const& x)const
        {
            int kMax = 20;
            double a = 0.5, b = 3, sum = 0, sum2 = 0;
            for(int k = 0; k <= kMax; ++k)
                sum2 += pow(a, k) * cos(PI() * pow(b, k));
            for(int i = 0; i < n; ++i)
            {
                double sum3 = 0;
                for(int k = 0; k <= kMax; ++k) sum3 += pow(a, k) *
                    cos(2 * PI() * pow(b, k) * (x[i] + 0.5));
                sum += sum3;
            }
            return sum - n * sum2;
        }
        string name()const{return "Weierstrass" + to_string(n);}
        Vector<pair<double, double> > getBox()const
            {return Vector<pair<double, double> >(n, make_pair(-5, 5));}
        Vector<double> getAnswer()const{return Vector<double>(n, 0);}
    };
    struct Bukin6: public BaseF
    {//Interesting structure; from Jamil & Yang; non-separable
        double operator()(Vector<double> const& x)const
        {
            assert(x.getSize() == 2);
            return 100 * sqrt(abs(x[1] - 0.01 * x[0] * x[0])) +
                0.01 * abs(x[0] + 10);
        }
        string name()const{return "Bukin6";}
        Vector<pair<double, double> > getBox()const
        {
            Vector<pair<double, double> > box(2);
            box[0] = make_pair(-15, -5);
            box[1] = make_pair(-3, 3);
            return box;
        }
        Vector<double> getAnswer()const
        {
            Vector<double> answer(2);
            answer[0] = -10;
            answer[1] = 1;
            return answer;
        }
    };
    struct Damavandi: public BaseF
    {//DeceptiveQuadratic; from Jamil & Yang; non-separable
        double operator()(Vector<double> const& x)const
        {
            double temp = 1, sum = 0;
            for(int i = 0; i < 2; ++i)
            {
                sum += 2 + (x[i] - 7) * (x[i] - 7);
                double temp2 = PI() * (x[i] - 2);
                temp *= sin(temp2)/temp2;
            }
            return (1 + temp * temp * temp * temp * temp) * sum;
        }
        string name()const{return "Damavandi";}
        Vector<pair<double, double> > getBox()const
            {return Vector<pair<double, double> >(2, make_pair(0, 14));}
        Vector<double> getAnswer()const
            {return Vector<double>(2, 2);}
    };
    struct Easom: public BaseF
    {//Interesting structure; from Jamil & Yang; non-separable
        double operator()(Vector<double> const& x)const
        {
            assert(x.getSize() == 2);
            double temp1 = x[0] - PI(), temp2 = x[1] - PI();
            return -cos(x[0]) * cos(x[1]) * exp(-temp1 * temp1 -temp2 * temp2);
        }
        string name()const{return "Easom";}
        Vector<pair<double, double> > getBox()const
            {return Vector<pair<double, double> >(2, make_pair(-100, 100));}
        Vector<double> getAnswer()const
            {return Vector<double>(2, PI());}
    };
    struct GulfRND: public BaseF
    {//From More et al; paper has mistake "mi" should be "-" not m * (i + 1)
    //Real-world function; box from Jamil & Yang; non-separable
        double operator()(Vector<double> const& x)const
        {
            assert(x.getSize() == 3);
            Vector<double> fx(3);
            for(int i = 0; i < fx.getSize(); ++i)
            {
                double ti = (i + 1.0)/100, yi = 25 + pow(-50 * log(ti), 2.0/3);
                fx[i] = exp(-pow(abs(yi - x[1]), x[2])/x[0]) - ti;
            }
            return dotProduct(fx, fx);
        }
        string name()const{return "GulfRND";}
        Vector<pair<double, double> > getBox()const
        {
            Vector<pair<double, double> > box(3);
            box[0] = make_pair(0.1, 100);
            box[1] = make_pair(0, 25.6);
            box[2] = make_pair(0, 5);
            return box;
        }
        Vector<double> getX0()const
        {
            Vector<double> x0(3);
            x0[0] = 5;
            x0[1] = 2.5;
            x0[2] = 0.15;
            return x0;
        }
        Vector<double> getAnswer()const{return getX0() * 10;}
    };
    struct Price2: public BaseF
    {//Almost flat global landscape + noisy; from Jamil & Yang; non-separable
        double operator()(Vector<double> const& x)const
        {
            assert(x.getSize() == 2);
            double temp1 = sin(x[0]), temp2 = sin(x[1]);
            return 1 + temp1 * temp1 + temp2 * temp2 - 0.1 * exp(-dotProduct(x, x));
        }
        string name()const{return "Price2";}
        Vector<pair<double, double> > getBox()const
            {return Vector<pair<double, double> >(2, make_pair(-10, 10));}
        Vector<double> getAnswer()const{return Vector<double>(2, 0);}
    };
    struct Trefethen: public BaseF
    {//Interesting structure; from Jamil & Yang; non-separable
     //High-prec answer from Bornemann et al
        double operator()(Vector<double> const& x)const
        {
            assert(x.getSize() == 2);
            return exp(sin(50 * x[0])) + sin(60 * exp(x[1])) +
                sin(70 * sin(x[0])) + sin(sin(80 * x[1])) -
                sin(10 * (x[0] + x[1])) + dotProduct(x, x)/4;
        }
        string name()const{return "Trefethen";}
        Vector<pair<double, double> > getBox()const
            {return Vector<pair<double, double> >(2, make_pair(-10, 10));}
        Vector<double> getAnswer()const
        {
            Vector<double> answer(2);
            answer[0] = -0.024403079743737212;
            answer[1] = 0.21061242727591697;
            return answer;
        }
    };
    struct Trig: public BaseF
    {//From Dennis & Schnabel, they give no solution, but x = 0 works
    //Box from Jamil & Yang; non-separable
        int n;
        Trig(int theN = 2): n(theN) {}
        double operator()(Vector<double> const& x)const
        {
            double cSum = 0;
            for(int i = 0; i < n; ++i) cSum += cos(x[i]);
            Vector<double> fx(n, n - cSum);
            for(int i = 0; i < n; ++i)
                fx[i] += (i + 1) * (1 - cos(x[i])) - sin(x[i]);
            return dotProduct(fx, fx);
        }
        string name()const{return "Trig" + to_string(n);}
        Vector<pair<double, double> > getBox()const
            {return Vector<pair<double, double> >(n, make_pair(0, PI()));}
        Vector<double> getAnswer()const{return Vector<double>(n, 0);}
    };
    struct Pinter: public BaseF
    {//Quadratic with much noise; From Jamil & Yang; non-separable
        int n;
        Pinter(int theN = 2): n(theN) {}
        double operator()(Vector<double> const& x)const
        {
            double sum = 0;
            for(int i = 0; i < n; ++i)
            {
                sum += (i + 1) + x[i] * x[i];
                double xm1 = x[(i - 1 + n) % n], xp1 = x[(i + 1) % n],
                    A = xm1 * sin(x[i]) + sin(xp1), temp = sin(A),
                    B = xm1 * xm1 - 2 * x[i] + 3 * xp1 - cos(x[i]) + 1;
                sum += 20 * (i + 1) * temp * temp;
                sum += (i + 1) * log10(1 + (i + 1) * B * B);
            }
            return sum;
        }
        string name()const{return "Pinter" + to_string(n);}
        Vector<pair<double, double> > getBox()const
            {return Vector<pair<double, double> >(n, make_pair(-10, 10));}
        Vector<double> getAnswer()const{return Vector<double>(n, 0);}
    };
    struct Salomon: public BaseF
    {//NoisyQuadratic; From Jamil & Yang; non-separable
        int n;
        Salomon(int theN = 2): n(theN) {}
        double operator()(Vector<double> const& x)const
        {
            double temp = norm(x);
            return 1 - cos(2 * PI() * temp) + 0.1 * temp;
        }
        string name()const{return "Salomon" + to_string(n);}
        Vector<pair<double, double> > getBox()const
            {return Vector<pair<double, double> >(n, make_pair(-100, 100));}
        Vector<double> getAnswer()const{return Vector<double>(n, 0);}
    };
    struct SchaeferF6: public BaseF
    {//Oscilatory; From Jamil & Yang; non-separable
        int n;
        SchaeferF6(int theN = 2): n(theN) {}
        double operator()(Vector<double> const& x)const
        {
            double sum = 0;
            for(int i = 0; i < n - 1; ++i)
            {
                double temp = x[i] * x[i] + x[i + 1] * x[i + 1],
                    temp1 = sin(sqrt(temp)),
                    temp2 = 1 + 0.001 * temp;
                sum += 0.5 + (temp1 * temp1 - 0.5)/(temp2 * temp2);
            }
            return sum;
        }
        string name()const{return "SchaeferF6" + to_string(n);}
        Vector<pair<double, double> > getBox()const
            {return Vector<pair<double, double> >(n, make_pair(-100, 100));}
        Vector<double> getAnswer()const{return Vector<double>(n, 0);}
    };
    static int evalCount;
    struct MetaF
    {
        shared_ptr<BaseF> f;
        template<typename F> MetaF(shared_ptr<F> const& theF): f(theF){}
        double operator()(Vector<double> const& x)const
        {
            ++evalCount;
            return (*f)(x);
        }
        string getName()const{return f->name();}
        Vector<pair<double, double> > getBox()const{return f->getBox();}
        pair<Vector<double>, double> getAnswer()const
            {return make_pair(f->getAnswer(), (*f)(f->getAnswer()));}
    };
    static Vector<MetaF> getFunctions()
    {
        Vector<MetaF> result;
        result.append(MetaF(make_shared<Bukin6>()));
        result.append(MetaF(make_shared<Damavandi>()));
        result.append(MetaF(make_shared<Easom>()));
        result.append(MetaF(make_shared<GulfRND>()));
        result.append(MetaF(make_shared<Price2>()));
        result.append(MetaF(make_shared<Trefethen>()));
        int ds[] = {2, 6, 10, 30, 100};
        for(int i = 0; i < sizeof(ds)/sizeof(ds[0]); ++i)
        {
            result.append(MetaF(make_shared<Ackley>(ds[i])));
            result.append(MetaF(make_shared<FletcherPowell>(ds[i])));//very expensive omit for now
            result.append(MetaF(make_shared<Griewank>(ds[i])));
            result.append(MetaF(make_shared<Rastrigin>(ds[i])));
            result.append(MetaF(make_shared<SchwefelDoubleSum>(ds[i])));
            result.append(MetaF(make_shared<SchwefelMax>(ds[i])));
            result.append(MetaF(make_shared<StepFunction>(ds[i])));
            result.append(MetaF(make_shared<Weierstrass>(ds[i])));//very expensive omit for now
            result.append(MetaF(make_shared<Trig>(ds[i])));
            result.append(MetaF(make_shared<Pinter>(ds[i])));
            result.append(MetaF(make_shared<Salomon>(ds[i])));
            result.append(MetaF(make_shared<SchaeferF6>(ds[i])));
        }
        return result;
    }
};
int TestFunctionsGlobalBoxMin::evalCount = 0;

template<typename FUNCTION, typename SAMPLER> pair<Vector<double>, double>
    randomMinimize(FUNCTION const& f, SAMPLER const& s, int maxEvals = 1000000)
{
    assert(maxEvals > 0);
    pair<Vector<double>, double> best;
    for(int i = 0; i < maxEvals; ++i)
    {
        Vector<double> xNext = s();
        double yNext = f(xNext);
        if(i == 0 || !isfinite(best.second) || yNext < best.second)
        {
            best.first = xNext;
            best.second = yNext;
        }
    }
    return best;
}
template<typename FUNCTION, typename SAMPLER> pair<Vector<double>, double>
    randomBeforeLocalMinimize(FUNCTION const& f, SAMPLER const& s,
    int maxEvals = 1000000, double yPrecision = highPrecEps)
{
    assert(maxEvals > 1);
    return hybridLocalMinimize(randomMinimize(f, s, maxEvals * 0.9).first,
        f, maxEvals * 0.1, yPrecision);
}

/*
class ILSHybrid
{
    typedef pair<Vector<double>, double> P;
    template<typename FUNCTION, typename M_SAMPLER> struct Move
    {
        FUNCTION const& f;
        M_SAMPLER s;
        int maxLocalEvals;
        typedef Vector<double> X;
        X localSearchBest(X const& x){return Cholesky11CMA_ESMinimize(f, x,
            box, maxLocalEvals9).first;}
        double getScore(X const& x){return f(x);}
        void bigMove(X& x)
        {//base jump on scale of best point
            X xOld = x;
            x = s(x);//don't want escape
            if(!isfinite(norm(x))) x = xOld;
        }
    };
public:
    template<typename FUNCTION, typename M_SAMPLER> static P minimize(
        Vector<double> const& x0, M_SAMPLER const& s,
        FUNCTION const& f, int maxEvals = 1000000)
    {
        --maxEvals;
        int maxLocalEvals = sqrt(maxEvals);
        Move<FUNCTION, M_SAMPLER> move = {f, s, maxLocalEvals};
        Vector<double> x = iteratedLocalSearch(move, x0,
            maxEvals/maxLocalEvals);
        return make_pair(x, f(x));
    }
};*/

template<typename FUNCTION, typename SAMPLER> pair<Vector<double>, double>
    differentialEvolutionBeforeLocalMinimize(FUNCTION const& f, SAMPLER const& s,
    Vector<pair<double, double> > const& box,
    int maxEvals = 1000000, double yPrecision = highPrecEps)
{
    assert(maxEvals > 1);
    return hybridLocalMinimize(differentialEvolutionMinimize(f, s, box,
        maxEvals * 0.9).first, f, maxEvals * 0.1, yPrecision);
}

template<typename FUNCTION, typename SAMPLER> pair<Vector<double>, double>
    geneticLocalSearchBeforeLocalMinimize(FUNCTION const& f, SAMPLER const&s,
    int maxEvals = 1000000, double yPrecision = highPrecEps)
{
    assert(maxEvals > 1);
    return hybridLocalMinimize(geneticLocalSearchContinuous(f, s,
        maxEvals * 0.9).first, f, maxEvals * 0.1, yPrecision);
}

//can define more general constraint then box later if needed
template<typename FUNCTION> pair<Vector<double>, double>
    Cholesky11CMA_ESMinimize(FUNCTION const& f, Vector<double> initialGuess,
    Vector<pair<double, double> > const& box, int maxEvals = 1000000)
{
    assert(maxEvals > 0);
    pair<Vector<double>, double> xy(initialGuess, f(initialGuess));
    --maxEvals;
    int D = initialGuess.getSize();
    Matrix<double> A = Matrix<double>::identity(D);
    double scale = max(1.0, norm(xy.first))/10, pst = 2.0/11, pt = 11.0/25,
        ca = sqrt(1 - 2.0/(D * D + 6)), cp = 1.0/12, d = 1 + 1.0/D,
        psAve = pst;
    while(maxEvals-- > 0)
    {
        Vector<double> z(D), xNext = xy.first;
        for(int i = 0; i < D; ++i) z[i] += GlobalRNG().normal01();
        xNext += A * z * scale;
        xNext = boxTrim(xNext, box);
        double yNext = f(xNext), ls = 0;
        if(yNext <= xy.second) ls = 1;
        psAve = (1 - cp) * psAve + cp * ls;
        scale *= exp((psAve - pst * (1 - psAve)/(1 - pst))/d);
        if(yNext <= xy.second)
        {
            xy.first = xNext;
            xy.second = yNext;
            if(psAve <= pt)
            {
                double z2 = dotProduct(z, z), ca2 = ca * ca;
                A = ca * (A + (sqrt(1 + (1 - ca2) * z2/ca2) - 1)/z2 *
                    outerProduct(A * z, z));
            }
        }
    }
    return xy;
}
template<typename FUNCTION> pair<Vector<double>,
    double> Cholesky11CMA_ESBeforeLocalMinimize(FUNCTION const& f,
    Vector<double> initialGuess, Vector<pair<double, double> > const& box,
    int maxEvals = 1000000, double yPrecision = highPrecEps)
{
    assert(maxEvals > 1);
    return hybridLocalMinimize(Cholesky11CMA_ESMinimize(f, initialGuess, box,
        maxEvals * 0.9).first, f, maxEvals * 0.1, yPrecision);
}

template<typename FUNCTION, typename M_SAMPLER> pair<Vector<double>,
    double> simulatedAnnealingMinimize(FUNCTION const& f, M_SAMPLER const& ms,
    Vector<double> x0, int maxEvals = 1000000)
{
    assert(maxEvals > 0 && isfinite(norm(x0)));
    Vector<double> x = selfTunedSimulatedAnnealing(ContinuousMMove<FUNCTION,
        M_SAMPLER>(f, ms), x0, maxEvals);
    return make_pair(x, f(x));
}
template<typename FUNCTION, typename M_SAMPLER> pair<Vector<double>,
    double> simulatedAnnealingBeforeLocalMinimize(FUNCTION const& f,
    M_SAMPLER const& s, Vector<double> initialGuess,
    int maxEvals = 1000000, double yPrecision = highPrecEps)
{
    assert(maxEvals > 1);
    return hybridLocalMinimize(simulatedAnnealingMinimize(f, s, initialGuess,
        maxEvals * 0.9).first, f, maxEvals * 0.1, yPrecision);
}

template<typename FUNCTION, typename SAMPLER> pair<Vector<double>,
    double> simulatedAnnealingSMinimize(FUNCTION const& f,
    SAMPLER const& s, Vector<double> x0, int maxEvals = 1000000)
{
    assert(maxEvals > 0 && isfinite(norm(x0)));
    IncrementalWrapper<FUNCTION> iw(f, x0);
    return simulatedAnnealingSMinimizeIncremental(iw, s, x0, maxEvals);
}
template<typename FUNCTION, typename SAMPLER> pair<Vector<double>,
    double> simulatedAnnealingSBeforeLocalMinimize(FUNCTION const& f,
    SAMPLER const& s, Vector<double> initialGuess,
    int maxEvals = 1000000, double yPrecision = highPrecEps)
{
    assert(maxEvals > 1);
    return hybridLocalMinimize(simulatedAnnealingSMinimize(f, s, initialGuess,
        maxEvals * 0.9).first, f, maxEvals * 0.1, yPrecision);
}

template<typename FUNCTION, typename M_SAMPLER> pair<Vector<double>,
    double> markovianBeforeLocalMinimize(FUNCTION const& f,
    M_SAMPLER const& s, Vector<double> initialGuess,
    int maxEvals = 1000000, double yPrecision = highPrecEps)
{
    assert(maxEvals > 1);
    return hybridLocalMinimize(markovianMinimize(f, s, initialGuess,
        maxEvals * 0.9).first, f, maxEvals * 0.1, yPrecision);
}


template<typename FUNCTION, typename SAMPLER> pair<Vector<double>, double>
    RCDBeforeLocalMinimize(FUNCTION const& f, SAMPLER const& s,
    int maxEvals = 1000000, double yPrecision = highPrecEps)
{
    return hybridLocalMinimize(RCDGeneral(f, s, s(),
        maxEvals * 0.9).first, f, maxEvals * 0.1, yPrecision);
}

template<typename FUNCTION, typename SAMPLER> pair<Vector<double>, double>
    incrementalSABeforeLocalMinimizeGeneral(FUNCTION const &f,
    SAMPLER const& s, int maxEvals = 1000000)
{
    IncrementalWrapper<FUNCTION> iw(f, s());
    double y = incrementalSABeforeLocalMinimize(iw, s, maxEvals);
    return make_pair(iw.xBound, y);
}

template<typename POINT, typename FUNCTION> void debugResultGlobalBox(
    Vector<pair<POINT, double> > const& results,
    FUNCTION const& f, Vector<Vector<string> > & matrix, int start)
{
    debugResultHelperBatch<TestFunctionsGlobalBoxMin>(results, f, matrix, start);
}
template<typename TESTCASE, typename SAMPLER, typename STEP_SAMPLER>
    void testAllSolversGlobalHelper(TESTCASE const& f, SAMPLER const& s,
    STEP_SAMPLER const& ms, bool useBox, bool smallBudget, int nRepeats,
    Vector<Vector<string> >& matrix)
{
    int D = f.getAnswer().first.getSize(), start = 0;
    BoxSampler sb(f.getBox());
    Vector<Vector<double> > initialGuesses;
    for(int i = 0; i < nRepeats; ++i) initialGuesses.append(sb());

    Vector<pair<double, double> > box = useBox ? f.getBox() : makeAgnosticBox(D);

    Vector<pair<Vector<double>, double> > results(nRepeats);

    if(smallBudget)
    {
        int nEvals = 100;
        DEBUG("small eval case");
        //benchmarks
        DEBUG("Random");
        matrix.lastItem().append("Random");
        start = clock();
        for(int i = 0; i < nRepeats; ++i) results[i] = randomMinimize(f, s, nEvals);
        debugResultGlobalBox(results, f, matrix, start);
        //better samplers
        DEBUG("RCD");
        matrix.lastItem().append("RCD");
        start = clock();
        for(int i = 0; i < nRepeats; ++i) results[i] = RCDGeneral(f, s, initialGuesses[i], nEvals);
        debugResultGlobalBox(results, f, matrix, start);
        DEBUG("Markovian");
        matrix.lastItem().append("Markovian");
        start = clock();
        for(int i = 0; i < nRepeats; ++i) results[i] = markovianMinimize(f, ms, initialGuesses[i], nEvals);
        debugResultGlobalBox(results, f, matrix, start);
        DEBUG("SmallBudgetHybrid");
        matrix.lastItem().append("SmallBudgetHybrid");
        start = clock();
        for(int i = 0; i < nRepeats; ++i) results[i] = smallBudgetHybridMinimize(f, s, ms, initialGuesses[i], nEvals);
        debugResultGlobalBox(results, f, matrix, start);
        return;
    }
    else DEBUG("large eval case");

    //benchmarks
    DEBUG("RandomBLM");
    matrix.lastItem().append("RandomBLM");
    start = clock();
    for(int i = 0; i < nRepeats; ++i) results[i] = randomBeforeLocalMinimize(f, s);
    debugResultGlobalBox(results, f, matrix, start);
    //better samplers
    DEBUG("RCDBLM");
    matrix.lastItem().append("RCDBLM");
    start = clock();
    for(int i = 0; i < nRepeats; ++i) results[i] = RCDBeforeLocalMinimize(f, s);
    debugResultGlobalBox(results, f, matrix, start);
    DEBUG("IncrementalSABLM");
    matrix.lastItem().append("IncrementalSABLM");
    start = clock();
    for(int i = 0; i < nRepeats; ++i) results[i] = incrementalSABeforeLocalMinimizeGeneral(f, s);
    debugResultGlobalBox(results, f, matrix, start);
    DEBUG("MarkovianBLM");
    matrix.lastItem().append("MarkovianBLM");
    start = clock();
    for(int i = 0; i < nRepeats; ++i) results[i] = markovianBeforeLocalMinimize(f, ms, initialGuesses[i]);
    debugResultGlobalBox(results, f, matrix, start);
    DEBUG("HybridBLM");
    matrix.lastItem().append("HybridBLM");
    start = clock();
    for(int i = 0; i < nRepeats; ++i) results[i] = hybridBeforeLocalMinimize(f, s, ms, box);
    debugResultGlobalBox(results, f, matrix, start);
    //metaheuristics
    /*DEBUG("ILS");MUST FIRST PASS BOX AND FIX CMA TO TRY THIS
    matrix.lastItem().append("ILS");
    start = clock();
    for(int i = 0; i < nRepeats; ++i) results[i] = ILSHybrid::minimize(initialGuesses[i], ms, f);
    debugResultGlobalBox(results, f, matrix, start);*/
    //only test with box for now
    DEBUG("Chol11CMA_ESBLM");
    matrix.lastItem().append("Chol11CMA_ESBLM");
    start = clock();
    for(int i = 0; i < nRepeats; ++i) results[i] = Cholesky11CMA_ESBeforeLocalMinimize(f, initialGuesses[i], f.getBox());
    debugResultGlobalBox(results, f, matrix, start);
    DEBUG("SimulatedAnnealingBLM");
    matrix.lastItem().append("SimulatedAnnealingBLM");
    start = clock();
    for(int i = 0; i < nRepeats; ++i) results[i] = simulatedAnnealingBeforeLocalMinimize(f, ms, initialGuesses[i]);
    debugResultGlobalBox(results, f, matrix, start);
    DEBUG("SimulatedAnnealingSBLM");
    matrix.lastItem().append("SimulatedAnnealingSBLM");
    start = clock();
    for(int i = 0; i < nRepeats; ++i) results[i] = simulatedAnnealingSBeforeLocalMinimize(f, s, initialGuesses[i]);
    debugResultGlobalBox(results, f, matrix, start);
    DEBUG("DifferentialEvolutionBLM");
    matrix.lastItem().append("DifferentialEvolutionBLM");
    start = clock();
    for(int i = 0; i < nRepeats; ++i) results[i] = differentialEvolutionBeforeLocalMinimize(f, s, box);
    debugResultGlobalBox(results, f, matrix, start);
    DEBUG("GeneticLSBLM");
    matrix.lastItem().append("GeneticLSBLM");
    start = clock();
    for(int i = 0; i < nRepeats; ++i) results[i] = geneticLocalSearchBeforeLocalMinimize(f, s);
    debugResultGlobalBox(results, f, matrix, start);
}
template<typename TESTCASE> void testAllSolversGlobalBox(TESTCASE const& f,
    Vector<Vector<string> >& matrix, bool useBox, bool smallBudget, int nRepeats)
{
    BoxSampler sb(f.getBox());
    if(useBox)
    {//box-aware data
        BoxConstrainedStepSampler ms(f.getBox());
        testAllSolversGlobalHelper(f, sb, ms, useBox, smallBudget, nRepeats, matrix);
    }
    else
    {//agnostic data
        UnboundedSampler s(sb());//to remove unfair advantage of 0
        AgnosticStepSampler ms;
        testAllSolversGlobalHelper(f, s, ms, useBox, smallBudget, nRepeats, matrix);
    }
}

void testAllFunctionsGlobalBox(bool useBox = true, bool smallBudget = false)
{//ignore starting points given in test set
    Vector<Vector<string> > matrix;
    Vector<TestFunctionsGlobalBoxMin::MetaF> fs =
        TestFunctionsGlobalBoxMin::getFunctions();
    //int nRepeats = 1;
    //int nRepeats = 10;
    int nRepeats = 30;
    for(int i = 0; i < fs.getSize(); ++i) if(fs[i].getBox().getSize() <= 10)
    {
        string name = fs[i].getName();
        DEBUG(name);
        matrix.append(Vector<string>());
        matrix.lastItem().append(name);
        testAllSolversGlobalBox(fs[i], matrix, useBox, smallBudget, nRepeats);
    }
    string name = useBox ? "reportMinGlobalBox" : "reportMinGlobalNoBox";
    if(smallBudget) name += "smallBudget";
    createMinReport(name, matrix);
}

int main()
{
    testAllFunctionsGlobalBox(true);
    testAllFunctionsGlobalBox(false);
    return 0;
}

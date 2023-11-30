#ifndef IGMDK_RANDOM_H
#define IGMDK_RANDOM_H
#include <ctime>
#include <cassert>
#include <cmath>
#include <limits>
#include "../Utils/Utils.h"
#include "../Utils/Vector.h"
#include "../Utils/Stack.h"
#include "../Utils/Bits.h"//for rightmost 0 count
#include <cstdint>
#include "Generators.h"
using namespace std;
namespace igmdk{

template<typename GENERATOR = QualityXorshift64> struct Random
{
    GENERATOR g;
    enum{PASSWORD = 19870804};
    Random(unsigned long long seed = time(0) ^ PASSWORD): g(seed){}
    unsigned long long next(){return g.next();}
    unsigned long long maxNextValue(){return g.maxNextValue();}
    unsigned long long mod(unsigned long long n)
    {
        assert(n > 0);
        return next() % n;
    }
    int sign(){return mod(2) ? 1 : -1;}
    long long inRange(long long a, long long b)
    {
        assert(b >= a);
        return a + mod(b - a + 1);
    }
    double uniform01(){return g.uniform01();}
    bool bernoulli(double p){return uniform01() <= p;}
    int binomial(double p, int n)
    {
        int result = 0;
        for(int i = 0; i < n; ++i) result += bernoulli(p);
        return result;
    }
    int geometric(double p)
    {
        assert(p > 0);
        int result = 1;
        while(!bernoulli(p)) ++result;
        return result;
    }
    int geometric05(){return rightmost0Count(next()) + 1;}
    int poisson(double l)
    {
        assert(l > 0);
        int result = -1;
        for(double p = 1; p > exp(-l); p *= uniform01()) ++result;
        return result;
    }
    double uniform(double a, double b){return a + (b - a) * uniform01();}
    double exponential(double a){return -log(uniform01())/a;}
    double cauchy(double m, double q)
        {return (tan((uniform01() - 0.5) * PI()) + m) * q;}
    double weibull1(double b){return pow(exponential(1), 1/b);}
    double normal01()
    {
        for(;;)
        {
            double a = 2 * uniform01() - 1, b = 2 * uniform01() - 1,
                c = a * a + b * b;
            if(c < 1)
            {
                double temp = sqrt(-2 * log(c)/c);
                return a * temp;//can return b * temp as 2nd iid sample
            }
        }
    }
    double normal(double m, double stdev){return m + stdev * normal01();}
    double logNormal(double m, double q){return exp(normal(m, q));}
    double gamma1(double b)
    {
        assert(b > 0);
        if(b >= 1)
        {
            for(double third = 1.0/3, d = b - third, x, v, u, xs;;)
            {
                do
                {
                    x = normal01();
                    v = 1 + x * third/sqrt(d);
                }while(v <= 0);
                v *= v * v; u = uniform01(), xs = x * x;
                if(u > 0.0331 * xs * xs || log(u) < xs/2 +
                    d * (1 - v + log(v))) return d * v;
            }
        }
        else return pow(uniform01(), 1/b) * gamma1(b + 1);
    }
    double erlang(double m, int k){return gamma1(k) * m/k;}
    double chiSquared(int k){return 2 * gamma1(k/2.0);}
    double t(int v){return sqrt(v/chiSquared(v)) * normal01();}
    double beta(double p, double q)
    {
        double G1 = gamma1(p);
        return G1/(G1 + gamma1(q));
    }
    double F(int v1, int v2)
        {return v2 * chiSquared(v1)/(v1 * chiSquared(v2));}
    double triangular01(double middle)
    {
        assert(0 < middle && middle < 1);
        double u = uniform01();
        return sqrt(u <= middle ? middle * u : (1 - middle) * (1 - u));
    }
    double triangular(double a, double b, double middle)
        {return a + (b - a) * triangular01((middle - a)/(b - a));}
    double Levy(double scale = 0.455)
    {
        double temp = normal(0, 1/sqrt(scale));
        return 1/(temp * temp);
    }
    double symmetrizedLevy(double scale = 0.455){return sign() * Levy(scale);}
    template<typename ITEM> void randomPermutation(ITEM* numbers, int size)
    {
        for(int i = 0; i < size; ++i)
            swap(numbers[i], numbers[inRange(i, size - 1)]);
    }
    Vector<int> sortedSample(int k, int n)
    {
        assert(k <= n && k > 0);
        Vector<int> result;
        for(int considered = 0, selected = 0; selected < k; ++considered)
            if(bernoulli(double(k - selected)/(n - considered)))
            {//select
                result.append(considered);
                ++selected;
            }
        return result;
    }
    Vector<int> randomCombination(int k, int n)
    {
        Vector<int> result = sortedSample(k, n);
        randomPermutation(result.getArray(), k);
        return result;
    }
    double uniformOrderStatistic(int i, int n){return beta(i, n - i + 1);}
    Vector<double> randomUnitVector(int n)
    {
        Vector<double> result(n);
        for(int i = 0; i < n; ++i) result[i] = uniform01() * sign();
        result *= 1/norm(result);
        return result;
    }
    pair<double, double> pointInUnitCircle()
    {
        double x = uniform(-1, 1), y = uniform(-1, 1);
        while(x * x + y * y > 1)
        {//regenerate if repeated
            x = uniform(-1, 1);
            y = uniform(-1, 1);
        }
        return make_pair(x, y);
    }
};
Random<>& GlobalRNG()
{
    static Random<> r;//runs only once
    return r;
}

template<typename ITEM> void permuteDeterministically(ITEM* a, int size)
{
    Random<> drng(0);
    drng.randomPermutation(a, size);
}

class AliasMethod
{
    int n;
    Vector<int> aliases;//-1 means no donor
    Vector<double> wealth;//original or after donation
public:
    AliasMethod(Vector<double> const& probabilities):
        n(probabilities.getSize()), aliases(n, -1), wealth(n, 0)
    {
        Stack<int> smaller, greater;
        for(int i = 0; i < n; ++i)
        {//separate into poor and rich
            (wealth[i] = n * probabilities[i]) < 1 ?
                smaller.push(i) : greater.push(i);
        }
        while(!smaller.isEmpty() && !greater.isEmpty())
        {//reassign wealth until no poor remain
            int rich = greater.getTop(), poor = smaller.pop();
            aliases[poor] = rich;
            wealth[rich] -= 1 - wealth[poor];
            if(wealth[rich] < 1) smaller.push(greater.pop());
        }
    }
    int next()
    {//-1 check handles wealth round-off accumulation
        int x = GlobalRNG().mod(n);
        return aliases[x] == -1 || GlobalRNG().uniform01() < wealth[x] ?
            x : aliases[x];
    }
};

template<typename ITEM> class SumHeap
{
    Vector<ITEM> heap;
    int parent(int i){return (i - 1)/2;}
    int leftChild(int i){return 2 * i + 1;}
    int rightChild(int i){return 2 * i + 2;}
public:
    ITEM total(){return heap[0];}

    ITEM get(int i)
    {
        ITEM result = heap[i];
        int c = leftChild(i);
        if(c < heap.getSize())
        {
            result -= heap[c];
            c = rightChild(i);
            if(c < heap.getSize()) result -= heap[c];
        }
        return result;
    }

    void add(ITEM value, int i = -1)
    {//-1 means no nodes yet
        if(i == -1)
        {
            i = heap.getSize();
            heap.append(0);
        }
        for(;; i = parent(i))
        {
            heap[i] += value;
            if(i < 1) break;
        }
    }

    int find(ITEM value)
    {
        assert(0 <= value && value <= total());
        ITEM totalLeftValue = 0;
        for(int i = 0, c;; i = c)
        {
            c = leftChild(i);
            if(c >= heap.getSize()) return i;//leaf node
            if(value > totalLeftValue + heap[c])
            {
                totalLeftValue += heap[c];
                c = rightChild(i);
                if(c >= heap.getSize() ||//check if value in parent
                    value > totalLeftValue + heap[c]) return i;
            }
        }
    }
    int next(){return find(GlobalRNG().uniform01()*total());}

    ITEM cumulative(int i)
    {
        ITEM sum = heap[i];
        while(i > 0)
        {//add value of every left sibling if
            int last = i;
            i = parent(i);
            int l = leftChild(i);
            if(l != last) sum += heap[l];
        }
        return sum;
    }
};

template<typename ITEM> struct ReservoirSampler
{
    int k, nProcessed;
    Vector<ITEM> selected;
    void processItem(ITEM const& item)
    {
        ++nProcessed;
        if(selected.getSize() < k) append(item);//select first k
        else
        {//replace random
            int kickedOut = GlobalRNG().mod(nProcessed);
            if(kickedOut < k) selected[kickedOut] = item;
        }
    }
    ReservoirSampler(int wantedSize): k(wantedSize), nProcessed(0){}
};

template<typename BIASED_COIN> class FairCoin
{
    BIASED_COIN b;
    FairCoin(BIASED_COIN const& theB = BIASED_COIN()): b(theB){}
    bool flip()
    {//HT is 1, TH is 0, ignore HH and TT
        bool flip1;
        do
        {
            flip1 = b.flip();
        }while(flip1 == b.flip());
        return flip1;
    }
};

void normalizeProbs(Vector<double>& probs)
{
    double sum = 0;
    for(int i = 0; i < probs.getSize(); ++i) sum += probs[i];
    probs *= 1/sum;
}

struct NormalSummary: public ArithmeticType<NormalSummary>
{
    double mean, variance;
    double stddev()const{return sqrt(variance);}
    double error95()const{return 2 * stddev();}
    explicit NormalSummary(double theMean = 0, double theVariance = 0):
        mean(theMean), variance(theVariance){assert(variance >= 0);}
    NormalSummary operator+=(NormalSummary const& b)
    {//for sum add means and variances
        mean += b.mean;
        variance += b.variance;
        return *this;
    }
    NormalSummary operator-=(NormalSummary const& b)
    {//for difference subtract means but add variances
        mean -= b.mean;
        variance += b.variance;
        return *this;
    }
    NormalSummary operator*=(double a)
    {//scale both mean and variance
        mean *= a;
        variance *= a * a;
        return *this;
    }
};

struct IncrementalStatistics
{
    double sum, squaredSum, minimum, maximum;
    long long n;
    IncrementalStatistics(): n(0), sum(0), squaredSum(0),
        minimum(numeric_limits<double>::infinity()), maximum(-minimum){}
    double getMean()const{return sum/n;}
    double getVariance()const{return n < 2 ? 0 :
        max(0.0, (squaredSum - sum * getMean())/(n - 1.0));}
    double stdev()const{return sqrt(getVariance());}
    void addValue(double x)
    {//update incremental variables
        ++n;
        maximum = max(maximum, x);
        minimum = min(minimum, x);
        sum += x;
        squaredSum += x * x;
    }
    NormalSummary getStandardErrorSummary()const
        {return NormalSummary(getMean(), getVariance()/n);}
};

template<typename FUNCTION> IncrementalStatistics MonteCarloSimulate(
    FUNCTION const& f, long long nSimulations = 10000)
{
    IncrementalStatistics s;
    for(long long i = 0; i < nSimulations; ++i) s.addValue(f());
    return s;
}

}//end namespace
#endif

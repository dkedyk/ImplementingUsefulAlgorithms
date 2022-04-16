#ifndef IGMDK_METAHEURISTICS_H
#define IGMDK_METAHEURISTICS_H
#include "../Utils/Utils.h"
#include "../RandomNumberGeneration/Random.h"
#include "../RandomNumberGeneration/Statistics.h"
namespace igmdk{

template<typename EVALUATOR> class PartitionProblemSwapMove
{
    EVALUATOR const& e;
public:
    typedef pair<Vector<Vector<int> >,
        typename EVALUATOR::INCREMENTAL_STATE> X;
    struct MOVE{int index, from, to;};
    PartitionProblemSwapMove(EVALUATOR const& theE): e(theE){}
    pair<MOVE, double> proposeMove(X const& x)const
    {
        assert(x.first.getSize() >= 2);//else code breaks
        Vector<int> fromTo =
            GlobalRNG().randomCombination(2, x.first.getSize() + 1);
        //to can be to last bin, from not
        if(fromTo[0] == x.first.getSize()) swap(fromTo[0], fromTo[1]);
        Vector<int> const& bin1 = x.first[fromTo[0]];
        assert(bin1.getSize() > 0);//else bad problem
        int index = GlobalRNG().mod(bin1.getSize());
        MOVE m = {index, fromTo[0], fromTo[1]};
        Vector<Vector<int> >& x2 = (Vector<Vector<int> >&)x.first;
        bool openNew = false;
        if(fromTo[1] == x2.getSize())
        {
            openNew = true;
            x2.append(Vector<int>());
        }
        pair<MOVE, double> result = make_pair(m, e.evalStep(index, m.from,
            m.to, x.first, x.second));
        if(openNew) x2.removeLast();
        return result;
    }
    void applyMove(X& x, MOVE const& m)const
    {//update incremental state
        assert(m.from >= 0 && m.from < x.first.getSize() &&
            m.to >= 0 && m.to <= x.first.getSize() && m.index >= 0 &&
            m.index < x.first[m.from].getSize());
        if(m.to == x.first.getSize()) x.first.append(Vector<int>());
        x.second = e.updateIncrementalState(m.index, m.from, m.to, x.first,
            x.second);
        //do swap
        Vector<int>& bin1 = x.first[m.from];
        x.first[m.to].append(bin1[m.index]);
        bin1[m.index] = bin1.lastItem();
        bin1.removeLast();
        //remove bin1 if it became empty
        if(bin1.getSize() == 0)
        {
            bin1 = x.first.lastItem();
            x.first.removeLast();
        }
    }
    double getScore(X const& x)const{return e(x.first);}
};

template<typename EVALUATOR> class SubsetProblemFlipMove
{
    EVALUATOR const& e;
public:
    typedef pair<Vector<bool>, typename EVALUATOR::INCREMENTAL_STATE> X;
    typedef int MOVE;
    SubsetProblemFlipMove(EVALUATOR const& theE): e(theE){}
    pair<MOVE, double> proposeMove(X const& x)const
    {
        int i = GlobalRNG().mod(x.first.getSize());
        return make_pair(i, e.evalStep(i, x.first, x.second));
    }
    void applyMove(X& x, MOVE const& i)const
    {//update incremental state
        assert(i >= 0 && i < x.first.getSize());
        x.second = e.updateIncrementalState(i, x.first, x.second);
        //do flip
        x.first[i] = !x.first[i];
    }
    double getScore(X const& x)const{return e(x.first);}
};

template<typename EVALUATOR> class SymmetricPermutationProblemReverseMove
{
    EVALUATOR const& e;
public:
    typedef pair<Vector<int>, typename EVALUATOR::INCREMENTAL_STATE> X;
    typedef pair<int, int> MOVE;
    SymmetricPermutationProblemReverseMove(EVALUATOR const& theE): e(theE){}
    pair<MOVE, double> proposeMove(X const& x)const
    {
        Vector<int> m = GlobalRNG().randomCombination(2, x.first.getSize());
        if(m[0] > m[1]) swap(m[0], m[1]);
        return make_pair(make_pair(m[0], m[1]),
            e.evalReverse(m[0], m[1], x.first, x.second));
    }
    void applyMove(X& x, MOVE const& m)const
    {
        assert(m.first < m.second && m.first >= 0 &&
            m.second < x.first.getSize());
        x.second = e.updateIncrementalStateReverse(m.first, m.second, x.first,
            x.second);
        x.first.reverse(m.first, m.second);
    }
    double getScore(X const& x)const{return e(x.first);}
};
template<typename EVALUATOR> class PermutationProblemSwapMove
{
    EVALUATOR const& e;
public:
    typedef pair<Vector<int>, typename EVALUATOR::INCREMENTAL_STATE> X;
    typedef pair<int, int> MOVE;
    PermutationProblemSwapMove(EVALUATOR const& theE): e(theE){}
    pair<MOVE, double> proposeMove(X const& x)const
    {
        Vector<int> m = GlobalRNG().randomCombination(2, x.first.getSize());
        if(m[0] > m[1]) swap(m[0], m[1]);//unnecessary, but cleaner this way
        return make_pair(make_pair(m[0], m[1]),
            e.evalSwap(m[0], m[1], x.first, x.second));
    }
    void applyMove(X& x, MOVE const& m)const
    {
        assert(m.first < m.second && m.first >= 0 &&
            m.second < x.first.getSize());
        x.second = e.updateIncrementalStateSwap(m.first, m.second, x.first,
            x.second);
        swap(x.first[m.first], x.first[m.second]);
    }
    double getScore(X const& x)const{return e(x.first);}
};

template<typename PROBLEM> typename PROBLEM::X localSearch(PROBLEM const& p,
    typename PROBLEM::X x, int maxMoves, int maxStall = -1)
{//the second value of the below is the score
    typedef pair<typename PROBLEM::MOVE, double> MOVE;
    for(int i = 0; maxMoves-- && (maxStall == -1 || i < maxStall); ++i)
    {
        MOVE m = p.proposeMove(x);
        if(m.second <= 0)
        {
            i = -1;//reset counter on accept
            p.applyMove(x, m.first);
        }
    }
    return x;
}

template<typename INSTANCE> Vector<int>
    solveSymmetricPermutationLocalSearchReverse(INSTANCE const& instance,
    Vector<int> const& initial, int maxMoves)
{//use 10nlg(n) as local search stop limit
    int nPossibleMoves = initial.getSize() * initial.getSize()/2;
    return localSearch(SymmetricPermutationProblemReverseMove<INSTANCE>(
        instance), make_pair(initial, instance.getIncrementalState(initial)),
        maxMoves, int(10 * nPossibleMoves * log(nPossibleMoves))).first;
}

template<typename INSTANCE> Vector<int> solvePermutationLocalSearchSwap(
    INSTANCE const& instance, Vector<int> const& initial, int maxMoves)
{//use 10nlg(n) as local search stop limit
    int nPossibleMoves = initial.getSize() * initial.getSize()/2;
    return localSearch(PermutationProblemSwapMove<INSTANCE>(
        instance), make_pair(initial, instance.getIncrementalState(initial)),
        maxMoves, int(10 * nPossibleMoves * log(nPossibleMoves))).first;
}

template<typename INSTANCE> Vector<bool> solveSubsetLocalSearchFlip(
    INSTANCE const& instance, Vector<bool> const& initial, int maxMoves)
{//use 10nlg(n) as local search stop limit
    return localSearch(SubsetProblemFlipMove<INSTANCE>(instance),
        make_pair(initial, instance.getIncrementalState(initial)), maxMoves,
        int(10 * initial.getSize() * log(initial.getSize()))).first;
}

template<typename INSTANCE> Vector<Vector<int> >
    solvePartitionLocalSearchSwap(INSTANCE const& instance,
    Vector<Vector<int> > const& initial, int maxMoves)
{//use 10nlg(n) as local search stop limit
    int nPossibleMoves = initial.getSize() * initial.getSize()/2;
    return localSearch(PartitionProblemSwapMove<INSTANCE>(instance),
        make_pair(initial, instance.getIncrementalState(initial)), maxMoves,
        int(10 * nPossibleMoves * log(nPossibleMoves))).first;
}

template<typename PROBLEM> typename PROBLEM::X simulatedAnnealing(
    PROBLEM const& p, typename PROBLEM::X x, double T, double coolingFactor,
    double TCap, int maxMoves)
{
    assert(maxMoves > 0 && isfinite(T) && T > TCap && TCap >= 0 &&
        coolingFactor < 1 && coolingFactor > 0 && isfinite(p.getScore(x)));
    typedef pair<typename PROBLEM::MOVE, double> MOVE;
    for(; maxMoves--; T *= coolingFactor)
    {//comparison false if m.second = NaN
        MOVE m = p.proposeMove(x);
        if(m.second <= 0 || max(m.second, TCap) < T *
           GlobalRNG().exponential(1)) p.applyMove(x, m.first);
    }
    return x;
}
template<typename PROBLEM> typename PROBLEM::X
    selfTunedSimulatedAnnealing(PROBLEM const& p, typename PROBLEM::X x,
    int maxMoves = 1000000)
{
    assert(maxMoves >= 9);//min for estimators to work
    int nEstimate = int(sqrt(maxMoves));//reasonable for a good estimate
    typedef pair<typename PROBLEM::MOVE, double> MOVE;
    //T estimation stage
    Vector<double> values;
    for(int i = 0; i < nEstimate; ++i)
    {//make some random moves to get magnitude of changes
        MOVE m = p.proposeMove(x);
        double change = abs(m.second);
        if(isfinite(change) && change > 0)
        {
            values.append(change);
            p.applyMove(x, m.first);
        }
    }
    maxMoves -= nEstimate;
    if(values.getSize() < 3)//bad problem, use local search
        return localSearch(p, x, maxMoves, -1);
    double prFirst = 0.9, prLast = 1.0/maxMoves,
        deltaMax = quantile(values, 0.95), deltaMin = quantile(values, 0.05),
        T0 = -deltaMax/log(prFirst), TLast = -deltaMin/log(prLast),
        inverseRange = max(TLast/T0, numeric_limits<double>::min()),
        coolingFactor = pow(inverseRange, 1.0/maxMoves);
    return simulatedAnnealing(p, x, T0, coolingFactor, deltaMin, maxMoves);
}

template<typename INSTANCE> Vector<int>
    solveSymmetricPermutationSimulatedAnnealingReverse(INSTANCE const&instance,
    Vector<int> const& initial, int maxMoves)
{
    return selfTunedSimulatedAnnealing(SymmetricPermutationProblemReverseMove<
        INSTANCE>(instance), make_pair(initial,instance.getIncrementalState(
        initial)), maxMoves).first;
}
template<typename INSTANCE> Vector<int>
    solvePermutationSimulatedAnnealingSwap(INSTANCE const& instance,
    Vector<int> const& initial, int maxMoves)
{
    return selfTunedSimulatedAnnealing(PermutationProblemSwapMove<
        INSTANCE>(instance), make_pair(initial,instance.getIncrementalState(
        initial)), maxMoves).first;
}

template<typename INSTANCE> Vector<bool> solveSubsetSimulatedAnnealing(
    INSTANCE const& instance, Vector<bool> const& initial, int maxMoves)
{
    return selfTunedSimulatedAnnealing(SubsetProblemFlipMove<INSTANCE>(
        instance), make_pair(initial, instance.getIncrementalState(initial)),
        maxMoves).first;
}

template<typename INSTANCE> Vector<Vector<int> >
    solvePartitionSimulatedAnnealing(INSTANCE const& instance,
    Vector<Vector<int> > const& initial, int maxMoves)
{
    return selfTunedSimulatedAnnealing(PartitionProblemSwapMove<INSTANCE>(
        instance), make_pair(initial, instance.getIncrementalState(initial)),
        maxMoves).first;
}

template<typename PROBLEM> typename PROBLEM::X iteratedLocalSearch(
    PROBLEM& p, typename PROBLEM::X x, int maxBigMoves)
{
    typename PROBLEM::X bestX = x;
    double bestScore = p.getScore(bestX);
    while(maxBigMoves--)
    {
        x = p.localSearchBest(x);
        //update best
        double xScore = p.getScore(x);
        if(xScore < bestScore)
        {
            bestScore = xScore;
            bestX = x;
        }
        p.bigMove(x);
    }
    return bestX;
}

template<typename EVALUATOR> struct SymmetricPermutationILSFromRandReverseMove
{
    EVALUATOR const& e;
    typedef Vector<int> X;
    int lsMoves;
    SymmetricPermutationILSFromRandReverseMove(EVALUATOR const& theE,
        int lsMaxMoves): lsMoves(lsMaxMoves), e(theE) {}
    X localSearchBest(X const& x)
        {return solveSymmetricPermutationLocalSearchReverse(e, x, lsMoves);}
    double getScore(X const& x){return e(x);}
    void bigMove(X& x)
        {GlobalRNG().randomPermutation(x.getArray(), x.getSize());}
};
template<typename INSTANCE>
Vector<int> solveSymmetricPermutationIteratedLocalSearch(INSTANCE const&
    instance, Vector<int> const& initial, int lsMaxMoves, int bigMoves)
{
    SymmetricPermutationILSFromRandReverseMove<INSTANCE> move(instance,
        lsMaxMoves/bigMoves);
    return iteratedLocalSearch(move, initial, bigMoves);
}

template<typename EVALUATOR> struct SubsetILSFromRandFlipMove
{
    EVALUATOR const& e;
    typedef Vector<bool> X;
    int lsMoves;
    SubsetILSFromRandFlipMove(EVALUATOR const& theE, int lsMaxMoves):
        lsMoves(lsMaxMoves), e(theE) {}
    X localSearchBest(X const& x)
        {return solveSubsetLocalSearchFlip(e, x, lsMoves);}
    double getScore(X const& x){return e(x);}
    void bigMove(X& x)//simple restart
        {for(int i = 0; i < x.getSize(); ++i) x[i] = GlobalRNG().mod(2);}
};
template<typename INSTANCE> Vector<bool> solveSubsetIteratedLocalSearch(
    INSTANCE const& instance, Vector<bool> const& initial, int lsMaxMoves,
    int bigMoves)
{
    SubsetILSFromRandFlipMove<INSTANCE> move(instance, lsMaxMoves/bigMoves);
    return iteratedLocalSearch(move, initial, bigMoves);
}

template<typename PROBLEM> pair<typename PROBLEM::X, double>
    geneticLocalSearch(PROBLEM const& p, int populationSize, int nLocalMoves,
    int maxEvals)
{
    assert(maxEvals > populationSize && populationSize > 1);
    int n = 2 * (populationSize/2) + 1;//must be odd, first is best
    //and the rest evolve in pairs
    Vector<pair<double, typename PROBLEM::X> > population(n),
        populationNew(n);
    for(int i = 0; i < n; ++i) population[i].first =
        p.evaluate(population[i].second = p.generate());
    PairFirstComparator<double, typename PROBLEM::X> c;
    while((maxEvals -= (n - 1) * nLocalMoves) >= 0)
    {//elitism - keep best
        populationNew[0] = population[argMin(population.getArray(), n, c)];
        //tournament selection
        for(int i = 1; i < n; ++i)
        {
            int j = GlobalRNG().mod(n), k = GlobalRNG().mod(n), winner =
                c(population[j], population[k]) ? j : k;
            populationNew[i] = population[winner];
        }
        for(int i = 1; i + 1 < n; i += 2)
        {//crossover picked parents to create children
            p.crossover(populationNew[i].second, populationNew[i + 1].second);
            //local search children and evaluate
            for(int j = 0; j < 2; ++j)
            {
                populationNew[i + j].second =
                    p.localSearch(populationNew[i + j].second, nLocalMoves);
                populationNew[i + j].first =
                    p.evaluate(populationNew[i + j].second);
            }
        }
        population = populationNew;
    }
    int best = argMin(population.getArray(), n, c);
    return make_pair(population[best].second, population[best].first);
}

template<typename FUNCTION> class GLSSubset
{
    FUNCTION const& f;
    int n;
public:
    typedef Vector<bool> X;
    GLSSubset(FUNCTION const& theF, int theN): f(theF), n(theN) {}
    X generate()const
    {
        Vector<bool> subset(n);
        for(int i = 0; i < n; ++i) subset[i] = GlobalRNG().mod(2);
        return subset;
    }
    void crossover(X& x1, X& x2)const
    {//uniform crossover
        assert(x1.getSize() == x2.getSize());
        for(int k = 0; k < x1.getSize(); ++k) if(GlobalRNG().mod(2))
            swap(x1[k], x2[k]);
    }
    X localSearch(X const& x, int nLocalMoves)const
        {return solveSubsetLocalSearchFlip(f, x, nLocalMoves);}
    double evaluate(X const& x)const{return f(x);}
};
template<typename FUNCTION> Vector<bool>
    geneticLocalSearchSubset(FUNCTION const& f, int n,
    int maxEvals = 10000000)
{
    assert(maxEvals >= n);//basic sanity check, really need 4 * nLocalMoves
    //based on incremental cost to get same generation and local with large n
    int nLocalMoves = max(n/5, int(pow(maxEvals, 1.0/3))),
        populationSize = int(sqrt(maxEvals * 1.0/nLocalMoves));
    //one generation + evaluation equivalent
    return geneticLocalSearch(GLSSubset<FUNCTION>(f, n), populationSize,
        nLocalMoves, maxEvals).first;
}

template<typename FUNCTION> class GLSPermutation
{
    FUNCTION const& f;
    int n;
    bool isSymmetric;
public:
    typedef Vector<int> X;
    GLSPermutation(FUNCTION const& theF, int theN, bool theIsSymmetric):
        f(theF), n(theN), isSymmetric(theIsSymmetric)
        {assert(theN >= 4);}//need 4 for crossover to work
    X generate()const
    {
        X p(n);
        for(int i = 0; i < n; ++i) p[i] = i;
        GlobalRNG().randomPermutation(p.getArray(), n);
        return p;
    }
    void crossover(X& x1, X& x2)const
    {//order crossover of length in [2, n - 2] from random start rotation
        assert(x1.getSize() == n && x2.getSize() == n);
        int start = GlobalRNG().mod(n), length = 1 + GlobalRNG().mod(n - 3);
        Vector<bool> x1FromX2(n, false), x2FromX1(n, false);
        for(int i = 0; i < n; ++i)
        {//x1 second part from x2; x2 first part from x1
            int j = (start + i) % n;//implicit rotation start from start
            if(i < length) x2FromX1[x2[j]] = true;
            else x1FromX2[x1[j]] = true;
        }
        X x1Copy = x1;//need temp due to overlap
        for(int i = 0, i1 = (start + length) % n; i < n; ++i)
        {//fill x1 second part from x2 in its order
            int element = x2[(start + i) % n];//iterate over x1 from start
            if(x1FromX2[element])
            {
                x1[i1] = element;
                i1 = (i1 + 1) % n;//contiguous advance
            }
        }
        for(int i = 0, i2 = start; i < n; ++i)//over x2 from start
        {//fill x2 second part from original x1 in its order
            int element = x1Copy[(start + i) % n];
            if(x2FromX1[element])
            {
                x2[i2] = element;
                i2 = (i2 + 1) % n;//contiguous advance
            }
        }
    }
    X localSearch(X const& x, int nLocalMoves)const
    {//for symmetric prefer reverse move
        return isSymmetric ? solveSymmetricPermutationLocalSearchReverse(f, x,
            nLocalMoves) : solvePermutationLocalSearchSwap(f, x, nLocalMoves);
    }
    double evaluate(X const& x)const{return f(x);}
};
template<typename FUNCTION> Vector<int>
    geneticLocalSearchPermutation(FUNCTION const& f, int n, bool isSymmetric,
    int maxEvals = 10000000)
{
    assert(maxEvals >= n);//basic sanity check, really need 4 * nLocalMoves
    //based on incremental cost to get same generation and local with large n
    int nLocalMoves = max(n/5, int(pow(maxEvals, 1.0/3))),
        populationSize = int(sqrt(maxEvals * 1.0/nLocalMoves));
    //one generation + evaluation equivalent
    return geneticLocalSearch(GLSPermutation<FUNCTION>(f, n, isSymmetric),
        populationSize, nLocalMoves, maxEvals).first;
}

template<typename FUNCTION> class GLSPartition
{
    FUNCTION const& f;
    int n;
    void removeEmptyBins(Vector<Vector<int> >& partition)const
    {
        for(int i = n - 1; i >= 0; --i) if(partition[i].getSize() == 0)
        {
            partition[i] = partition.lastItem();
            partition.removeLast();
        }
    }
public:
    typedef Vector<Vector<int> > X;
    GLSPartition(FUNCTION const& theF, int theN): f(theF), n(theN) {}
    X generate()const
    {//random assignment
        Vector<Vector<int> > partition(n);
        for(int i = 0; i < n; ++i) partition[GlobalRNG().mod(n)].append(i);
        removeEmptyBins(partition);
        return partition;
    }
    void crossover(X& x1, X& x2)const
    {//convert to assignments
        X* xs[2] = {&x1, &x2};//c++ doesn't allow arrays of references
        Vector<int> assignments[2] = {Vector<int>(n, -1),Vector<int>(n, -1)};
        for(int k = 0; k < 2; ++k)
        {
            X& x = *xs[k];
            for(int bin = 0; bin < x.getSize(); ++bin)
                for(int j = 0; j < x[bin].getSize(); ++j)
                {
                    int item = x[bin][j];
                    assert(item >= 0 && item < n);
                    assignments[k][item] = bin;
                }
        }//uniform crossover on assignments
        for(int item = 0; item < n; ++item)
        {//check for bad input first--every element must be assigned
            assert(assignments[0][item] != -1 && assignments[1][item] != -1);
            if(GlobalRNG().mod(2))
                swap(assignments[0][item], assignments[1][item]);
        }//convert back to bin vectors
        for(int k = 0; k < 2; ++k)
        {
            X& x = *xs[k];
            x = Vector<Vector<int> >(n);
            for(int item = 0; item < n; ++item)
                x[assignments[k][item]].append(item);
            removeEmptyBins(x);
        }
    }
    X localSearch(X const& x, int nLocalMoves)const
        {return solvePartitionLocalSearchSwap(f, x, nLocalMoves);}
    double evaluate(X const& x)const{return f(x);}
};
template<typename FUNCTION> Vector<Vector<int> >
    geneticLocalSearchPartition(FUNCTION const& f, int n,
    int maxEvals = 10000000)
{
    assert(maxEvals >= n);//basic sanity check, really need 4 * nLocalMoves
    //based on incremental cost to get same generation and local with large n
    int nLocalMoves = max(n/5, int(pow(maxEvals, 1.0/3))),
        populationSize = int(sqrt(maxEvals * 1.0/nLocalMoves));
    //one generation + evaluation equivalent
    return geneticLocalSearch(GLSPartition<FUNCTION>(f, n), populationSize,
        nLocalMoves, maxEvals).first;
}

}//end namespace
#endif

#include "OptTestAuto.h"
#include "MetaHeuristics.h"
#include "SearchAlgorithms.h"
#include "../Utils/Vector.h"
#include "../MiscAlgs/LRUCache.h"
#include "TSP.h"
#include "Knapsack.h"
#include "Satisfiability3.h"
#include "BinPacking.h"
#include "../ExternalMemoryAlgorithms/CSV.h"
#include "../NumericalOptimization/NumericalOptimization.h"
using namespace igmdk;

template<typename PROBLEM> pair<typename PROBLEM::X, double>
    distributionEstimationMinimize(PROBLEM const& p, int n, double k,
    int maxEvals, int nLocalMoves)
{
    assert(maxEvals >= n * nLocalMoves && n > 0 && k > 0 && k < n);
    typename PROBLEM::DISTRIBUTIONS ds = p.getInitialDistributions();
    //initialize best
    typename PROBLEM::X best;
    double bestScore = numeric_limits<double>::infinity();
    maxEvals -=  nLocalMoves;//for best
    //setup for loop
    Vector<typename PROBLEM::X> populationTemp(k);
    Heap<pair<double, int>, PairFirstComparator<double, int,
        ReverseComparator<double> > > bestSamples;//max heap
    for(int iteration = 1; (maxEvals -= (n - 1) * nLocalMoves) >= 0;
        ++iteration)
    {//generate top of population; 0th iteration was default initialization
        for(int i = 0; i < n; ++i)
        {
            typename PROBLEM::X x = best;
            double score = bestScore;
            if(i != 0 || !isfinite(score))
            {//generate
                x = p.generate(ds);
                //local search
                x = p.localSearch(x, nLocalMoves);
                //evaluate
                score = p.evaluate(x);
            }
            int j = -1;//default no replace
            if(i < k) j = i;//heap not full
            else if(score < bestSamples.getMin().first)
                j = bestSamples.deleteMin().second;//new solution better
            if(j != -1)
            {
                populationTemp[j] = x;
                bestSamples.insert(make_pair(score, j));
            }
        }//update best and clear heap
        int bestIndex;//eventually assigned min
        while(!bestSamples.isEmpty())
        {
            bestScore = bestSamples.getMin().first;
            bestIndex = bestSamples.deleteMin().second;
        }
        best = populationTemp[bestIndex];
        //update distributions
        p.updateDistributions(ds, populationTemp, iteration);
    }
    return make_pair(best, bestScore);
}

template<typename EVALUATOR> struct SubsetDistributionEstimationProblem
{
    EVALUATOR const& e;
    typedef Vector<double> DISTRIBUTIONS;
    typedef Vector<bool> X;
    int n;
    SubsetDistributionEstimationProblem(EVALUATOR const& theE, int theN):
        e(theE), n(theN){}
    DISTRIBUTIONS getInitialDistributions()const
        {return Vector<double>(n, 0.5);}
    X localSearch(X const& x, int nLocalMoves)const
        {return solveSubsetLocalSearchFlip(e, x, nLocalMoves);}
    void updateDistributions(DISTRIBUTIONS& ds, Vector<X> const& sample,
        int iteration)const
    {
        assert(ds.getSize() == n);
        double rate = RMRate(iteration + 1);//start with 1/2
        for(int i = 0; i < n; ++i)
        {
            int oneCount = 0;
            for(int j = 0; j < sample.getSize(); ++j)
                oneCount += sample[j][i];
            double delta = 1.0 * oneCount/sample.getSize() - ds[i];
            ds[i] += delta * rate;
        }
    }
    X generate(DISTRIBUTIONS const& ds)const
    {
        X x(n);
        for(int i = 0; i < n; ++i) x[i] = GlobalRNG().bernoulli(ds[i]);
        return x;
    }
    double evaluate(X const& x)const{return e(x);}
};

template<typename FUNCTION> Vector<bool>
    distributionEstimationMinimizeSubset(FUNCTION const& f, int n,
    int maxEvals = 10000000)
{
    assert(maxEvals >= 100 * n/5);//basic sanity check
    int nLocalMoves = max(n/5, int(pow(maxEvals, 1.0/3))),
        populationSize = max(100, int(sqrt(maxEvals * 1.0/nLocalMoves)));
    return distributionEstimationMinimize(
        SubsetDistributionEstimationProblem<FUNCTION>(f, n),
        populationSize, populationSize/5, maxEvals, nLocalMoves).first;
}

template<typename EVALUATOR>
struct SymmetricPermutationDistributionEstimationProblem
{
    EVALUATOR const& e;
    typedef Matrix<double> DISTRIBUTIONS;
    typedef Vector<int> X;
    int n;
    SymmetricPermutationDistributionEstimationProblem(EVALUATOR const& theE,
        int theN): e(theE), n(theN){}
    DISTRIBUTIONS getInitialDistributions()const
    {//symmetric, not normalized
        Matrix<double> inverseDistances(n, n);
        for(int r = 0; r < n; ++r)
            for(int c = 0; c < r; ++c) if(r != c)
                inverseDistances(c, r) = inverseDistances(r, c) = 1;
        return inverseDistances;
    }
    X localSearch(X const& x, int nLocalMoves)const
        {return solveSymmetricPermutationLocalSearchReverse(e,x, nLocalMoves);}
    void updateDistributions(DISTRIBUTIONS& ds, Vector<X> const& sample,
        int iteration)const
    {
        assert(ds.getRows() == n && ds.getColumns() == n);
        Matrix<double> counts(n, n);
        for(int j = 0; j < sample.getSize(); ++j)
            for(int i = 1; i < n; ++i)
            {
                int r = sample[j][i - 1], c = sample[j][i];
                ++counts(r, c);
                ++counts(c, r);
            }//scale to have update of same frobenius norm
        double scale = normFrobenius(ds)/normFrobenius(counts);
        ds += (counts * scale - ds) * RMRate(iteration + 1);//start with 1/2
    }
    X generate(DISTRIBUTIONS const& ds)const
    {//random start to avoid bias
        Vector<int> permutation(n);
        Vector<bool> visited(n, false);
        permutation[0] = GlobalRNG().mod(n);
        for(int i = 1; i < n; ++i)
        {
            int from = permutation[i - 1];
            visited[from] = true;
            Vector<int> toMap;
            Vector<double> probs;
            double sum = 0;
            for(int to = 0; to < n; ++to) if(from != to && !visited[to])
            {
                toMap.append(to);
                probs.append(ds(from, to));
                sum += ds(from, to);
            }
            if(toMap.getSize() > 1)
            {
                probs *= 1/sum;
                AliasMethod a(probs);
                permutation[i] = toMap[a.next()];
            }//last city only 1 choice
            else permutation[i] = toMap[0];
        }
        return permutation;
    }
    double evaluate(X const& x)const{return e(x);}
};

template<typename FUNCTION>
Vector<int> distributionEstimationMinimizeSymmetricPermutation(
    FUNCTION const& f, int n, int maxEvals = 10000000)
{
    assert(maxEvals >= n/5 * max(100, n));//basic sanity check
    //update O(n^2) parameters so have at least max(100, n) samples
    //also want n local moves
    int nLocalMoves = max(n/5, int(pow(maxEvals, 1.0/3))),
        populationSize = max(max(100, n),
        int(sqrt(maxEvals * 1.0/nLocalMoves)));
    return distributionEstimationMinimize(
        SymmetricPermutationDistributionEstimationProblem<FUNCTION>(f, n),
        populationSize, populationSize/5, maxEvals, nLocalMoves).first;
}

template<typename PROBLEM> class SimpleILSHelper
{
    PROBLEM const& p;
    int lsMoves, randomMoves;
public:
    SimpleILSHelper(PROBLEM const& theP, int lsMaxMoves): lsMoves(lsMaxMoves),
        randomMoves(sqrt(lsMoves)), p(theP) {}
    typedef typename PROBLEM::X X;
    X localSearchBest(X const& x)const{return localSearch(p, x, lsMoves);}
    double getScore(X const& x)const{return p.getScore(x);}
    void bigMove(X& x)const
    {
        typedef pair<typename PROBLEM::MOVE, double> MOVE;
        for(int i = 0; i < randomMoves; ++i)
            p.applyMove(x, p.proposeMove(x).first);
    }
};
template<typename PROBLEM> typename PROBLEM::X simpleIteratedLocalSearch(
    PROBLEM const& p, typename PROBLEM::X x, int maxMoves)
{
    int bigMoves = pow(maxMoves, 1.0/3);
    SimpleILSHelper<PROBLEM> h(p, maxMoves/bigMoves);
    return iteratedLocalSearch(h, x, bigMoves);
}

template<typename PROBLEM> Vector<int>
    solveSymmetricPermutationSILSReverse(PROBLEM const& p,
    Vector<int> const& initial, int maxMoves)
{
    return simpleIteratedLocalSearch(SymmetricPermutationProblemReverseMove<
        PROBLEM>(p), make_pair(initial, p.getIncrementalState(initial)),
        maxMoves).first;
}

template<typename PROBLEM> Vector<bool> solveSubsetSILSFlip(PROBLEM const& p,
    Vector<bool> const& initial, int maxMoves)
{
    return simpleIteratedLocalSearch(SubsetProblemFlipMove<PROBLEM>(p),
        make_pair(initial, p.getIncrementalState(initial)), maxMoves).first;
}

template<typename PROBLEM> Vector<Vector<int> >
    solvePartitionSILSSwap(PROBLEM const& p,
    Vector<Vector<int> > const& initial, int maxMoves)
{
    return simpleIteratedLocalSearch(PartitionProblemSwapMove<PROBLEM>(p),
        make_pair(initial, p.getIncrementalState(initial)), maxMoves).first;
}

void KnapsackCompete()
{
    DEBUG("KnapsackCompete");
    Vector<Vector<string> > matrix;
    Vector<string> row(2);
    int nRepeats = 1;
    int n = 1000;
    DEBUG(n);
    row[0] = "Problem Size";
    row[1] = to_string(n);
    matrix.append(row);
    IncrementalStatistics lsr, stsar, de, gls, ils, sils;

    Vector<bool> initial(n);
    for(int i = 0; i < n; ++i) initial[i] = GlobalRNG().mod(2);
    KnapsackRandomInstance ki(n);
    typedef KnapsackPenaltyProxy<KnapsackRandomInstance> P;
    P instance(ki);

    DEBUG("Initial");
    row[0] = "Initial";
    double score = instance(initial);
    DEBUG(score);
    row[1] = toStringDouble(score);
    matrix.append(row);
    DEBUG("Greedy");
    row[0] = "Greedy";
    score = instance(approximateKnapsackGreedy(ki));
    DEBUG(score);
    row[1] = toStringDouble(score);
    matrix.append(row);

    for(int i = 0; i < nRepeats; ++i)
    {
        DEBUG(i);
        lsr.addValue(instance(solveSubsetLocalSearchFlip(instance, initial, 10000000)));
        DEBUG(lsr.getMean());
        stsar.addValue(instance(solveSubsetSimulatedAnnealing(instance, initial, 10000000)));
        DEBUG(stsar.getMean());
        sils.addValue(instance(solveSubsetSILSFlip(instance, initial, 10000000)));
        DEBUG(sils.getMean());
        gls.addValue(instance(geneticLocalSearchSubset<P>(
            instance, n, 10000000)));
        DEBUG(gls.getMean());
        de.addValue(instance(distributionEstimationMinimizeSubset<P>(
            instance, n, 10000000)));
        DEBUG(de.getMean());
        ils.addValue(instance(solveSubsetIteratedLocalSearch<P>(
            instance, initial, 100000, 100)));
        DEBUG(ils.getMean());
    }
    DEBUG("LS Flips");
    row[0] = "LS Flips";
    score = lsr.getMean();
    DEBUG(score);
    row[1] = toStringDouble(score);
    matrix.append(row);
    DEBUG("STSA Flips");
    row[0] = "STSA Flips";
    score = stsar.getMean();
    DEBUG(score);
    row[1] = toStringDouble(score);
    matrix.append(row);
    DEBUG("Genetic");
    row[0] = "Genetic";
    score = gls.getMean();
    DEBUG(score);
    row[1] = toStringDouble(score);
    matrix.append(row);
    /*DEBUG("DE");
    row[0] = "DE";
    score = de.getMean();
    DEBUG(score);
    row[1] = toStringDouble(score);
    matrix.append(row);*/
    //below not random enough
    DEBUG("ILS");
    row[0] = "ILS";
    score = ils.getMean();
    DEBUG(score);
    row[1] = toStringDouble(score);
    matrix.append(row);
    DEBUG("B&B");
    row[0] = "B&B";
    score = instance(solveKnapsackBranchAndBound(ki, 1000000));
    DEBUG(score);
    row[1] = toStringDouble(score);
    matrix.append(row);

    //also include RBFS after get good bounds and maybe A* to see if safe
    createCSV(matrix, "Knapsack100RandResults.csv");
}

void Sat3Compete()
{
    DEBUG("Sat3Compete");
    Vector<Vector<string> > matrix;
    Vector<string> row(2);
    int nRepeats = 1;
    int n = 100, nc = 1000;
    DEBUG(n);
    DEBUG(nc);
    row[0] = "Problem Size";
    row[1] = to_string(n) + "|" + to_string(nc);
    matrix.append(row);
    IncrementalStatistics lsr, stsar, de, gls, ils, sils;

    Vector<bool> initial(n);
    for(int i = 0; i < n; ++i) initial[i] = GlobalRNG().mod(2);
    Satisfiability3RandomInstance s3i(n, nc);
    typedef SatisfiabilityProxy<Satisfiability3RandomInstance> P;
    P instance(s3i);

    DEBUG("Initial");
    row[0] = "Initial";
    double score = instance(initial);
    DEBUG(score);
    row[1] = toStringDouble(score);
    matrix.append(row);

    for(int i = 0; i < nRepeats; ++i)
    {
        DEBUG(i);
        lsr.addValue(instance(solveSubsetLocalSearchFlip(instance, initial, 10000000)));
        DEBUG(lsr.getMean());
        stsar.addValue(instance(solveSubsetSimulatedAnnealing(instance, initial, 10000000)));
        DEBUG(stsar.getMean());
        sils.addValue(instance(solveSubsetSILSFlip(instance, initial, 10000000)));
        DEBUG(sils.getMean());
        gls.addValue(instance(geneticLocalSearchSubset<P>(
            instance, n, 10000000)));
        DEBUG(gls.getMean());
        de.addValue(instance(distributionEstimationMinimizeSubset<P>(
            instance, n, 10000000)));
        DEBUG(de.getMean());
        ils.addValue(instance(solveSubsetIteratedLocalSearch<P>(
            instance, initial, 100000, 100)));
        DEBUG(ils.getMean());
    }
    DEBUG("LS Flips");
    row[0] = "LS Flips";
    score = lsr.getMean();
    DEBUG(score);
    row[1] = toStringDouble(score);
    matrix.append(row);
    DEBUG("STSA Flips");
    row[0] = "STSA Flips";
    score = stsar.getMean();
    DEBUG(score);
    row[1] = toStringDouble(score);
    matrix.append(row);
    DEBUG("Genetic");
    row[0] = "Genetic";
    score = gls.getMean();
    DEBUG(score);
    row[1] = toStringDouble(score);
    matrix.append(row);
    /*DEBUG("DE");
    row[0] = "DE";
    score = de.getMean();
    DEBUG(score);
    row[1] = toStringDouble(score);
    matrix.append(row);*/
    //below not random enough
    DEBUG("ILS");
    row[0] = "ILS";
    score = ils.getMean();
    DEBUG(score);
    row[1] = toStringDouble(score);
    matrix.append(row);

    //also include RBFS after get good bounds and maybe A* to see if safe
    createCSV(matrix, "Sat3100RandResults.csv");
}

void BinPackCompete()
{
    DEBUG("BinPackCompete");
    Vector<Vector<string> > matrix;
    Vector<string> row(2);
    int nRepeats = 1;
    int n = 1000;
    DEBUG(n);
    row[0] = "Problem Size";
    row[1] = to_string(n);
    matrix.append(row);
    IncrementalStatistics lsr, stsar, de, gls, ils, sils;
    //start with 1 bin per item
    Vector<Vector<int> > initial(n);
    for(int i = 0; i < n; ++i) initial[i].append(i);
    BinPackingRandomInstance bi(n);
    typedef BinPackingProxy<BinPackingRandomInstance> P;
    P instance(bi);

    DEBUG("Initial");
    row[0] = "Initial";
    double score = instance(initial);
    DEBUG(score);
    row[1] = toStringDouble(score);
    matrix.append(row);

    DEBUG("Greedy");
    row[0] = "Greedy";
    score = instance(firstFitDecreasing(bi, n));
    DEBUG(score);
    row[1] = toStringDouble(score);
    matrix.append(row);

    for(int i = 0; i < nRepeats; ++i)
    {
        DEBUG(i);
        lsr.addValue(instance(solvePartitionLocalSearchSwap(instance, initial, 10000000)));
        DEBUG(lsr.getMean());
        stsar.addValue(instance(solvePartitionSimulatedAnnealing(instance, initial, 10000000)));
        DEBUG(stsar.getMean());
        sils.addValue(instance(solvePartitionSILSSwap(instance, initial, 10000000)));
        DEBUG(sils.getMean());
        gls.addValue(instance(geneticLocalSearchPartition<P>(
            instance, n, 10000000)));
        DEBUG(gls.getMean());
    }
    DEBUG("LS Swaps");
    row[0] = "LS Swaps";
    score = lsr.getMean();
    DEBUG(score);
    row[1] = toStringDouble(score);
    matrix.append(row);
    DEBUG("STSA Swaps");
    row[0] = "STSA Swaps";
    score = stsar.getMean();
    DEBUG(score);
    row[1] = toStringDouble(score);
    matrix.append(row);
    DEBUG("Genetic");
    row[0] = "Genetic";
    score = gls.getMean();
    DEBUG(score);
    row[1] = toStringDouble(score);
    matrix.append(row);
    //below not random enough
    /*DEBUG("ILS");
    row[0] = "ILS";
    score = ils.getMean();
    DEBUG(score);
    row[1] = toStringDouble(score);
    matrix.append(row);*/

    //also include RBFS after get good bounds and maybe A* to see if safe
    createCSV(matrix, "BinPack100RandResults.csv");
}

bool isPermutationFrom0(Vector<int> const& p)
{//assumed every component 0 to n - 1, must be present
    int n = p.getSize();
    Bitset<> components(n);
    for(int i = 0; i < n; ++i)
    {
        if(p[i] < 0 || p[i] >= n) return false;//element out of range
        components.set(p[i], true);
    }
    components.flip();
    return components.isZero();
}

void TSPCompete()
{
    DEBUG("TSPCompete");
    Vector<Vector<string> > matrix;
    Vector<string> row(2);
    int nRepeats = 1;
    int n = 100;
    DEBUG(n);
    row[0] = "Problem Size";
    row[1] = to_string(n);
    matrix.append(row);
    row[0] = "Expected";
    double expected = sqrt(n/2.0);
    DEBUG("Expected");
    DEBUG(expected);
    row[1] = toStringDouble(expected);
    matrix.append(row);
    IncrementalStatistics lsr, stsar, stsar2, lss, stsas, de, ge, sils;

    Vector<int> initial = GlobalRNG().randomCombination(n, n);
    assert(isPermutationFrom0(initial));
    TSPRandomInstance instance(n);
    typedef TSPProxy<TSPRandomInstance> P;
    P instanceProxy(instance);
    Vector<int> solution;


    DEBUG("Initial");
    row[0] = "Initial";
    double score = instance(initial);
    DEBUG(score);
    row[1] = toStringDouble(score);
    matrix.append(row);
    for(int i = 0; i < nRepeats; ++i)
    {
        solution = solveSymmetricPermutationLocalSearchReverse(instanceProxy, initial, 10000000);
        assert(isPermutationFrom0(solution));
        lsr.addValue(instance(solution));
        DEBUG(lsr.getMean());
        solution = solveSymmetricPermutationSimulatedAnnealingReverse(instanceProxy, initial, 10000000);
        assert(isPermutationFrom0(solution));
        stsar.addValue(instance(solution));
        DEBUG(stsar.getMean());
        solution = solveSymmetricPermutationSILSReverse(instanceProxy, initial, 10000000);
        assert(isPermutationFrom0(solution));
        sils.addValue(instance(solution));
        DEBUG(sils.getMean());
        solution = solvePermutationLocalSearchSwap(instanceProxy, initial, 10000000);
        assert(isPermutationFrom0(solution));
        lss.addValue(instance(solution));
        DEBUG(lss.getMean());
        solution = solvePermutationSimulatedAnnealingSwap(instanceProxy, initial, 10000000);
        assert(isPermutationFrom0(solution));
        stsas.addValue(instance(solution));
        DEBUG(stsas.getMean());
        solution = distributionEstimationMinimizeSymmetricPermutation(instanceProxy, n, 10000000);
        assert(isPermutationFrom0(solution));
        de.addValue(instance(solution));
        DEBUG(de.getMean());
        solution = geneticLocalSearchPermutation(instanceProxy, n, true, 10000000);
        assert(isPermutationFrom0(solution));
        ge.addValue(instance(solution));
        DEBUG(ge.getMean());
    }

    DEBUG("LS Reversals");
    row[0] = "LS Reversals";
    score = lsr.getMean();
    DEBUG(score);
    row[1] = toStringDouble(score);
    matrix.append(row);
    DEBUG("STSA Reversals");
    row[0] = "STSA Reversals";
    score = stsar.getMean();
    DEBUG(score);
    row[1] = toStringDouble(score);
    matrix.append(row);
    DEBUG("GE Reversals");
    row[0] = "GE Reversals";
    score = ge.getMean();
    DEBUG(score);
    row[1] = toStringDouble(score);
    matrix.append(row);
    DEBUG("ILS Reversals");
    row[0] = "ILS Reversals";
    score = instance(solveSymmetricPermutationIteratedLocalSearch(
        instanceProxy, initial, 10000000, 100));
    DEBUG(score);
    row[1] = toStringDouble(score);
    matrix.append(row);

    DEBUG("LS Swaps");
    row[0] = "LS Swaps";
    score = lss.getMean();
    DEBUG(score);
    row[1] = toStringDouble(score);
    matrix.append(row);
    DEBUG("STSA Swaps");
    row[0] = "STSA Swaps";
    score = stsas.getMean();
    DEBUG(score);
    row[1] = toStringDouble(score);
    matrix.append(row);

    //below not random enough
    DEBUG("B&B");
    row[0] = "B&B";
    pair<Vector<int>, bool> resultBB = solveTSPBranchAndBound(instance, 10000000);
    DEBUG(resultBB.second);
    resultBB.first.debug();
    assert(isPermutationFrom0(resultBB.first));
    score = instance(resultBB.first);
    DEBUG(score);
    row[1] = toStringDouble(score);
    matrix.append(row);
    DEBUG("RTAS");
    row[0] = "RTAS";
    score = instance(solveTSPRTAS(instance));
    DEBUG(score);
    row[1] = toStringDouble(score);
    matrix.append(row);
    /*DEBUG("AStar");
    pair<Vector<int>, bool> result = solveTSPAStar(instance, 1000000);
    DEBUG(result.second);
    if(result.second)
    {
        DEBUG("A* solution");
        result.first.debug();
        assert(isPermutationFrom0(result.first));
        score = instance(result.first);
        DEBUG(score);
    }
    DEBUG("RBFS");
    result = solveTSPRBFS(instance, 10000000);
    DEBUG(result.second);
    if(result.second)
    {
        DEBUG("RBFS solution");
        result.first.debug();
        assert(isPermutationFrom0(result.first));
        score = instance(result.first);
        DEBUG(score);
    }*/
    createCSV(matrix, "TSP100RandResults.csv");
}


template<typename PROBLEM> struct AStartTSPProblemCheap
{
    PROBLEM const& problem;
    typedef int STATE_ID;//state is current permutation
    static STATE_ID nullState(){return -1;}
    typedef EHash<BUHash> HASHER;
    mutable Vector<int> stateMap;//all states ever created
    //don't know best first node for no-return case
    STATE_ID start()const{return nullState();}
    AStartTSPProblemCheap(PROBLEM const& theProblem): problem(theProblem){}
    Vector<int> convertStatePath(Vector<STATE_ID> const& path)const
    {//skip start state as irrelevant here
        Vector<int> result;
        for(int i = 0; i < path.getSize(); ++i)
            if(path[i] != start()) result.append(stateMap[path[i]]);
        return result;
    }
    Vector<int> findRemainder(Vector<int> const& path)const
    {
        int n = problem.points.getSize();
        Vector<bool> included(n, false);
        for(int i = 0; i < path.getSize(); ++i) included[path[i]] = true;
        Vector<int> result;
        for(int i = 0; i < n; ++i) if(!included[i]) result.append(i);
        return result;
    }
    template<typename VISITOR>
    bool isGoal(STATE_ID const& id, VISITOR const& v)const
    {
        Vector<int> path = convertStatePath(v.getPath(id));
        return path.getSize() == problem.points.getSize();
    }
    template<typename VISITOR>
    Vector<STATE_ID> nextStates(STATE_ID const& id, VISITOR const& v)const
    {
        Vector<STATE_ID> result;
        Vector<int> remainder =
            findRemainder(convertStatePath(v.getPath(id)));
        for(int i = 0; i < remainder.getSize(); ++i)
        {
            STATE_ID to = stateMap.getSize();
            stateMap.append(remainder[i]);
            result.append(to);
        }
        return result;
    }
    template<typename VISITOR> double remainderLowerBound(STATE_ID const& id,
        STATE_ID const& to, VISITOR const& v)const
    {
        Vector<int> path = convertStatePath(v.getPath(id));
        if(to != start())
            path.append(stateMap[to]);//for remainder order doesn't matter
        return findMSTCost(problem, findRemainder(path));
    }//return cost of direct link of the elements
    template<typename VISITOR>
    double distance(STATE_ID const& j, STATE_ID const& to, VISITOR const& dummy)const
    {
        return j == start() ? 0 :
            problem.evalStep(stateMap[j], stateMap[to]);
    }
};
template<typename INSTANCE> pair<Vector<int>, bool> solveTSPAStarCheap(INSTANCE
    const& instance, int maxSetSize)
{
    typedef AStartTSPProblemCheap<INSTANCE> P;
    P p(instance);
    pair<Vector<typename P::STATE_ID>, bool> result =
        AStar<P>::solve(p, maxSetSize);
    return make_pair(p.convertStatePath(result.first), result.second);
}

void testAStarTSP()
{
    for(int n = 5; n <= 50; n += 5)
    {
        DEBUG(n);
        double expected = sqrt(n/2.0);
        DEBUG(expected);
        Vector<int> initial = GlobalRNG().randomCombination(n, n);
        TSPRandomInstance instance(n);
        DEBUG("random");
        double score = instance(initial);
        DEBUG(score);
        pair<Vector<int>, bool> result = solveTSPAStarCheap(instance, 10000000);
        DEBUG(result.second);
        DEBUG(result.first.getSize());

        if(result.second)
        {
            DEBUG("A* Cheap solution");
            result.first.debug();
            score = instance(result.first);
            DEBUG(score);
        }
        result = solveTSPAStar(instance, 10000000);
        DEBUG(result.second);
        DEBUG(result.first.getSize());

        if(result.second)
        {
            DEBUG("A* solution");
            result.first.debug();
            score = instance(result.first);
            DEBUG(score);
        }
        //system("PAUSE");
    }
}

void testRBFSTSP()
{
    for(int n = 5; n <= 30; n += 5)
    {
        DEBUG(n);
        double expected = sqrt(n/2.0);
        DEBUG(expected);
        Vector<int> initial = GlobalRNG().randomCombination(n, n);
        TSPRandomInstance instance(n);
        DEBUG("random");
        double score = instance(initial);
        DEBUG(score);
        pair<Vector<int>, bool> result = solveTSPRBFS(instance, 10000000);
        DEBUG(result.second);
        if(result.second)
        {
            DEBUG("RFBS solution");
            result.first.debug();
            score = instance(result.first);
            DEBUG(score);
        }
    }
}

int main()
{
    testAllAutoOpt();
    testAStarTSP();
    return 0;
    BinPackCompete();
    return 0;
    Sat3Compete();
    return 0;
    KnapsackCompete();
    return 0;
    TSPCompete();
    return 0;
    testRBFSTSP();
    return 0;
    return 0;
}

#ifndef IGMDK_MRMR_H
#define IGMDK_MRMR_H

#include "../MachineLearning/Classification.h"
#include "../ComputationalGeometry/KDTree.h"
#include "../ComputationalGeometry/Point.h"

namespace igmdk{

template<typename VECTOR> struct MaxDistance
{
    static double iDistanceIncremental(
        VECTOR const& lhs, VECTOR const& rhs, int i)
        {return abs(lhs[i] - rhs[i]);}
    static double distanceIncremental(VECTOR const& lhs,
        VECTOR const& rhs, double bound = numeric_limits<double>::max())
    {
        assert(lhs.getSize() == rhs.getSize());
        double sum = 0;
        for(int i = 0; i < lhs.getSize() && sum < bound; ++i)
            sum = max(sum, iDistanceIncremental(lhs, rhs, i));
        return sum;
    }
    struct Distance
    {
        double operator()(VECTOR const& lhs, VECTOR const& rhs)const
            {return distanceIncremental(lhs, rhs);}
    };
    struct DistanceIncremental
    {
        double operator()(VECTOR const& lhs, VECTOR const& rhs)const
            {return distanceIncremental(lhs, rhs);}
        double operator()(VECTOR const& lhs, VECTOR const& rhs, int i)const
            {return iDistanceIncremental(lhs, rhs, i);}
        double operator()(double bound, VECTOR const& lhs, VECTOR const& rhs)
            const{return distanceIncremental(lhs, rhs, bound);}
    };
};

double digamma(double x) {return log(x) - 1/(2 * x) - 1/(12 * x * x);}
double estimateMI(Vector<pair<double, double> > const& data, int k = 6)
{
    assert(data.getSize() > k);
    double minX = data[0].first, maxX = data[0].first, minY = data[0].second,
        maxY = data[0].second;
    for(int i = 1; i < data.getSize(); ++i)
    {
        minX = min(minX, data[i].first);
        maxX = max(maxX, data[i].first);
        minY = min(minY, data[i].second);
        maxY = max(maxY, data[i].second);
    }
    double deltaX = maxX - minX, deltaY = maxY - minY;
    if(deltaX == 0) deltaX = 1;
    if(deltaY == 0) deltaY = 1;
    KDTree<Point2, bool> index(2);

    for(int i = 0; i < data.getSize(); ++i)
    {
        index.insert(Point2((data[i].first - minX)/deltaX,
            (data[i].second - minY)/deltaY), true);
    }
    MaxDistance<Point2>::DistanceIncremental d;
    double sum = 0;
    for(int i = 0; i < data.getSize(); ++i)
    {
        Point2 p((data[i].first - minX)/deltaX,(data[i].second - minY)/deltaY);
        Vector<KDTree<Point2, bool>::NodeType*> neighbors =
            index.kNN(p, k + 1, d);
        double dist = d(p, neighbors.lastItem()->key);
        int nx = 1, ny = 1;
        for(int j = 0; j < data.getSize(); ++j) if(i != j)
        {
            if(abs(data[j].first - data[i].first) < dist) ++nx;
            if(abs(data[j].second - data[i].second) < dist) ++ny;
        }
        sum += digamma(nx) + digamma(ny);
    }
    return digamma(k) + digamma(data.getSize()) - sum/data.getSize();
}
template<typename DATA> Vector<Bitset<> > makeFSListMRMR(DATA const& data)
{
    assert(data.getSize() > 0);
    int D = getD(data);
    Vector<double> singleFScores;
    Vector<Vector<double> > doubleFScores;
    for(int i = 0; i < D; ++i)
    {
        Bitset<> mask(D);
        mask.set(i);
        Vector<pair<double, double> > subData;
        for(int k = 0; k < data.getSize(); ++k)
            subData.append(make_pair(data.getX(k, i), data.getY(k)));
        singleFScores.append(estimateMI(subData));
        Vector<double> doubleFScoresI;
        for(int j = 0; j < i; ++j)
        {
            mask.set(j);
            for(int k = 0; k < data.getSize(); ++k)
                subData[k].second = data.getX(k, j);
            doubleFScoresI.append(estimateMI(subData));
            mask.set(j, 0);
        }
        doubleFScores.append(doubleFScoresI);
    }
    Vector<Bitset<> > result;
    Bitset<> last(D);
    for(int i = 0; i < D; ++i)
    {
        int bestJ = -1;
        double bestScore;
        for(int j = 0; j < D; ++j) if(!last[j])
        {
            double score = singleFScores[j], total = 0;
            int n = 0;
            for(int k = 0; k < D; ++k) if(last[k])
            {
                total += doubleFScores[max(j, k)][min(j, k)] -singleFScores[k];
                ++n;
            }
            if(n > 0) score += total/n;
            if(bestJ == -1 || score > bestScore)
            {
                bestJ = j;
                bestScore = score;
            }
        }
        last.set(bestJ);
        result.append(last);
    }
    return result;
}

template<typename SUBSET_LEARNER> struct MRMRLearner
{
    typedef FeatureSubsetLearner<SUBSET_LEARNER> MODEL;
    MODEL model;
public:
    template<typename DATA>
    MRMRLearner(DATA const& data): model(data,
        pickBestSubset(SCVRiskFunctor<MODEL, Bitset<>,DATA>(data),
        makeFSListMRMR(data))){}
    int predict(NUMERIC_X const& x)const{return model.predict(x);}
};

}//end namespace
#endif


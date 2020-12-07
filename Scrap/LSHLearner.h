#ifndef IGMDK_LSH_LEARNER_H
#define IGMDK_LSH_LEARNER_H

#include "../MachineLearning/Classification.h"
#include "LSH.h"

namespace igmdk{

template<typename DISTANCE = EuclideanDistance<NUMERIC_X>::Distance>
struct MeanNN
{
    Vector<NUMERIC_X> means;
    DISTANCE d;
public:
    template<typename DATA>
    MeanNN(DATA const& data, DISTANCE const& theD = DISTANCE()): d(theD)
    {
        int theNClasses = findNClasses(data);
        assert(data.getSize() > 0);
        int D = getD(data);
        NUMERIC_X zero(D);
        Vector<int> counts(theNClasses);
        means = Vector<NUMERIC_X>(theNClasses, zero);
        for(int i = 0; i < data.getSize(); ++i)
        {
            int y = data.getY(i);
            means[y] += data.getX(i);
            ++counts[y];
        }
        for(int i = 0; i < theNClasses; ++i)
            if(counts[i] > 0) means[i] *= 1.0/counts[i];
    }
public:
    int predict(NUMERIC_X const& x)const
    {
        int best = -1;
        double bestD;
        for(int i = 0; i < means.getSize(); ++i)
        {
            double dist = d(means[i], x);
            if(best == -1 || dist < bestD)
            {
                best = i;
                bestD = dist;
            }
        }
        return best;
    }
};

struct LSHE2Learner
{
    typedef E2LSHHasher HASHER;
    typedef HASHER::RESULT_TYPE RESULT_TYPE;
    //can compress buckets further using byte code for buckets, or vector
    //of chars, if max count is 255
    mutable Vector<ChainingHashTable<RESULT_TYPE, Vector<int> > > buckets;
    HASHER g;
    int nClasses;
    LSHE2Learner(double w, int k, int l, int D, int theNClasses):
        g(k, l, D, w), buckets(l), nClasses(theNClasses) {}
    void insert(NUMERIC_X const& x, int label)
    {
        for(int i = 0; i < buckets.getSize(); ++i)
        {
            RESULT_TYPE hash = g(x, i);
            Vector<int>* xBucket = buckets[i].find(hash);
            if(!xBucket)
            {
                buckets[i].insert(hash, Vector<int>(nClasses));
                xBucket = buckets[i].find(hash);
            }
            ++(*xBucket)[label];
        }
    }
    int predict(NUMERIC_X const& x)const
    {
        Vector<int> votes(nClasses);
        for(int i = 0; i < buckets.getSize(); ++i)
        {
            RESULT_TYPE hash = g(x, i);
            Vector<int>* xBucket = buckets[i].find(hash);
            if(xBucket) votes += *xBucket;
        }
        for(int i = 0; i < nClasses; ++i)
        {
            //DEBUG(votes[i]);
        }
        //system("PAUSE");
        int ama = argMax(votes.getArray(), votes.getSize()),
            ami = argMin(votes.getArray(), votes.getSize());
        return ama == ami ? -1 : ama;
    }
};

struct LSHE2LearnerWrap2
{
    LSHE2Learner l;
    MeanNN<> l2;
    template<typename DATA>
    static pair<double, int> findW(DATA const& data, int l)
    {
        assert(l >= 1);
        //sample for w and then cross-val drop 2^10?
        EuclideanDistance<NUMERIC_X>::Distance d;
        MeanNN<> l3(data);
        IncrementalStatistics s2;
        for(int i = 0; i < data.getSize(); ++i)
            s2.addValue(d(data.getX(i), l3.means[data.getY(i)]));
        DEBUG(s2.getMean());
        pair<double, int> best(-1, -1);
        double bestScore;
        for(int i = 1; i <= 1024; i *= 2)
        {
            for(int k = max(1, l/4); k <= l; ++k)
            {
                pair<double, int> next(s2.getMean()/i, k);
                double score = crossValidation<LSHE2LearnerWrap2>(next, data);
                if(best.second == -1 || score > bestScore)
                {
                    bestScore = score;
                    best = next;
                    DEBUG(bestScore);
                    DEBUG(best.first);
                    DEBUG(best.second);
                }
            }
        }
        DEBUG(best.first);
        DEBUG(best.second);
        return best;
    }
    template<typename DATA>
    LSHE2LearnerWrap2(DATA const& data, pair<double, double> const& wk, int ll = 32):
        l(wk.first, wk.second, ll, getD(data), findNClasses(data)),
        l2(data)
    {
        for(int i = 0; i < data.getSize(); ++i)
            l.insert(data.getX(i), data.getY(i));
    }
    int predict(NUMERIC_X const& x)const
    {
        int result = l.predict(x);
        if(result == -1) result = l2.predict(x);
        return result;
    }
};

struct LSHE2LearnerWrap
{
    LSHE2LearnerWrap2 model;
    template<typename DATA>
    LSHE2LearnerWrap(DATA const& data): model(data, LSHE2LearnerWrap2::findW(data, 32)){}
    int predict(NUMERIC_X const& x)const{return model.predict(x);}
};

typedef ScaledLearner<NoParamsLearner<LSHE2LearnerWrap, int>, int> LSHRange01;

}//end namespace
#endif


#ifndef IGMDK_NAIVE_BAYES_H
#define IGMDK_NAIVE_BAYES_H
#include "ClassificationCommon.h"
#include "../Utils/Vector.h"
#include "../HashTable/ChainingHashTable.h"
#include <cmath>
namespace igmdk{

class NaiveBayes
{
    struct Feature
    {
        int count;
        LinearProbingHashTable<int, int> valueCounts;
        Feature(): count(0) {}
        void add(int value)
        {
            ++count;
            int* valueCount = valueCounts.find(value);
            if(valueCount) ++*valueCount;
            else valueCounts.insert(value, 1);
        }
        double prob(int value)
        {
            int* valueCount = valueCounts.find(value);
            return (valueCount ? 1 + *valueCount : 1)/(1.0 + count);
        }
    };
    typedef ChainingHashTable<int, Feature> FEATURE_COUNTS;
    typedef ChainingHashTable<int, FEATURE_COUNTS> CLASS_COUNTS;
    mutable CLASS_COUNTS counts;
public:
    typedef Vector<pair<int, int> > SPARSE_CATEGORICAL_X;
    static SPARSE_CATEGORICAL_X convertToSparse(CATEGORICAL_X const& x)
    {
        SPARSE_CATEGORICAL_X result;
        for(int i = 0 ; i < x.getSize(); ++i)
            result.append(make_pair(i, x[i]));
        return result;
    }
    void learn(SPARSE_CATEGORICAL_X const& x, int label)
    {
        for(int i = 0; i < x.getSize(); ++i)
        {
            FEATURE_COUNTS* classCounts = counts.find(label);
            if(!classCounts) classCounts = &counts.insert(label,
                FEATURE_COUNTS())->value;
            Feature* f = classCounts->find(x[i].first);
            if(!f) f = &classCounts->insert(x[i].first, Feature())->value;
            f->add(x[i].second);
        }
    }
    int predict(SPARSE_CATEGORICAL_X const& x)const
    {
        double maxLL;
        int bestClass = -1;
        for(CLASS_COUNTS::Iterator i = counts.begin(); i != counts.end();
            ++i)
        {
            double ll = 0;
            for(int j = 0; j < x.getSize(); j++)
            {
                Feature* f = i->value.find(x[j].first);
                if(f) ll += log(f->prob(x[j].second));
            }
            if(bestClass == -1 || maxLL < ll)
            {
                maxLL = ll;
                bestClass = i->key;
            }
        }
        return bestClass;
    }
};

struct NumericalBayes
{
    NaiveBayes model;
    DiscretizerEqualWidth disc;
    template<typename DATA> NumericalBayes(DATA const& data): disc(data)
    {
        for(int i = 0; i < data.getSize(); ++i) model.learn(NaiveBayes::
            convertToSparse(disc(data.getX(i))), data.getY(i));
    }
    int predict(NUMERIC_X const& x)const
        {return model.predict(NaiveBayes::convertToSparse(disc(x)));}
};

}//end namespace
#endif


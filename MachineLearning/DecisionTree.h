#ifndef IGMDK_DECISION_TREE_H
#define IGMDK_DECISION_TREE_H
#include "ClassificationCommon.h"
#include "../Utils/Utils.h"
#include "../Sorting/Sort.h"
#include "../Utils/GCFreelist.h"
#include "../Utils/Bitset.h"
#include "../RandomNumberGeneration/Statistics.h"
#include <cmath>
namespace igmdk{

struct DecisionTree
{
    struct Node
    {
        union
        {
            int feature;//for internal nodes
            int label;//for leaf nodes
        };
        double split;
        Node *left, *right;
        bool isLeaf(){return !left;}
        Node(int theFeature, double theSplit): feature(theFeature),
            split(theSplit), left(0), right(0) {}
    }*root;
    Freelist<Node> f;
    double H(double p){return p > 0 ? p * log(1/p) : 0;}
    template<typename DATA> struct Comparator
    {
        int feature;
        DATA const& data;
        double v(int i)const{return data.data.getX(i, feature);}
        bool operator()(int lhs, int rhs)const{return v(lhs) < v(rhs);}
        bool isEqual(int lhs, int rhs)const{return v(lhs) == v(rhs);}
    };
    void rDelete(Node* node)
    {
        if(node)
        {
            rDelete(node->left);
            f.remove(node->left);
            rDelete(node->right);
            f.remove(node->right);
        }
    }
    typedef pair<Node*, int> RTYPE;
    template<typename DATA> RTYPE rHelper(DATA& data, int left, int right,
        int nClasses, double pruneZ, int depth, bool rfMode)
    {
        int D = data.getX(left).getSize(), bestFeature = -1,
            n = right - left + 1;
        double bestSplit, bestRem, h = 0;
        Comparator<DATA> co = {-1, data};
        Vector<int> counts(nClasses, 0);
        for(int j = left; j <= right; ++j) ++counts[data.getY(j)];
        for(int j = 0; j < nClasses; ++j) h += H(counts[j] * 1.0/n);
        int majority = argMax(counts.getArray(), nClasses),
            nodeAccuracy = counts[majority];
        Bitset<> allowedFeatures;
        if(rfMode)
        {//sample features for random forest
            allowedFeatures = Bitset<>(D);
            allowedFeatures.setAll(0);
            Vector<int> p = GlobalRNG().sortedSample(sqrt(D), D);
            for(int j = 0; j < p.getSize(); ++j)allowedFeatures.set(p[j], 1);
        }
        if(h > 0) for(int i = 0; i < D; ++i)//find best feature and split
            if(allowedFeatures.getSize() == 0 || allowedFeatures[i])
            {
                co.feature = i;
                quickSort(data.permutation.getArray(), left, right, co);
                int nRight = n, nLeft = 0;
                Vector<int> countsLeft(nClasses, 0), countsRight = counts;
                for(int j = left; j < right; ++j)
                {//incrementally roll counts
                    int label = data.getY(j);
                    ++nLeft;
                    ++countsLeft[label];
                    --nRight;
                    --countsRight[label];
                    double fLeft = data.getX(j, i), hLeft = 0,
                        fRight = data.getX(j + 1, i), hRight = 0;
                    if(fLeft != fRight)
                    {//don't split equal values
                        for(int l = 0; l < nClasses; ++l)
                        {
                            hLeft += H(countsLeft[l] * 1.0/nLeft);
                            hRight += H(countsRight[l] * 1.0/nRight);
                        }
                        double rem = hLeft * nLeft + hRight * nRight;
                        if(bestFeature == -1 || rem < bestRem)
                        {
                            bestRem = rem;
                            bestSplit = (fLeft + fRight)/2;
                            bestFeature = i;
                        }
                    }
                }
            }
        if(depth <= 1 || h == 0 || bestFeature == -1)
            return RTYPE(new(f.allocate())Node(majority, 0), nodeAccuracy);
        //split examples into left and right
        int i = left - 1;
        for(int j = left; j <= right; ++j)
            if(data.getX(j, bestFeature) < bestSplit)
                swap(data.permutation[j], data.permutation[++i]);
        if(i < left || i > right)
            return RTYPE(new(f.allocate())Node(majority, 0), nodeAccuracy);
        Node* node = new(f.allocate())Node(bestFeature, bestSplit);
        //recursively compute children
        RTYPE lData = rHelper(data, left, i, nClasses, pruneZ, depth - 1,
            rfMode), rData = rHelper(data, i + 1, right, nClasses, pruneZ,
            depth - 1, rfMode);
        node->left = lData.first;
        node->right = rData.first;
        int treeAccuracy = lData.second + rData.second, nTreeWins =
            treeAccuracy - nodeAccuracy, nDraws = n - nTreeWins;
        //try to prune
        if(!rfMode &&
            signTestAreEqual(nDraws/2.0, nDraws/2.0 + nTreeWins, pruneZ))
        {
            rDelete(node);
            node->left = node->right = 0;
            node->label = majority;
            node->split = 0;
            treeAccuracy = nodeAccuracy;
        }
        return RTYPE(node, treeAccuracy);
    }
    Node* constructFrom(Node* node)
    {
        Node* tree = 0;
        if(node)
        {
            tree = new(f.allocate())Node(*node);
            tree->left = constructFrom(node->left);
            tree->right = constructFrom(node->right);
        }
        return tree;
    }
public:
    template<typename DATA> DecisionTree(DATA const& data, double pruneZ = 1,
        bool rfMode = false, int maxDepth = 50): root(0)
    {
        assert(data.getSize() > 0);
        int left = 0, right = data.getSize() - 1;
        PermutedData<DATA> pData(data);
        for(int i = 0; i < data.getSize(); ++i) pData.addIndex(i);
        root = rHelper(pData, left, right, findNClasses(
            data), pruneZ, maxDepth, rfMode).first;
    }
    DecisionTree(DecisionTree const& other)
        {root = constructFrom(other.root);}
    DecisionTree& operator=(DecisionTree const& rhs)
        {return genericAssign(*this, rhs);}
    int predict(NUMERIC_X const& x)const
    {
        assert(root);//check for bad data
        Node* current = root;
        while(!current->isLeaf()) current = x[current->feature] <
            current->split ? current->left : current->right;
        return current->label;
    }
};

}//end namespace
#endif


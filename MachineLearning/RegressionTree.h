#ifndef IGMDK_REGRESSION_TREE_H
#define IGMDK_REGRESSION_TREE_H
#include "../Utils/GCFreelist.h"
#include "../Utils/Bitset.h"
#include "../RandomNumberGeneration/Statistics.h"
#include <cmath>

namespace igmdk{

struct RegressionTree
{
    struct Node
    {
        union
        {
            double split;//for internal nodes
            double label;//for leaf nodes
        };
        int feature;//for internal nodes
        Node *left, *right;
        bool isLeaf(){return !left;}
        Node(int theFeature, double theSplit): feature(theFeature),
            split(theSplit), left(0), right(0) {}
    }* root;
    Freelist<Node> f;
    double SSE(double sum, double sum2, int n)const
        {return sum2 - sum * sum/n;}
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
    double classifyHelper(NUMERIC_X const& x, Node* current)const
    {
        while(!current->isLeaf()) current = x[current->feature] <
            current->split ? current->left : current->right;
        return current->label;
    }
    template<typename DATA> Node* rHelper(DATA& data, int left, int right,
        double pruneZ, int depth, bool rfMode)
    {
        int D = data.getX(left).getSize(), bestFeature = -1,
            n = right - left + 1;
        double bestSplit, bestScore, sumY = 0, sumY2 = 0;
        Comparator<DATA> co = {-1, data};
        for(int j = left; j <= right; ++j)
        {
            double y = data.getY(j);
            sumY += y;
            sumY2 += y * y;
        }
        double ave = sumY/n, sse = max(0.0, SSE(sumY, sumY2, n));
        Bitset<> allowedFeatures;
        if(rfMode)
        {//sample features for random forest
            allowedFeatures = Bitset<>(D);
            allowedFeatures.setAll(0);
            Vector<int> p = GlobalRNG().sortedSample(sqrt(D), D);
            for(int j = 0; j < p.getSize(); ++j)allowedFeatures.set(p[j], 1);
        }
        if(sse > 0) for(int i = 0; i < D; ++i)//find best feature and split
            if(allowedFeatures.getSize() == 0 || allowedFeatures[i])
            {
                co.feature = i;
                quickSort(data.permutation.getArray(), left, right, co);
                double sumYLeft = 0, sumYRight = sumY, sumY2Left = 0,
                    sumY2Right = sumY2;
                int nRight = n, nLeft = 0;
                for(int j = left; j < right; ++j)
                {//incrementally roll counts
                    int y = data.getY(j);
                    ++nLeft;
                    sumYLeft += y;
                    sumY2Left += y * y;
                    --nRight;
                    sumYRight -= y;
                    sumY2Right -= y * y;
                    double fLeft = data.getX(j, i), score =
                        SSE(sumYLeft, sumY2Left, nLeft) +
                        SSE(sumYRight, sumY2Right, nRight),
                        fRight = data.getX(j + 1, i);
                    if(fLeft != fRight && //don't split equal values
                        (bestFeature == -1 || score < bestScore))
                    {
                        bestScore = score;
                        bestSplit = (fLeft + fRight)/2;
                        bestFeature = i;
                    }
                }
            }
        if(n < 3 || depth <= 1 || sse <= 0 || bestFeature == -1)
            return new(f.allocate())Node(-1, ave);
        //split examples into left and right
        int i = left - 1;
        for(int j = left; j <= right; ++j)
            if(data.getX(j, bestFeature) < bestSplit)
                swap(data.permutation[j], data.permutation[++i]);
        if(i < left || i > right) return new(f.allocate())Node(-1, ave);
        Node* node = new(f.allocate())Node(bestFeature, bestSplit);
        //recursively compute children
        node->left = rHelper(data, left, i, pruneZ, depth - 1, rfMode);
        node->right = rHelper(data, i + 1, right, pruneZ, depth - 1, rfMode);
        //try to prune
        double nodeWins = 0, treeWins = 0;
        for(int j = left; j <= right; ++j)
        {
            double y = data.getY(j), eNode = ave - y, eTree =
                classifyHelper(data.getX(j), node) - y;
            if(eNode * eNode == eTree * eTree)
            {
                nodeWins += 0.5;
                treeWins += 0.5;
            }
            else if(eNode * eNode < eTree * eTree) ++nodeWins;
            else ++treeWins;
        }
        if(!rfMode && signTestAreEqual(nodeWins, treeWins, pruneZ))
        {
            rDelete(node);
            node->left = node->right = 0;
            node->label = ave;
            node->feature = -1;
        }
        return node;
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
    template<typename DATA> RegressionTree(DATA const& data, double pruneZ =
        0.25, int maxDepth = 50, bool rfMode = false): root(0)
    {
        assert(data.getSize() > 0);
        int left = 0, right = data.getSize() - 1;
        PermutedData<DATA> pData(data);
        for(int i = 0; i < data.getSize(); ++i) pData.addIndex(i);
        root = rHelper(pData, left, right, pruneZ, maxDepth, rfMode);
    }
    RegressionTree(RegressionTree const& other)
        {root = constructFrom(other.root);}
    RegressionTree& operator=(RegressionTree const& rhs)
        {return genericAssign(*this, rhs);}
    double predict(NUMERIC_X const& x)const
        {return root ? classifyHelper(x, root) : 0;}
};

}//end namespace
#endif


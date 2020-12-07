#ifndef IGMDK_SCRAP_CLUSTERING_H
#define IGMDK_SCRAP_CLUSTERING_H

#include "../MachineLearning/Clustering.h"
#include "../NumericalMethods/Matrix.h"
#include "../RandomNumberGeneration/Statistics.h"
#include <cmath>
namespace igmdk{

struct EMClustering
{
    static double normalLL(Vector<double> x, Vector<double> const& m,
        Cholesky<double> const& l)
    {
        x -= m;
        return -(x.getSize() * log(2 * PI()) + l.logDet() +
            dotProduct(l.solve(x), x))/2;
    }
    static double llim(){return log(numeric_limits<double>::min())/2;}
    static pair<Vector<double>, double> findILL(NUMERIC_X const& x,
        Vector<double> const& w, Vector<Vector<double> > const& m,
        Vector<Cholesky<double> > const &ls)
    {
        int k = w.getSize();
        Vector<double> temp(k);
        for(int j = 0; j < k; ++j)
            temp[j] = w[j] > 0 ? log(w[j]) + normalLL(x, m[j], ls[j]) : 0;
        double b = argMax(temp.getArray(), k);
        for(int j = 0; j < k; ++j) temp[j] -= b;
        return make_pair(temp, b);
    }
    template<typename DATA> static double findLL(DATA const& data,
        Vector<double> const& w, Vector<Vector<double> > const& m,
        Vector<Cholesky<double> > const &ls)
    {
        double ll = 0;
        for(int i = 0; i < data.getSize(); ++i)
        {
            pair<Vector<double>, double> temp =
                findILL(data.getX(i), w, m, ls);
            double kSum = temp.second;
            for(int j = 0; j < w.getSize(); ++j) kSum += exp(temp.first[j]);
            ll += log(kSum);
        }
        return ll;
    }
    template<typename DATA> ClusterResult operator()(DATA const& data, int k,
        int maxIterations = 1000)const
    {
        int n = data.getSize(), D = getD(data);
        assert(k > 0 && k <= n && n > 0);
        //initial values
        Vector<int> assignments =
            RepeatedKMeans()(data, k, maxIterations).assignments;
        Vector<Vector<double> > m = findCentroids(data, k, assignments),
            g(n, Vector<double>(k));
        Vector<double> w(k, 1.0/k);
        double pooledVar = 0;
        for(int i = 0; i < n; ++i)
        {
            Vector<double> diff = (data.getX(i) - m[assignments[i]]);
            pooledVar += dotProduct(diff, diff);
        }
        pooledVar /= (n - k) * D;
        Vector<Cholesky<double> > ls(k,
            Cholesky<double>(Matrix<double>::identity(D) * pooledVar));
        double ll = findLL(data, w, m, ls);
        while(maxIterations--)
        {//E step
            for(int i = 0; i < n; ++i)
            {
                pair<Vector<double>, double> temp =
                    findILL(data.getX(i), w, m, ls);
                for(int j = 0; j < k; ++j) g[i][j] = exp(temp.first[j]);
                normalizeProbs(g[i]);
            }
            //M step
            bool isNumericIssue = false;
            for(int j = 0; j < k; ++j)
            {//nj
                double nj = 0;
                for(int i = 0; i < n; ++i) nj += g[i][j];
                if(nj < 1){isNumericIssue = true; break;}
                //w
                w[j] = nj/n;
                //m
                for(int i = 0; i < n; ++i) m[j] += data.getX(i) * g[i][j];
                m[j] *= 1/nj;
                //average in pooled variance
                Matrix<double> covar = Matrix<double>::identity(D) *
                    pooledVar;
                for(int i = 0; i < n; ++i)
                {
                    Vector<double> xm = data.getX(i) - m[j];
                    covar += outerProduct(xm, xm) * g[i][j];
                }
                covar *= 1/(1 + nj);
                ls[j] = Cholesky<double>(covar);
                if(ls[j].failed){isNumericIssue = true; break;}
            }
            if(isNumericIssue) break;
            double newLL = findLL(data, w, m, ls);
            if(!isfinite(newLL) || !isELess(ll, newLL, 0.00001)) break;
            ll = newLL;
        }
        double BIC = -2 * ll + (k * (1 + D + D * (D + 1)/2) - 1) * log(n);
        for(int i = 0; i < n; ++i)
            assignments[i] = argMax(g[i].getArray(), k);
        return ClusterResult(assignments, BIC);
    }
};
typedef FindKClusterer<NoParamsClusterer<EMClustering> > EMBIC;
struct EMSmart
{
    KMeansGeneral km;
    EMBIC em;
    template<typename DATA> bool useKMeans(DATA const& data)const
        {return getD(data) > 100;}//for efficiency and numerics
    template<typename DATA> ClusterResult operator()(DATA const& data,
        int k)const
    {
        if(useKMeans(data)) return km(data, k);
        else return em(data, k);
    }
    template<typename DATA> ClusterResult operator()(DATA const& data)const
    {
        if(useKMeans(data)) return km(data);
        return em(data);
    }
};


template<typename DISTANCE = EuclideanDistance<NUMERIC_X>::Distance>
struct HierarchicalClustering
{
    struct Node
    {
        union{int i, sequenceNumber;};
        int size;
        Node *left, *right;
        bool isLeaf(){return !left;}
        Node(int theI): i(theI), left(0), right(0), size(1){}
        bool operator<(Node const& rhs)const//larger wins
            {return sequenceNumber > rhs.sequenceNumber;}
    } *root;
    Freelist<Node> f;
    void createAssigmentsHelper(Node* node, Vector<int>& assignments,
        int nextC)const
    {
        if(node->isLeaf()) assignments[node->i] = nextC;
        else
        {
            createAssigmentsHelper(node->left, assignments, nextC);
            createAssigmentsHelper(node->right, assignments, nextC);
        }
    }
    Vector<int> createAssigments(int k, int n)const
    {
        Vector<int> assignments(n);
        Heap<Node*, PointerComparator<Node> > h;
        h.insert(root);
        int nextC = 0;
        while(!h.isEmpty())
        {
            Node* node = h.deleteMin();
            if(node->isLeaf() || k < 2)
                createAssigmentsHelper(node, assignments, nextC++);
            else
            {
                --k;
                h.insert(node->left);
                h.insert(node->right);
            }
        }
        return assignments;
    }
    template<typename DATA> HierarchicalClustering(
        DATA const& data, DISTANCE const& d = DISTANCE()): root(0)
    {
        int n = data.getSize(), sequenceNumber = n;
        Vector<Node*> nodeIndex(n);
        for(int i = 0; i < n; ++i) nodeIndex[i] = new(f.allocate())Node(i);
        IndexedArrayHeap<double> iHeap;
        for(int i = 1; i < n; ++i) for(int j = 0; j < i; ++j)
            iHeap.insert(d(data.getX(i), data.getX(j)), i * n + j);
        while(!iHeap.isEmpty())
        {
            IndexedArrayHeap<double>::ITEM_TYPE minP = iHeap.deleteMin();
            int index = minP.first, j = index % n,
                i = index/n;//j < i by construction
            //update smaller index pair distance, put node in small index map
            Node* winnerNode = new(f.allocate())Node(sequenceNumber++);
            int sizeI = nodeIndex[i]->size, sizeJ = nodeIndex[j]->size;
            winnerNode->left = nodeIndex[j];
            winnerNode->right = nodeIndex[i];
            winnerNode->size = sizeI + sizeJ;
            nodeIndex[j] = winnerNode;
            nodeIndex[i] = 0;
            for(int k = 0; k < n; ++k)//update all i and j distances
                if(i != k && j != k && nodeIndex[k])
                {
                    int indexIK = max(i, k) * n + min(i, k),
                        indexJK = max(j, k) * n + min(j, k);
                    double dijk = (sizeI * *iHeap.find(indexIK) +
                        sizeJ * *iHeap.find(indexJK))/(sizeI + sizeJ);
                    iHeap.changeKey(dijk, indexJK);
                    iHeap.deleteKey(indexIK);
                }
        }
        root = nodeIndex[0];
    }
    struct Functor
    {
        template<typename DATA> ClusterResult operator()(DATA const& data,
            int k, HierarchicalClustering const& h) const
        {
            Vector<int> assignments = h.createAssigments(k, data.getSize());
            return ClusterResult(assignments, -clusterSilhouette(data,
                assignments, DISTANCE()));
        }
        template<typename DATA> ClusterResult operator()(DATA const& data,
            int k)const
            {return operator()(data, k, HierarchicalClustering(data));}
        template<typename DATA> ClusterResult operator()(DATA const& data)const
        {//Wrong! - needs to recompute -- need to fix!
            return findClustersAndK(data, Functor(),
                HierarchicalClustering(data));
        }
    };
};

template<typename DISTANCE = EuclideanDistance<NUMERIC_X>::Distance>
struct DBSCAN
{
    template<typename DATA, typename TREE> static pair<Vector<int>, double>
        findClusters(TREE const& t, DATA const& data, int nMin, double eps)
    {//determine core points and their border points
        int n = data.getSize();
        UnionFind uf(n);
        Vector<bool> isCore(n);
        Vector<int> assignments(n, -1);
        for(int i = 0; i < n; ++i)
        {
            Vector<typename TREE::NodeType*> epsNeighbors =
                t.distanceQuery(data.getX(i), eps);
            if(epsNeighbors.getSize() - 1 >= nMin)
            {
                isCore[i] = true;
                for(int j = 0; j < epsNeighbors.getSize(); ++j)
                {
                    int v = epsNeighbors[j]->value;
                    uf.join(i, v);
                    assignments[v] = i;//reuse array for membership
                }
            }
        }//determine core/noise clusters by whether their roots are core
        for(int i = 0; i < n; ++i) if(isCore[i]) isCore[uf.isRoot(i)] = true;
        int k = 0;
        for(int i = 0; i < n; ++i)
            if(uf.isRoot(i) && isCore[i]) assignments[i] = k++;
        //find classes of border points +
        //code -1 as k + 1 for compatibility with analysis code
        int noise = 0;
        for(int i = 0; i < n; ++i)
        {
            if(assignments[i] == -1)
            {
                assignments[i] = k;
                ++noise;
            }
            else if(!uf.isRoot(i))
                assignments[i] = assignments[uf.find(assignments[i])];
        }
        return make_pair(assignments, noise * 1.0/n);
    }
    template<typename DATA> ClusterResult operator()(
        DATA const& data, double noisePercentage = 0.05)const
    {//estimate params
        int n = data.getSize(), nMin = log(n) + 1;
        assert(n > 0);
        typedef VpTree<typename DATA::X_TYPE, int, DISTANCE> TREE;
        TREE t;
        for(int i = 0; i < n; ++i) t.insert(data.getX(i), i);
        IncrementalStatistics s;
        DISTANCE d;
        for(int i = 0; i < n; ++i)
        {
            Vector<typename TREE::NodeType*> nMinNNs = t.kNN(data.getX(i),
                nMin + 1);//+1 for self
            s.addValue(d(data.getX(i), nMinNNs.lastItem()->key));
        }
        double stdev = sqrt(s.getVariance()), z = -3;
        pair<Vector<int>, double> result;
        do
        {
            result = findClusters(t, data, nMin, s.getMean() + z * stdev);
            z += 0.5;
        }while(result.second >= noisePercentage && z <= 3);
        return ClusterResult(result.first);
    }
};

}//end namespace
#endif


#ifndef IGMDK_SPARSE_MATRIX_H
#define IGMDK_SPARSE_MATRIX_H
#include "../Utils/Vector.h"
#include "../Sorting/Sort.h"
#include "../Utils/Utils.h"
#include "../RandomTreap/Treap.h"
#include "../HashTable/LinearProbingHashTable.h"
#include "Matrix.h"
namespace igmdk{

template<typename ITEM = double>
class SparseMatrix: public ArithmeticType<SparseMatrix<ITEM> >
{
    int rows;
    typedef pair<int, ITEM> Item;
    typedef Vector<Item> SparseVector;
    Vector<SparseVector> itemColumns;
    int findPosition(int r, int c)const
    {//return index of the given element or -1 if it doesn't exist
        assert(0 <= r && r < rows && 0 <= c && c < getColumns());
        SparseVector const& column = itemColumns[c];
        return binarySearch(column.getArray(), 0, column.getSize() - 1,
            Item(r, 0), PairFirstComparator<int, ITEM>());
    }
public:
    SparseMatrix(int theRows, int theColumns): rows(theRows),
        itemColumns(theColumns){}
    int getRows()const{return rows;}
    int getColumns()const{return itemColumns.getSize();}
    ITEM operator()(int r, int c)const
    {//absent entries are 0
        int position = findPosition(r, c);
        return position == -1 ? 0 : itemColumns[c][position].second;
    }
    void set(int r, int c, ITEM const& item)
    {//first try to find and update
        SparseVector& column = itemColumns[c];
        int position = findPosition(r, c);
        if(position != -1) column[position].second = item;
        else if(item != 0)//if can't need to do vector insertion by shifting
        {//down the rest of the column
            Item temp(r, item);
            column.append(temp);
            position = column.getSize() - 1;
            for(;position > 0 && column[position - 1].first > r; --position)
                column[position] = column[position - 1];
            column[position] = temp;
        }
    }
    static SparseVector addSparseVectors(SparseVector const& a,
        SparseVector const& b)
    {//take elements from both, adding where both exist
        SparseVector result;
        for(int aj = 0, bj = 0; aj < a.getSize() || bj < b.getSize();)
        {
            bool considerBoth = aj < a.getSize() && bj < b.getSize();
            int j = aj < a.getSize() ? a[aj].first : b[bj].first;
            ITEM item = aj < a.getSize() ? a[aj++].second : b[bj++].second;
            if(considerBoth)//made init with a, now consider b
                if(j == b[bj].first) item += b[bj++].second;
                else if(j > b[bj].first)
                {
                    j = b[bj].first;
                    --aj;//undo aj increment
                    item = b[bj++].second;
                }
            if(item != 0) result.append(Item(j, item));//just in case
        }
        return result;
    }
    SparseMatrix& operator+=(SparseMatrix const& rhs)
    {//add column-by-column
        assert(rows == rhs.rows && getColumns() == rhs.getColumns());
        SparseMatrix result(rows, getColumns());
        for(int c = 0; c < getColumns(); ++c) result.itemColumns[c] =
            addSparseVectors(itemColumns[c], rhs.itemColumns[c]);
        return *this = result;
    }
    SparseMatrix& operator*=(ITEM a)
    {
        for(int c = 0; c < getColumns(); ++c)
            for(int j = 0; j < itemColumns[c].getSize(); ++j)
                itemColumns[c][j].second *= a;
        return *this;
    }
    friend SparseMatrix operator*(SparseMatrix const& A, ITEM a)
    {
        SparseMatrix result(A);
        return result *= a;
    }
    SparseMatrix operator-()const
    {
        SparseMatrix result(*this);
        return result *= -1;
    }
    SparseMatrix& operator-=(SparseMatrix const& rhs)
        {return *this += -rhs;}

    static SparseMatrix identity(int n)
    {
        SparseMatrix result(n, n);
        for(int c = 0; c < n; ++c) result.itemColumns[c].append(Item(c, 1));
        return result;
    }
    SparseMatrix transpose()const
    {
        SparseMatrix result(getColumns(), rows);
        for(int c = 0; c < getColumns(); ++c)
            for(int j = 0; j < itemColumns[c].getSize(); ++j)
                result.itemColumns[itemColumns[c][j].first].append(
                    Item(c, itemColumns[c][j].second));
        return result;
    }
    static ITEM dotSparseVectors(SparseVector const& a,
        SparseVector const& b)
    {//add to sum when both present
        ITEM result = 0;
        for(int aj = 0, bj = 0; aj < a.getSize() && bj < b.getSize();)
            if(a[aj].first == b[bj].first)
                result += a[aj++].second * b[bj++].second;
            else if(a[aj].first < b[bj].first) ++aj;
            else ++bj;
        return result;
    }
    SparseMatrix& operator*=(SparseMatrix const& rhs)
    {//O(n^2) * space factor(1 to n)
        assert(getColumns() == rhs.rows);
        SparseMatrix result(rows, rhs.getColumns()), bT = rhs.transpose();
        //compute sum of outer product sums
        typedef typename Key2DBuilder<>::WORD_TYPE W;
        LinearProbingHashTable<W, ITEM> outerSums;
        Key2DBuilder<> kb(max(result.rows, result.getColumns()),
            result.rows >= result.getColumns());
        for(int k = 0; k < rhs.rows; ++k)
            for(int aj = 0; aj < itemColumns[k].getSize(); ++aj)
                for(int btj = 0; btj < bT.itemColumns[k].getSize(); ++btj)
                {
                    int r = itemColumns[k][aj].first,
                        c = bT.itemColumns[k][btj].first;
                    W key = kb.to1D(r, c);
                    ITEM* rcSum = outerSums.find(key), rcValue =
                        itemColumns[k][aj].second *
                        bT.itemColumns[k][btj].second;
                    if(rcSum) *rcSum = rcValue;
                    else outerSums.insert(key, rcValue);
                }
        //convert outer sum hash table into final data structure
        for(typename LinearProbingHashTable<W, ITEM>::Iterator iter =
            outerSums.begin(); iter != outerSums.end(); ++iter)
        {//n must be the larger of r, c to make sense!
            pair<unsigned int, unsigned int> rc = kb.to2D(iter->key);
            result.itemColumns[rc.second].append(Item(rc.first, iter->value));
        }//sort each column to fix order
        for(int c = 0; c < result.getColumns(); ++c) quickSort(
            result.itemColumns[c].getArray(), 0,
            result.itemColumns[c].getSize() - 1,
            PairFirstComparator<int, ITEM>());
        return *this = result;
    }
    static Vector<ITEM> sparseToDense(SparseVector const& sv, int n)
    {//need n because don't know sparse tail
        assert(sv.getSize() == 0 || sv[sv.getSize() - 1].first < n);
        Vector<ITEM> v(n);
        for(int i = 0; i < sv.getSize(); ++i) v[sv[i].first] = sv[i].second;
        return v;
    }
    static SparseVector denseToSparse(Vector<ITEM> const& v)
    {
        SparseVector sv;
        for(int i = 0; i < v.getSize(); ++i)
            if(v[i] != 0) sv.append(Item(i, v[i]));
        return sv;
    }
    friend SparseVector operator*(SparseVector const& v, SparseMatrix const& A)
    {
        assert(v.getSize() == 0 || v.lastItem().first < A.getRows());
        SparseVector result;
        for(int c = 0; c < A.getColumns(); ++c)
        {//add one row at a time
            ITEM rc = dotSparseVectors(v, A.itemColumns[c]);
            if(rc != 0) result.append(Item(c, rc));
        }
        return result;
    }
    friend SparseVector operator*(SparseMatrix const& A, SparseVector const& v)
        {return v * A.transpose();}
    friend Vector<ITEM> operator*(Vector<ITEM> const& v, SparseMatrix const& A)
        {return sparseToDense(denseToSparse(v) * A, A.rows);}
    friend Vector<ITEM> operator*(SparseMatrix const& b, Vector<ITEM> const& v)
        {return sparseToDense(b * denseToSparse(v), b.getColumns());}
    friend double normInf(SparseMatrix const& A)
    {//first calculate transpose for better iteration
        SparseMatrix AT = A.transpose();
        double maxRowSum = 0;
        for(int r = 0; r < A.getRows(); ++r)
        {
            double rSum = 0;
            for(int cj = 0; cj < AT.itemColumns[r].getSize(); ++cj)
                rSum += abs(AT.itemColumns[r][cj].second);
            maxRowSum = max(maxRowSum, rSum);
        }
        return maxRowSum;
    }
    void debug()const
    {
        for(int i = 0; i < rows; ++i)
        {
            for(int j = 0; j < getColumns(); ++j)
            {
                cout << (*this)(i, j) << " ";
            }
            cout << endl;
        }
    }
};

template<typename MATRIX> MATRIX findJacobiPreconditioner(MATRIX A,
    bool isSPSD = true)
{
    if(!isSPSD) A = A.transpose() * A;
    int n = A.getRows();
    MATRIX diagInv(n, n);
    for(int r = 0; r < n; ++r) diagInv.set(r, r, 1/A(r, r));
    return diagInv;
}
template<typename MATRIX, typename PRECONDITIONER> pair<Vector<double>, double>
    conjugateGradientSolve(MATRIX const& A, Vector<double> b,
    PRECONDITIONER const& pInv, bool isSPSD = true,
    Vector<double> x = Vector<double>(), double eFactor = highPrecEps)
{
    MATRIX AT(1, 1);//dummy for SPSD case
    int n = A.getRows(), maxIter = n;
    if(x.getSize() == 0) x = Vector<double>(n, 0);
    assert(x.getSize() == n && b.getSize() == A.getColumns());
    if(!isSPSD)
    {
        AT = A.transpose();
        b = AT * b;
    }
    Vector<double> temp = A * x, r = b - (isSPSD ? temp : AT * temp),
        z = pInv * r, p = z;
    while(maxIter-- > 0 && norm(r) > eFactor * (1 + norm(b)))
    {
        Vector<double> ap = A * p;
        if(!isSPSD) ap = AT * ap;
        double rz = dotProduct(r, z), a = rz/dotProduct(p, ap);
        if(!isfinite(a)) break;
        x += p * a;
        r -= ap * a;
        z = pInv * r;
        p = z + p * (dotProduct(r, z)/rz);
    }
    return make_pair(x, norm(r));
}

template<typename MATRIX> BandMatrix<5> LanczosEigReduce(MATRIX const& A,
    Vector<Vector<double> >* vs = 0, int m = -1,
    Vector<double> v = Vector<double>())
{
    int n = A.getRows();
    if(m == -1) m = n;
    if(v.getSize() == 0) v = GlobalRNG().randomUnitVector(n);
    assert(A.getColumns() == n && v.getSize() == n);
    BandMatrix<5> result(m);
    double b = 1;
    Vector<double> prevV;
    for(int i = 0; i < m; ++i)
    {
        if(vs) vs->append(v);
        Vector<double> w = A * v;
        double a = dotProduct(w, v);
        result(i, i) = a;
        w -= v * a;
        if(i > 0)
        {
            w -= prevV * b;
            result(i, i - 1) = result(i - 1, i) = b;
        }
        b = norm(w);
        if(b < numeric_limits<double>::epsilon())
        {//can't continue so return what have
            BandMatrix<5> result2(i + 1);
            for(int j = 0; j <= i; ++j)
            {
                result2(j, j) = result(j, j);
                if(j > 0)
                    result2(j, j - 1) = result2(j - 1, j) = result(j, j - 1);
            }
            return result2;
        }
        prevV = v;
        v = w * (1/b);
    }
    return result;
}
template<typename MATRIX> pair<Vector<double>, Matrix<double> >
    LanczosEigenSymmetric(MATRIX const& A, int m = -1, int maxIter = 1000,
    double prec = highPrecEps)
{//rows of Q are eigenvectors
    int n = A.getRows();
    if(m == -1) m = n;
    Vector<Vector<double> > vs;
    BandMatrix<5> T = LanczosEigReduce(A, &vs, m);
    Matrix<double> QT(vs.getSize(), n);
    for(int c = 0; c < vs.getSize(); ++c)
        for(int r = 0; r < n; ++r) QT(c, r) = vs[c][r];
    Vector<double> eigve = QREigenTridiagonal(T, &QT, maxIter, prec);
    return make_pair(eigve, QT);
}

}//end namespace
#endif

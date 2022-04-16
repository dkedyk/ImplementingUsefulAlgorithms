#ifndef IGMDK_LDPC_H
#define IGMDK_LDPC_H
#include "../Utils/Bitset.h"
#include "../HashTable/LinearProbingHashTable.h"
#include "../NumericalMethods/NumericalMethods.h"
#include <cmath>

namespace igmdk{

class BooleanMatrix: public ArithmeticType<BooleanMatrix>
{
    int rows, columns;
    int index(int row, int column)const
    {
        assert(row >= 0 && row < rows && column >= 0 && column < columns);
        return row + column * rows;
    }
    Bitset<> items;
public:
    BooleanMatrix(int theRows, int theColumns): rows(theRows),
        columns(theColumns), items(theRows * theColumns)
        {assert(items.getSize() > 0);}
    int getRows()const{return rows;}
    int getColumns()const{return columns;}
    bool operator()(int row, int column)const
        {return items[index(row, column)];}
    void set(int row, int column, bool value = true)
        {items.set(index(row, column), value);}
    BooleanMatrix operator*=(bool scalar)
    {
        if(!scalar) items.setAll(false);
        return *this;
    }
    friend BooleanMatrix operator*(bool scalar, BooleanMatrix const& a)
    {
        BooleanMatrix result(a);
        return result *= scalar;
    }
    friend BooleanMatrix operator*(BooleanMatrix const& a, bool scalar)
        {return scalar * a;}
    BooleanMatrix& operator+=(BooleanMatrix const& rhs)
    {//+ and - are both xor
        assert(rows == rhs.rows && columns == rhs.columns);
        items ^= rhs.items;
        return *this;
    }
    BooleanMatrix& operator-=(BooleanMatrix const& rhs){return *this += rhs;}
    BooleanMatrix& operator*=(BooleanMatrix const& rhs)
    {//the usual row by column
        assert(columns == rhs.rows);
        BooleanMatrix result(rows, rhs.columns);
        for(int i = 0; i < rows; ++i)
            for(int j = 0; j < rhs.columns; ++j)
            {
                bool sum = false;
                for(int k = 0; k < columns; ++k)
                    sum ^= (*this)(i, k) * rhs(k, j);
                result.set(i, j, result(i, j) ^ sum);
            }
        return *this = result;
    }
    Bitset<> operator*(Bitset<> const& v)const
    {//matrix * vector
        assert(columns == v.getSize());
        Bitset<> result(rows);
        for(int i = 0; i < rows; ++i)
            for(int j = 0; j < columns; ++j)
                result.set(i, result[i] ^ ((*this)(i, j) * v[j]));
        return result;
    }//vector * matrix transposed
    friend Bitset<> operator*(Bitset<> const& v, BooleanMatrix const& m)
        {return m.transpose() * v;}
    static BooleanMatrix identity(int n)
    {
        BooleanMatrix result(n, n);
        for(int i = 0; i < n; ++i) result.set(i, i);
        return result;
    }
    BooleanMatrix transpose()const
    {
        BooleanMatrix result(columns, rows);
        for(int i = 0; i < rows; ++i)
            for(int j = 0; j < columns; ++j) result.set(j, i, (*this)(i, j));
        return result;
    }
    bool operator==(BooleanMatrix const& rhs)
    {
        if(rows != rhs.rows || columns != rhs.columns) return false;
        return items == rhs.items;
    }
    void debug()const
    {
        for(int i = 0; i < rows; ++i)
        {
            for(int j = 0; j < columns; ++j)
            {
                cout << (*this)(i, j) << " ";
            }
            cout << endl;
        }
    }
};

class LDPC
{
    BooleanMatrix a, g;//sparsity of A not exploited here
    struct H0Functor//for numerical solving for p
    {
        double hValue;
        double H(double p)const{return p > 0 ? p * log2(1/p) : 0;}
        double operator()(double p)const{return H(p) + H(1 - p) - hValue;}
        H0Functor(double theHValue): hValue(theHValue){}
    };
    double pFromCapacity(double capacity)const//solver guaranteed to succeed
        {return solveFor0(H0Functor(1 - capacity), 0, 0.5).first;}
    Bitset<> extractMessage(Bitset<> const &code)const
    {
        int k = getNewK(), n = code.getSize();
        Bitset<> message(k);
        for(int i = 0; i < k; ++i) message.set(i, code[i + n - k]);
        return message;
    }
    unsigned int uIndex(unsigned int r, unsigned int c)const
        {return r * a.getColumns() + c;}
public:
    int getNewK()const{return g.getColumns();}
    LDPC(int n, int k, int wc = 3): a(n - k, n), g(n, n - k)
    {
        int t = n - k, wr = n/(t/wc);
        assert(t % wc == 0 && n % wr == 0 && wc * n == wr * t && t < n);
        //create a
        for(int r = 0; r < t/wc; ++r)//the first section
            for(int c = 0; c < wr; ++c) a.set(r, c + r * wr);
        Vector<int> perm(n);
        for(int c = 0; c < n; ++c) perm[c] = c;
        for(int i = 1; i < wc; ++i)//other sections as permutations of the
        {//first
            GlobalRNG().randomPermutation(perm.getArray(), n);
            for(int r = 0; r < t/wc; ++r)
                for(int c = 0; c < wr; ++c)
                    a.set(r + i * t/wc, perm[c + r * wr]);
        }//create H from A
        BooleanMatrix h = a;
        int skip = 0;
        for(int r = 0; r < t; ++r)
        {//find column with 1, if not return
            int cNow = r - skip, c = cNow;
            for(; c < n; ++c) if(h(r, c)) break;
            if(c == n) ++skip;//all-0 row
            else if(c != cNow)//swap columns cNow and c
                for(int rb = 0; rb < t; ++rb)
                {
                    bool temp = h(rb, cNow);
                    h.set(rb, cNow, h(rb, c));
                    h.set(rb, c, temp);
                    //same for a
                    temp = a(rb, cNow);
                    a.set(rb, cNow, a(rb, c));
                    a.set(rb, c, temp);
                }
            for(int rb = 0; rb < t; ++rb)
                if(rb != r && h(rb, cNow))
                    for(c = cNow; c < n; ++c)
                        h.set(rb, c, h(rb, c) ^ h(r, c));
        }//remove 0 rows from H
        int tProper = t - skip, delta = 0;
        BooleanMatrix hNew(tProper, n);
        for(int r = 0; r < tProper; ++r)
        {//nonzero rows have correct identity part set
            while(!h(r + delta, r) && r < tProper) ++delta;
            for(int c = 0; c < n; ++c) hNew.set(r, c, h(r + delta, c));
        }//create g from h
        int kProper = n - tProper;
        g = BooleanMatrix(n, kProper);
        for(int r = 0; r < n; ++r)
            for(int c = 0; c < kProper; ++c)
                if(r < tProper) g.set(r, c, hNew(r, tProper + c));//x part
                else g.set(r, c, r - tProper == c);//identity part
        assert(a * g == BooleanMatrix(t, kProper));
    }
    Bitset<> encode(Bitset<> const& message)const
    {
        assert(message.getSize() == getNewK());
        return g * message;
    }
    pair<Bitset<>, bool> decode(Bitset<> const &code, int maxIter = 1000,
        double p = -1)const
    {
        int n = a.getColumns(), k = getNewK(), t = a.getRows();
        assert(code.getSize() == n && maxIter > 0);
        Bitset<> zero(k), corrected = code;
        if(a * code == zero) return make_pair(extractMessage(code), true);
        if(p == -1) p = pFromCapacity(1.0 * k/n);//find p if not given
        double const llr1 = log((1 - p)/p);//initialize l
        Vector<double> l(n);
        for(int i = 0; i < n; ++i) l[i] = llr1 * (code[i] ? 1 : -1);
        LinearProbingHashTable<unsigned int, double> nu;//initialize nu
        for(int r = 0; r < t; ++r) for(int c = 0; c < n; ++c) if(a(r, c))
            nu.insert(uIndex(r, c), 0);
        while(a * corrected != zero && maxIter-- > 0)//main loop
        {//update nu
            for(int r = 0; r < t; ++r)
            {
                double temp = 1;
                for(int c = 0; c < n; ++c) if(a(r, c))
                    temp *= tanh((*nu.find(uIndex(r, c)) - l[c])/2);
                for(int c = 0; c < n; ++c) if(a(r, c))
                {
                    double *nuv = nu.find(uIndex(r, c)), product = temp/
                        tanh((*nuv - l[c])/2), value = -2 * atanh(product);
                    //set numerical infinities to heuristic 100
                    if(!isfinite(value)) value = 100 * (product > 0 ? -1 : 1);
                    *nuv = value;
                }
            }//update l and the correction
            for(int c = 0; c < n; ++c)
            {
                l[c] = llr1 * (code[c] ? 1 : -1);
                for(int r = 0; r < t; ++r) if(a(r, c))
                    l[c] += *nu.find(uIndex(r, c));
                corrected.set(c, l[c] > 0);
            }
        }
        bool succeeded = maxIter > 0;
        return make_pair(succeeded ? extractMessage(corrected) : code,
            succeeded);
    }
};

}//end namespace
#endif

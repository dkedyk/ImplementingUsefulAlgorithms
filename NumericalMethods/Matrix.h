#ifndef IGMDK_MATRIX_H
#define IGMDK_MATRIX_H
#include "../Utils/Vector.h"
#include "../Utils/Utils.h"
#include "../RandomNumberGeneration/Random.h"
#include "NumericalCommon.h"
#include <cmath>
#include <complex>
namespace igmdk{

template<typename ITEM> struct Matrix: public ArithmeticType<Matrix<ITEM> >
{
    int rows, columns;
    int index(int row, int column)const
    {
        assert(row >= 0 && row < rows && column >= 0 && column < columns);
        return row + column * rows;
    }
    Vector<ITEM> items;
public:
    int getRows()const{return rows;}
    int getColumns()const{return columns;}
    ITEM& operator()(int row, int column){return items[index(row, column)];}
    ITEM const& operator()(int row, int column)const
        {return items[index(row, column)];}
    Matrix(int theRows, int theColumns): rows(theRows), columns(theColumns),
        items(rows * columns) {}

    Matrix operator*=(ITEM const& scalar)
    {
        items *= scalar;
        return *this;
    }
    friend Matrix operator*(ITEM const& scalar, Matrix const& a)
    {
        Matrix result(a);
        return result *= scalar;
    }
    friend Matrix operator*(Matrix const& a, ITEM const& scalar)
        {return scalar * a;}

    Matrix& operator+=(Matrix const& rhs)
    {
        assert(rows == rhs.rows && columns == rhs.columns);
        items += rhs.items;
        return *this;
    }
    Matrix& operator-=(Matrix const& rhs)
    {
        assert(rows == rhs.rows && columns == rhs.columns);
        items -= rhs.items;
        return *this;
    }

    Matrix& operator*=(Matrix const& rhs)
    {
        assert(columns == rhs.rows);
        Matrix result(rows, rhs.columns);
        for(int i = 0; i < rows; ++i)//row
            for(int j = 0; j < rhs.columns; ++j)//by column
            {
                ITEM sum(0);
                for(int k = 0; k < rhs.rows; ++k)
                    sum += (*this)(i, k) * rhs(k, j);
                result(i, j) += sum;
            }
        return *this = result;
    }
    Vector<ITEM> operator*(Vector<ITEM> const& v)const
    {
        assert(columns == v.getSize());
        Vector<ITEM> result(rows);
        for(int i = 0; i < rows; ++i)
            for(int j = 0; j < columns; ++j)
                result[i] += (*this)(i, j) * v[j];
        return result;
    }
    friend Vector<ITEM> operator*(Vector<ITEM> const& v, Matrix const& m)
        {return m.transpose() * v;}

    static Matrix identity(int n)
    {
        Matrix result(n, n);
        for(int i = 0; i < n; ++i) result(i, i) = ITEM(1);
        return result;
    }
    Matrix transpose()const
    {
        Matrix result(columns, rows);
        for(int i = 0; i < rows; ++i)
            for(int j = 0; j < columns; ++j) result(j, i) = (*this)(i, j);
        return result;
    }
    bool operator==(Matrix const& rhs)
    {
        if(rows != rhs.rows || columns != rhs.columns) return false;
        for(int i = 0; i < rows; ++i)
            for(int j = 0; j < columns; ++j)
                if((*this)(i, j) != rhs(i, j)) return false;
        return true;
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
    void assignSubmatrix(Matrix const& sub, int r1, int c1)
    {
        assert(r1 >= 0 && c1 >= 0 && r1 + sub.getRows() <= rows &&
            c1 + sub.getColumns() <= columns);
        for(int r = r1; r - r1 < sub.getRows(); ++r)
            for(int c = c1; c - c1 < sub.getColumns(); ++c)
                (*this)(r, c) = sub(r - r1, c - c1);
    }
};

template<typename ITEM, typename MATRIX> Matrix<ITEM> submatrix(
    MATRIX const& A, int r1, int r2, int c1, int c2)
{
    assert(r1 >= 0 && r2 < A.getRows() && c1 >= 0 && c2 < A.getColumns() &&
        r1 <= r2 && c1 <= c2);
    Matrix<ITEM> result(r2 - r1 + 1, c2 - c1 + 1);
    for(int r = r1; r <= r2; ++r)
        for(int c = c1; c <= c2; ++c) result(r - r1, c - c1) = A(r, c);
    return result;
}
template<typename ITEM, typename MATRIX> Matrix<ITEM> toDense(MATRIX const& A)
    {return submatrix<ITEM>(A, 0, A.getRows() - 1, 0, A.getColumns() - 1);}

double normInf(Matrix<double> const& A)
{
    int m = A.getRows();
    Vector<double> rowSums(m);
    for(int r = 0; r < m; ++r)
        for(int c = 0; c < A.getColumns(); ++c) rowSums[r] += abs(A(r, c));
    return valMax(rowSums.getArray(), m);
}

double normFrobenius(Matrix<double> const& A)
{
    double sum = 0;
    for(int r = 0; r < A.getRows(); ++r)
        for(int c = 0; c < A.getColumns(); ++c) sum += A(r, c) * A(r, c);
    return sqrt(sum);
}

template<typename ITEM> Matrix<ITEM> outerProduct(Vector<ITEM> const& a,
    Vector<ITEM> const& b)
{
    Matrix<ITEM> result(a.getSize(), b.getSize());
    for(int r = 0; r < a.getSize(); ++r)
        for(int c = 0; c < b.getSize(); ++c)
            result(r, c) = a[r] * b[c];
    return result;
}

Vector<double> backsubstitution(Matrix<double> const& U, Vector<double> b)
{//overwrite b with solution of Ux = b
    int n = b.getSize();
    assert(U.getRows() == n && U.getColumns() == n);
    for(int i = n - 1; i >= 0; --i)
    {
        for(int j = i + 1; j < n; ++j) b[i] -= b[j] * U(i, j);
        b[i] /= U(i, i);
    }
    return b;
}

template<typename ITEM = double> struct LUP
{
    Matrix<ITEM> d;
    Vector<int> permutation;
    bool isSingular;
    LUP(Matrix<ITEM> const& a): d(a), isSingular(false)
    {
        assert(d.rows = d.columns);//first create identity P
        for(int i = 0; i < d.rows; ++i) permutation.append(i);
        for(int i = 0; i < d.rows; ++i)
        {
            ITEM p = 0;
            int entering = -1;
            for(int j = i; j < d.rows; ++j)
                if(abs(d(i, j)) > p)
                {
                    p = abs(d(i, j));
                    entering = i;
                }
            if(entering == -1)
            {
                isSingular = true;
                continue;
            }
            swap(permutation[i], permutation[entering]);
            for(int j = 0; j < d.rows; ++j) swap(d(i, j), d(entering, j));
            for(int j = i + 1; j < d.rows; ++j)
            {
                d(j, i) /= d(i, i);
                for(int k = i + 1; k < d.rows; ++k)
                    d(j, k) -= d(j, i) * d(i, k);
            }
        }
    }

    double logAbsDet()const//-inf if 0
    {
        double result = 0;
        for(int i = 0; i < d.rows; ++i) result += log(abs(d(i, i)));
        return result;
    }
    int signDet()const
    {
        int sign = 1;
        for(int i = 0; i < d.rows; ++i)
        {
            if(d(i, i) < 0) sign *= -1;
            if(permutation[i] % 2 != i % 2) sign *= -1;
        }
        return sign;
    }

    Vector<ITEM> solve(Vector<ITEM> const& b)const
    {
        Vector<ITEM> y(d.rows);
        for(int i = 0; i < d.rows; ++i)
        {
            y[i] = b[permutation[i]];
            for(int j = 0; j < i; ++j) y[i] -= y[j] * d(i, j);
        }
        return backsubstitution(d, y);
    }
};

template<typename DECOMPOSITION> Matrix<double> inverse(DECOMPOSITION const& d,
    int n)
{
    Vector<double> identityRow(n, 0);
    Matrix<double> result(n, n);
    for(int i = 0; i < n; ++i)
    {
        identityRow[i] = 1;
        Vector<double> column = d.solve(identityRow);
        identityRow[i] = 0;
        for(int j = 0; j < n; ++j) result(j, i) = column[j];
    }
    return result;
}



template<typename ITEM> struct Cholesky
{
    Matrix<ITEM> l;
    bool failed;
    Cholesky(Matrix<ITEM> const& a): l(a.rows, a.columns), failed(false)
    {//a must be symmetric and positive definite
        for(int c = 0; c < l.columns; ++c)
            for(int r = c; r < l.rows; ++r)
            {
                ITEM sum = a(r, c);
                for(int k = 0; k < c; ++k) sum -= l(r, k) * l(c, k);
                if(r == c)
                {
                    if(sum <= 0){failed = true; return;}
                    l(c, c) = sqrt(sum);
                }
                else l(r, c) = sum/l(c, c);
            }
    }

    double logDet()const//-inf if 0
    {
        assert(!failed);
        double result = 0;
        for(int i = 0; i < l.rows; ++i) result += log(abs(l(i, i)));
        return 2 * result;
    }

    Vector<ITEM> solve(Vector<ITEM> b)const
    {
        int n = b.getSize();
        assert(l.getRows() == n && l.getColumns() == n && !failed);
        for(int i = 0; i < n; ++i)
        {
            for(int j = 0; j < i; ++j) b[i] -= b[j] * l(i, j);
            b[i] /= l(i, i);
        }
        return backsubstitution(l.transpose(), b);
    }
};

class MultivariateNormal
{
    Vector<double> means;
    Matrix<double> L;
    static Matrix<double> makeL(Matrix<double> const& covariance)
    {
        Cholesky<double> c(covariance);
        assert(!c.failed);//covariance might be wrong or have numerical issues
        return c.l;
    }
public:
    MultivariateNormal(Vector<double> const& theMeans, Matrix<double> const&
        covariance): means(theMeans), L(makeL(covariance)) {}
    Vector<double> next()const
    {
        Vector<double> normals;
        for(int i = 0 ; i < means.getSize(); ++i)
            normals.append(GlobalRNG().normal01());
        return means + L * normals;
    }
};

template<int D = 3, typename ITEM = double> struct BandMatrix
{//below main diagonal d < 0
    Vector<Vector<ITEM> > diagonals;
    int diagIndex(int d)const{return D/2 + d;}//”0” is the main diagonal
    int rowIndex(int r, int d)const{return r + min(0, d);}
    int diag(int r, int c)const{return c - r;}//storage index
    int getRows()const{return diagonals[diagIndex(0)].getSize();}
    int getColumns()const{return getRows();}
    void setupBands(int n)
    {//initialize the banded parts
        assert(n > 0 && D % 2 == 1);
        for(int d = -D/2; d <= D/2; ++d)
        diagonals[diagIndex(d)] = Vector<ITEM>(n - abs(d));
    }
    BandMatrix(int n): diagonals(D) {setupBands(n); }
    BandMatrix(Matrix<ITEM> const& A): diagonals(D)
    {//copy the banded part only
        int n = A.getRows();
        setupBands(n);//this must be called first
        assert(n == A.getColumns());
        for(int r = 0; r < n; ++r)
        for(int c = max(0, r - D/2); c <= min(n - 1, r + D/2); ++c)
        (*this)(r, c) = A(r, c);
    }
    ITEM operator()(int r, int c)const
    {
        int d = diag(r, c);
        assert(r >= 0 && r < getRows() && c >= 0 && c < getColumns());
        return abs(d) <= D/2 ? diagonals[diagIndex(d)][rowIndex(r, d)] : 0;
    }
    ITEM& operator()(int r, int c)
    {
        int d = diag(r, c);
        assert(r >= 0 && r < getRows() && c >= 0 && c < getColumns() &&
            abs(d) <= D/2);
        return diagonals[diagIndex(d)][rowIndex(r, d)];
    }
    static BandMatrix identity(int n)
    {
        BandMatrix result(n);
        result.diagonals[result.diagIndex(0)] = Vector<ITEM>(n, 1);
        return result;
    }
    BandMatrix& operator+=(BandMatrix const& A)
    {
        assert(getRows() == A.getRows());
        diagonals += A.diagonals;
        return *this;
    }
    BandMatrix& operator*=(ITEM const& a)
    {
        diagonals *= a;
        return *this;
    }
    friend BandMatrix operator*(ITEM const& scalar, BandMatrix const& A)
    {
        BandMatrix result(A);
        return result *= scalar;
    }
    BandMatrix operator-()const
    {
        BandMatrix result(*this);
        return result *= -1;
    }
    BandMatrix& operator-=(BandMatrix const& A){return *this += -A;}
    void assignSubmatrix(Matrix<ITEM> const& sub, int r1, int c1)
    {//assign only the tridiagonal part
        assert(r1 >= 0 && c1 >= 0 && r1 + sub.getRows() <= getRows() &&
            c1 + sub.getColumns() <= getColumns());
        for(int r = r1; r - r1 < sub.getRows(); ++r)
            for(int c = c1; c - c1 < sub.getColumns(); ++c)
                if(abs(diag(r, c)) <= D/2) (*this)(r, c) = sub(r - r1, c - c1);
    }
    void debug()const
    {
        for(int i = 0; i < getRows(); ++i)
        {
            for(int j = 0; j < getColumns(); ++j)
            {
                cout << (*this)(i, j) << " ";
            }
            cout << endl;
        }
    }
};
template<typename ITEM> using TridiagonalMatrix = BandMatrix<3, ITEM>;

Vector<double> solveTridiag(TridiagonalMatrix<double> const& A,
    Vector<double> r)
{//U = above diag, D = diag, L = below diag; overwrite r with solution
    Vector<double> U = A.diagonals[2], D = A.diagonals[1];
    Vector<double> const& L = A.diagonals[0];
    int n = r.getSize();
    assert(n > 1 && D.getSize() == n);
    //forward elimination
    U[0] /= D[0];
    r[0] /= D[0];
    for(int i = 1; i < n; ++i)
    {
        double denom = D[i] - L[i - 1] * U[i - 1];
        if(i < n - 1) U[i] /= denom;//allow 0 denom; user handles inf or NaN
        r[i] = (r[i] - L[i - 1] * r[i - 1])/denom;
    }//back substitution
    for(int i = n - 2; i >= 0; --i) r[i] -= U[i] * r[i + 1];
    return r;
}

double sign(double x){return x > 0 ? 1 : -1;}
Vector<double> HouseholderReduction(Vector<double> x)
{
    x[0] += sign(x[0]) * norm(x);
    double normX = norm(x);
    if(normX > 0) x *= 1/normX;
    return x;
}
Vector<double> HouseholderReduction2(double a, double b)
{
    Vector<double> x;
    x.append(a);
    x.append(b);
    return HouseholderReduction(x);
}
template<typename MATRIX> void HouseholderLeftMult(MATRIX& M,
    Vector<double> const& h, int r, int c = 0, int c2 = -1)
{
    if(c2 == -1) c2 = M.getColumns() - 1;
    assert(r + h.getSize() <= M.getRows());
    Matrix<double> sub = submatrix<double>(M, r, r + h.getSize() - 1, c, c2);
    sub -= outerProduct(h * 2, h * sub);
    M.assignSubmatrix(sub, r, c);
}
template<typename MATRIX> void HouseholderRightMult(MATRIX& M,
    Vector<double> const& h, int c, int r = 0, int r2 = -1)
{
    if(r2 == -1) r2 = M.getRows() - 1;
    Matrix<double> sub = submatrix<double>(M, r, r2, c, c + h.getSize() - 1);
    sub -= outerProduct(sub * h , h * 2);
    M.assignSubmatrix(sub, r, c);
}

struct QRDecomposition
{
    Matrix<double> Q, R;
    void doHessenbergQR()
    {
        for(int c = 0; c < R.getRows() - 1; ++c)
        {
            Vector<double> x = HouseholderReduction2(R(c, c), R(c + 1, c));
            HouseholderLeftMult(R, x, c, c);
            HouseholderLeftMult(Q, x, c);
        }
    }
    QRDecomposition(Matrix<double> const& A, bool isHessenberg = false):
        Q(Matrix<double>::identity(A.getRows())), R(A)
    {
        int n = R.getRows();
        assert(R.getColumns() == n);
        if(isHessenberg) doHessenbergQR();//Hessenberg already
        else//regular case
            for(int c = 0; c < n - 1; ++c)
            {
                Vector<double> x(n - c);
                for(int j = c; j < n; ++j) x[j - c] = R(j, c);
                x = HouseholderReduction(x);
                HouseholderLeftMult(R, x, c, c);
                HouseholderLeftMult(Q, x, c);
            }
    }
    void rank1Update(Vector<double> const& u, Vector<double> const& v)
    {//first calculate Q_wR
        int n = u.getSize();
        assert(v.getSize() == n && Q.getRows() == n);
        Vector<double> w = Q * u;
        for(int r = n - 1; r > 0; --r)
        {//find orthogonal transformation to eliminate w[r]
            Vector<double> x = HouseholderReduction2(w[r - 1], w[r]);
            Matrix<double> W(2, 1);//apply it to temp submatrix
            W(0, 0) = w[r - 1];
            W(1, 0) = w[r];
            HouseholderLeftMult(W, x, 0);
            w[r - 1] = W(0, 0);//copy back the result
            w[r] = W(1, 0);
            //apply to R also
            HouseholderLeftMult(R, x, r - 1, r - 1);
            HouseholderLeftMult(Q, x, r - 1);
        }
        for(int c = 0; c < n; ++c) R(0, c) += w[0] * v[c];
        doHessenbergQR();//process Hessenberg form of R
    }
    Vector<double> operator*(Vector<double> const& x)//Ax without forming A
        {return Q.transpose() * (R * x);}
    Vector<double> solve(Vector<double> const & b)const
    {//solve Rx = Qb to solve Ax = b
        assert(Q.getRows() == b.getSize());
        return backsubstitution(R, Q * b);
    }

    double logAbsDet()const//-inf if 0
    {
        double result = 0;
        for(int i = 0; i < R.rows; ++i) result += log(abs(R(i, i)));
        return result;
    }
};
Vector<double> outerProductMultLeft(Vector<double> const& u,//u column, v row
    Vector<double> const& v, Vector<double> const& x)
    {return u * dotProduct(v, x);}
Vector<double> outerProductMultRight(Vector<double> const& u,//u column, v row
    Vector<double> const& v, Vector<double> const& x)
    {return v * dotProduct(x, u);}

template<typename MATRIX> double wilkinsonShift(MATRIX const& A, int n)
{//use stable formula
    double d = (A(n - 2, n - 2) - A(n - 1, n - 1))/2,
        temp = A(n - 1, n - 2) * A(n - 1, n - 2);
    return A(n - 1, n - 1) - temp/(d + sign(d) * sqrt(d * d + temp));
}

struct SVD
{
    Vector<double> svds;
    Matrix<double> U, V;
    SVD(Matrix<double> A, int maxIter = 100, double prec = highPrecEps):
        svds(A.getColumns()), U(Matrix<double>::identity(A.getRows())),
        V(Matrix<double>::identity(A.getColumns()))
    {//swap U and V if needed
        bool needSwap = A.getColumns() > A.getRows();
        if(needSwap)
        {
            svds = Vector<double>(A.getRows());
            A = A.transpose();
            swap(U, V);
        }
        toBidiagonalForm(A);//convert to bidiagonal form
        int n = A.getColumns();
        //work with square submatrix if needed
        if(n < A.getRows()) A = submatrix<double>(A, 0, n - 1, 0, n - 1);
        SVDBidiagonal(A, maxIter, prec);
        if(needSwap)//swap back U and V if needed
        {
            Matrix<double> temp = V.transpose();
            V = U.transpose();
            U = temp;
        }
    }
    void toBidiagonalForm(Matrix<double>& A)
    {
        int m = A.getRows(), n = A.getColumns();
        assert(m >= n);
        for(int c = 0; c < n; ++c)
        {
            Vector<double> x(m - c);
            for(int j = c; j < m; ++j) x[j - c] = A(j, c);
            x = HouseholderReduction(x);
            HouseholderLeftMult(A, x, c, c);
            HouseholderLeftMult(U, x, c);
            if(c < n - 2)
            {
                x = Vector<double>(n - (c + 1));
                for(int j = c + 1; j < n; ++j) x[j - (c + 1)] = A(c, j);
                x = HouseholderReduction(x);
                HouseholderRightMult(A, x, c + 1, c);
                HouseholderRightMult(V, x, c + 1);
            }
        }
    }
    void SVDBidiagonal(Matrix<double>& A, int maxIter = 100,
        double prec = highPrecEps)
    {
        assert(A.getColumns() == A.getRows());
        int n = A.getRows();
        while(maxIter--)
        {
            int lastNZ = 0;
            for(int c = 0; c < n - 1; ++c) if(isEEqual(A(c, c), 0) ||
                abs(A(c, c + 1)) <= prec *
                    (abs(A(c, c)) + abs(A(c + 1, c + 1)))) A(c, c + 1) = 0;
                else lastNZ = c + 1;
            //deflate
            n = lastNZ + 1;
            if(n < 2) break;
            Matrix<double> T22 =
                submatrix<double>(A, n - 2, n - 1, n - 2, n - 1);
            double dm = T22(0, 0), dn = T22(1, 1), fm = T22(0, 1),
                fmm1 = (n >= 3 ? A(n - 3, n - 2) : 0);
            T22(0, 0) = dm * dm + fmm1 * fmm1;
            T22(0, 1) = dm * fm;
            T22(1, 0) = dm * fm;
            T22(1, 1) = dn * dn + fm * fm;
            double shift = wilkinsonShift(T22, 2),
                t00 = A(0, 0) * A(0, 0), t01 = A(0, 0) * A(0, 1),
                y = t00 - shift, z = t01;
            for(int c = 0; c < n - 1; ++c)
            {
                Vector<double> x = HouseholderReduction2(y, z);
                HouseholderRightMult(A, x, c);
                HouseholderRightMult(V, x, c);
                y = A(c, c);
                z = A(c + 1, c);
                x = HouseholderReduction2(y, z);
                HouseholderLeftMult(A, x, c);
                HouseholderLeftMult(U, x, c);
                if(c < n - 2)
                {
                    y = A(c, c + 1);
                    z = A(c, c + 2);
                }
            }
        }
        n = A.getRows();
        for(int d = 0; d < n; ++d)
        {//ensure non-negative svds by moving minuses into V
            if(A(d, d) < 0)
            {
                A(d, d) *= -1;
                for(int r = 0; r < n; ++r) V(r, d) *= -1;
            }
            svds[d] = A(d, d);
        }
    }
    int rank(double precFor0 = highPrecEps)const
    {
        int non0 = 0;
        for(int i = 0; i < svds.getSize(); ++i) non0 += (svds[i] >= precFor0);
        return non0;
    }
    double norm2()const{return valMax(svds.getArray(), svds.getSize());}
    double norm2Inv()const{return 1/valMin(svds.getArray(), svds.getSize());}
    double condition2()const{return norm2() * norm2Inv();}
};

pair<double, double> estimateMatrixEquationError2Norm(Matrix<double> const& A,
    Vector<double> const& b, Vector<double> const& x)
{//abs, rel
    SVD svd(A);
    double nr = norm(A * x - b);
    return make_pair(nr * svd.norm2Inv(), nr * svd.condition2()/norm(b));
}

pair<Matrix<double>, Matrix<double> > toHessenbergForm(Matrix<double> A)
{
    int n = A.getRows();
    assert(A.getColumns() == n);
    Matrix<double> Q = Matrix<double>::identity(n);
    for(int c = 0; c < n - 2; ++c)
    {
        Vector<double> x(n - (c + 1));
        for(int j = c + 1; j < n; ++j) x[j - (c + 1)] = A(j, c);
        x = HouseholderReduction(x);
        HouseholderLeftMult(A, x, c + 1, c);
        HouseholderLeftMult(Q, x, c + 1);
        HouseholderRightMult(A, x, c + 1);
    }
    return make_pair(A, Q);
}

template<typename MATRIX> Vector<double> QREigenTridiagonal(MATRIX A,
    Matrix<double>* Q = 0, int maxIter = 1000, double prec = highPrecEps)
{
    int n = A.getRows();
    assert(A.getColumns() == n);
    while(maxIter--)
    {
        int lastNZ = 0;
        for(int c = 0; c < n - 1; ++c) if(abs(A(c + 1, c)) <= prec *
            (abs(A(c, c)) + abs(A(c + 1, c + 1))))
            {
                A(c + 1, c) = 0;
                A(c, c + 1) = 0;
            }
            else lastNZ = c + 1;
        //deflate
        n = lastNZ + 1;
        if(n < 2) break;
        //shift
        double shift = wilkinsonShift(A, n);
        A -= shift * MATRIX::identity(A.getRows());
        //QR = A
        Vector<Vector<double> > reductions(n - 1);
        for(int c = 0; c < n - 1; ++c)
        {
            reductions[c] = HouseholderReduction2(A(c, c), A(c + 1, c));
            HouseholderLeftMult(A, reductions[c], c, c, min(n - 1, c + 2));
            if(Q) HouseholderLeftMult(*Q, reductions[c], c);
        }//A = RQ
        for(int c = 0; c < n - 1; ++c)
            HouseholderRightMult(A, reductions[c], c, max(0, c - 1), c + 1);
        //unshift
        A += shift * MATRIX::identity(A.getRows());
    }
    Vector<double> eig(A.getRows());
    for(int d = 0; d < A.getRows(); ++d) eig[d] = A(d, d);
    return eig;
}

pair<double, double> traceDet2x2(Matrix<double> const& A, int d)
{
    double trace = (A(d, d) + A(d + 1, d + 1)), det =
        A(d, d) * A(d + 1, d + 1) - A(d, d + 1) * A(d + 1, d);
    return make_pair(trace, det);
}
Vector<complex<double> > QREigenHessenberg(Matrix<double> A,
    int maxIter = 10000, double prec = highPrecEps)
{
    assert(A.getColumns() == A.getRows());
    int n = A.getRows(), lastDeflateIter = maxIter;
    while(maxIter--)
    {//process converged entries
        for(int c = 0; c < n - 1; ++c) if(abs(A(c + 1, c)) <= prec *
            (abs(A(c, c)) + abs(A(c + 1, c + 1)))) A(c + 1, c) = 0;
        //deflate
        while(n >= 3)
        {
            if(A(n - 1, n - 2) == 0) --n;
            else if(A(n - 2, n - 3) == 0) n -= 2;
            else break;
            lastDeflateIter = maxIter;
        }
        if(n < 3) break;
        pair<double, double> td = traceDet2x2(A, n - 2);
        if(lastDeflateIter - maxIter > 10)
        {//unlikely exceptional shift
            pair<double, double> ab = GlobalRNG().pointInUnitCircle();
            double r = abs(A(n - 1, n - 2));
            ab.first = A(n - 1, n - 1) + ab.first * r;
            ab.second *= r;
            td.first = 2 * ab.first;
            td.second = ab.first * ab.first + ab.second * ab.second;
            lastDeflateIter = maxIter;
        }//calculate Householder vector
        Vector<double> xyz(3);
        xyz[0] = A(0, 0) * A(0, 0) + A(0, 1) * A(1, 0) - td.first * A(0, 0) +
            td.second;
        xyz[1] = A(1, 0) * (A(0, 0) + A(1, 1) - td.first);
        xyz[2] = A(1, 0) * A(2, 1);
        for(int c = 0; c < n - 1; ++c)
        {//process it
            xyz = HouseholderReduction(xyz);
            HouseholderLeftMult(A, xyz, c, max(0, c - 1), n - 1);
            HouseholderRightMult(A, xyz, c, 0, min(c + 3, n - 1));
            if(c + 2 == n) break;
            xyz[0] = A(c + 1, c);//update it for next iteration
            xyz[1] = A(c + 2, c);
            if(c + 3 < n) xyz[2] = A(c + 3, c);
            else xyz.removeLast();
        }
    }//extract eigenvalues from blocks
    n = A.getRows();
    double NaN = numeric_limits<double>::quiet_NaN();
    Vector<complex<double> > eigs(n, complex<double>(NaN, NaN));
    for(int c = 0; c < n; ++c)
        if(c == 0 || A(c, c - 1) == 0)
        {//found block
            if(c + 1 >= n || A(c + 1, c) == 0)//1 x 1
                eigs[c] = complex<double>(A(c, c), 0);
            else if(c + 2 >= n || A(c + 2, c + 1) == 0)//2 x 2
            {
                pair<double, double> td = traceDet2x2(A, c);
                td.first /= 2;
                double temp = td.first * td.first - td.second,
                    temp2 = sqrt(abs(temp));
                if(temp > 0)
                {//stable formula for e = trace05 +- temp2
                    double eig1 = td.first + sign(td.first) * temp2;
                    eigs[c] = complex<double>(eig1, 0);
                    eigs[c + 1] = complex<double>(td.second/eig1, 0);
                }
                else
                {
                    eigs[c] = complex<double>(td.first, td.second);
                    eigs[c + 1] = complex<double>(td.first, -td.second);
                }
                ++c;
            }
            else c += 3;//didn't converge - not isolated
        }
    return eigs;
}

Vector<double> findEigenvector(Matrix<double> H, double eig,
    int maxIter = 10, double prec = highPrecEps)
{//inverse iteration
    int n = H.getRows();
    double HInf = normInf(H);
    H -= eig * Matrix<double>::identity(n);
    QRDecomposition qr(H, true);
    Vector<double> x = GlobalRNG().randomUnitVector(n);
    while(maxIter--)
    {
        x = qr.solve(x);
        for(int i = 0; i < n; ++i) if(!isfinite(x[i])) x[i] = 0;
        x *= 1/norm(x);
        if(normInf(H * x) <= prec * HInf) break;
    }
    return x;
}
Matrix<double> findEigenvectors(Matrix<double> const& H,
    Matrix<double> const& Q, Vector<complex<double> > const& eigs)
{
    int n = eigs.getSize();
    Matrix<double> result(n, n);
    for(int i = 0; i < n; ++i)
    {
        if(eigs[i].imag() == 0)
        {
            Vector<double> eigve = findEigenvector(H, eigs[i].real());
            for(int r = 0; r < n; ++r) result(r, i) = eigve[r];
        }
        else for(int r = 0; r < n; ++r)
            result(r, i) = numeric_limits<double>::quiet_NaN();
    }
    return Q.transpose() * result;
}
pair<Vector<complex<double> >, Matrix<double> > QREigen(Matrix<double> A,
    int maxIter = 1000, double prec = highPrecEps)
{
    pair<Matrix<double>, Matrix<double> > TQ = toHessenbergForm(A);
    Vector<complex<double> > eigs = QREigenHessenberg(TQ.first, maxIter, prec);
    return make_pair(eigs, findEigenvectors(TQ.first, TQ.second, eigs));
}

//update
bool isESymmetric(Matrix<double> const& A, double eRelAbs = defaultPrecEps)
{
    int n = A.getRows();
    if(n != A.getColumns()) return false;
    for(int r = 0; r < n; ++r)
    {
        if(!isfinite(A(r, r))) return false;
        for(int c = r + 1; c < n; ++c) if(!isfinite(A(r, c)) ||
            !isEEqual(A(r, c), A(c, r), eRelAbs)) return false;
    }
    return true;
}
pair<Vector<double>, Matrix<double> > QREigenSymmetric(Matrix<double> A,
    int maxIter = 1000, double prec = highPrecEps)
{//rows of Q are eigenvectors
    assert(isESymmetric(A));
    pair<Matrix<double>, Matrix<double> > TQ = toHessenbergForm(A);
    Vector<double> eigve = QREigenTridiagonal(TQ.first, &TQ.second, maxIter,
        prec);
    return make_pair(eigve, TQ.second);
}

}//end namespace
#endif

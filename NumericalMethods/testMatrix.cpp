#include "SparseMatrix.h"
#include "../Utils/Debug.h"
#include "NumericalMethodsTestAuto.h"
using namespace igmdk;

void testHessenbergPerm()
{//Permutation matrix
    Matrix<double> b2(3, 3);
    b2(0, 0) = 0;
    b2(0, 1) = 0;
    b2(0, 2) = 1;
    b2(1, 0) = 1;
    b2(1, 1) = 0;
    b2(1, 2) = 0;
    b2(2, 0) = 0;
    b2(2, 1) = 1;
    b2(2, 2) = 0;
    b2.debug();
    Vector<complex<double> > eig = QREigenHessenberg(b2);
    DEBUG("eig");
    eig.debug();
}

void testHessenbergComplex()
{//From Datta?
    Matrix<double> b2(3, 3);
    b2(0, 0) = 1;
    b2(0, 1) = 2;
    b2(0, 2) = 3;
    b2(1, 0) = 1;
    b2(1, 1) = 0;
    b2(1, 2) = 1;
    b2(2, 0) = 0;
    b2(2, 1) = -2;
    b2(2, 2) = 2;
    b2.debug();
    pair<Vector<complex<double> >, Matrix<double> > result = QREigen(b2);
    Vector<complex<double> > eig = result.first;
    DEBUG("eig");
    eig.debug();
    DEBUG("eigve");
    result.second.debug();
}

void testRepeatedNonsymmtric()
{
    int n = 3;
    Matrix<double> v(n, n), d = Matrix<double>::identity(n);
    d(0, 0) = 2;
    for(int i = 0; i < n; ++i)
    {
        Vector<double> vi = GlobalRNG().randomUnitVector(n);
        for(int j = 0; j < n; ++j) v(i, j) = vi[j];
    }
    LUP<double> lup(v);
    Matrix<double> b2 = v * d * inverse(lup, n);
    b2.debug();
    pair<Vector<complex<double> >, Matrix<double> > result = QREigen(b2);
    Vector<complex<double> > eig = result.first;
    DEBUG("eig");
    eig.debug();
    DEBUG("eigve");
    result.second.debug();
    Matrix<double> D(eig.getSize(), eig.getSize());
    for(int i = 0; i < eig.getSize(); ++i) D(i, i) = eig[i].real();
    QRDecomposition qr(result.second);
    DEBUG("V-1DV");
    (result.second * D * inverse(qr, n)).debug();
}

//to polish approximate Lancsoz values
pair<double, Vector<double> > polishSymEigenpair(SparseMatrix<double> H,
    double eig, Vector<double> x, int maxIter = 10, double prec = highPrecEps)
{//inverse iteration
    int n = H.getRows();
    double HInf = normInf(H);
    while(maxIter--)
    {
        DEBUG(eig);
        H -= SparseMatrix<double>::identity(n) * eig;
        x = conjugateGradientSolve(H, x, SparseMatrix<double>::identity(n)).first;
        for(int i = 0; i < n; ++i) if(!isfinite(x[i])) x[i] = 0;
        x *= 1/norm(x);
        if(normInf(H * x) <= prec * HInf) break;
        //update eigenvalue
        H += SparseMatrix<double>::identity(n) * eig;
        eig = dotProduct(x * H, x)/dotProduct(x, x);
    }
    return make_pair(eig, x);
}

void testTridiagonalSparse()
{
    //From Burden - correct
    SparseMatrix<double> b2(3, 3);
    b2.set(0, 0, 3);
    b2.set(0, 1, 1);
    b2.set(0, 2, 0);
    b2.set(1, 0, 1);
    b2.set(1, 1, 3);
    b2.set(1, 2, 1);
    b2.set(2, 0, 0);
    b2.set(2, 1, 1);
    b2.set(2, 2, 3);
    DEBUG("b2");
    b2.debug();
    Matrix<double> z2 = toDense<double>(b2);
    pair<Vector<double>, Matrix<double> > EQD = QREigenSymmetric(z2);
    DEBUG("dense eigva");
    EQD.first.debug();
    DEBUG("dense eigVe");
    EQD.second.debug();
    pair<Vector<double>, Matrix<double> > EQ = LanczosEigenSymmetric(b2, 2);
    Vector<double> eig = EQ.first;
    DEBUG("eigVa");
    eig.debug();
    DEBUG("eigVe");
    EQ.second.debug();
    for(int i = 0; i < EQ.second.rows; ++i)
    {//polish doesn't seem to work at all - probably too ill-conditioned for CG
        Vector<double> temp;
        for(int j = 0; j < EQ.second.columns; ++j)
        {
            temp.append(EQ.second(i, j));
        }
        pair<double, Vector<double> > polished = polishSymEigenpair(b2, EQ.first[i], temp);
        DEBUG(polished.first);
        polished.second.debug();
    }

    Matrix<double> D(eig.getSize(), eig.getSize());
    for(int i = 0; i < eig.getSize(); ++i) D(i, i) = eig[i];
    DEBUG("VTDV");
    (EQ.second.transpose() * D * EQ.second).debug();
}

int main()
{
    testTridiagonalSparse();
    return 0;
    matrixTestAllAuto();
    return 0;
    testHessenbergPerm();
    return 0;
    testHessenbergComplex();
    return 0;
    testRepeatedNonsymmtric();
    return 0;



}

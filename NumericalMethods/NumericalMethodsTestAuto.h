#ifndef IGMDK_NUMERICAL_METHODS_TEST_AUTO_H
#define IGMDK_NUMERICAL_METHODS_TEST_AUTO_H
#include "Matrix.h"
#include "SparseMatrix.h"
#include "NumericalMethods.h"
#include <cassert>
using namespace std;

namespace igmdk{

//presented
void testELessAuto()
{
    DEBUG("testELessAuto");
    double nan = numeric_limits<double>::quiet_NaN(),
        inf = numeric_limits<double>::infinity();
    double es[4] = {numeric_limits<double>::epsilon(), highPrecEps,
        defaultPrecEps, 0.1};
    for(int i = 0; i < 4; ++i)
    {
        double eps = es[i];
        //nan
        assert(isELess(nan, nan, nan) == false);
        assert(isELess(nan, nan, eps) == false);
        assert(isELess(nan, 1, eps) == false);
        assert(isELess(1, nan, eps) == false);
        assert(isELess(nan, inf, eps) == false);
        assert(isELess(inf, nan, eps) == false);
        //inf
        assert(isELess(inf, inf, eps) == false);
        assert(isELess(-inf, -inf, eps) == false);
        assert(isELess(-inf, inf, eps) == true);
        assert(isELess(-inf, 1, eps) == true);
        assert(isELess(-1, inf, eps) == true);
        //normal
        for(double x = numeric_limits<double>::min() * 10;
            x < numeric_limits<double>::max()/10; x *= 10)
        {
            double dx = eps * max(1.0, abs(x));
            assert(isELess(x, x + 2 * dx, eps) == true);
            assert(isELess(x, x + 0.5 * dx, eps) == false);
            assert(isELess(x - 2 * dx, x, eps) == true);
            assert(isELess(x - 0.5 * dx, x, eps) == false);
        }
    }
    DEBUG("testELessAuto passed");
}

Vector<complex<double> > slowFourierTransform(Vector<complex<double> > const&x)
{
    int n = x.getSize();
    typedef complex<double> C;
    Vector<C> result(n);
    for(int j = 0; j < n; ++j)
    {
        C c(0, 0);
        for(int i = 0; i < n; ++i)
            c += x[i] * exp((-2 * j * PI() * i/n) * C(0, 1));
        result[j] = c;
    }
    return result;
}
void FFTTestHelper(Vector<complex<double> > const& x,
    double eps = defaultPrecEps)
{//normInf works as is because complex abs function return the radius
    double normXDiff = normInf(x - IFFTGeneral(FFTGeneral(x))),
        normYDiff = normInf(slowFourierTransform(x) - FFTGeneral(x));
    if(normXDiff >= eps || normYDiff >= eps)
    {
        /*DEBUG("failed for x=");
        DEBUG(x.getSize());
        x.debug();
        DEBUG(normXDiff);
        DEBUG(normYDiff);
        DEBUG("IFFTGeneral(FFTGeneral(x))");
        IFFTGeneral(FFTGeneral(x)).debug();
        DEBUG("FFTGeneral(x)");
        FFTGeneral(x).debug();
        DEBUG("slowFourierTransform(x)");
        slowFourierTransform(x).debug();*/
        assert(false);
    }
}
void FFTTestAuto()
{
    DEBUG("FFTTestAuto");
    int nMax = 100, nn = 1000;
    for(int n = 2; n <= nMax; ++n)//fails for 1 in bits
    {
        for(int j = 0; j < nn; ++j)
        {
            Vector<complex<double> > x(n);
            for(int i = 0; i < n; ++i) x[i] = complex<double>(
                GlobalRNG().uniform(-1, 1), GlobalRNG().uniform(-1, 1));
            FFTTestHelper(x);
        }
    }
    DEBUG("FFTTestAuto passed");
}

void testQRRank1UpdateAuto()
{
    DEBUG("testQRRank1UpdateAuto");
    Matrix<double> a(3, 3);
    a(0, 0) = 1;
    a(0, 1) = 2;
    a(0, 2) = 0;
    a(1, 0) = 3;
    a(1, 1) = 4;
    a(1, 2) = 4;
    a(2, 0) = 5;
    a(2, 1) = 6;
    a(2, 2) = 3;
    Vector<double> u(3, 1), v(3, 2);
    Matrix<double> b = a - outerProduct(u, v);
    QRDecomposition qrb(b);
    qrb.rank1Update(u, v);
    //DEBUG("BQR");
    //(qrb.Q.transpose() * qrb.R).debug();
    double diff = normInf(a - (qrb.Q.transpose() * qrb.R));
    //DEBUG(diff);
    assert(diff < highPrecEps);
    DEBUG("testQRRank1UpdateAuto passed");
}

void testSVDAuto()
{
    DEBUG("testSVDAuto");
    Matrix<double> b2(3, 2);//Datta example transposed
    b2(0, 0) = 1;
    b2(0, 1) = 2;
    b2(1, 0) = 2;
    b2(1, 1) = 3;
    b2(2, 0) = 3;
    b2(2, 1) = 4;
    b2 = b2.transpose();
    //b2.debug();
    SVD s(b2);
    /*DEBUG("svd");
    s.svds.debug();
    DEBUG("U");
    s.U.debug();
    DEBUG("V");
    s.V.debug();*/

    Matrix<double> E(s.U.getRows(), s.V.getColumns());
    for(int i = 0; i < s.svds.getSize(); ++i) E(i, i) = s.svds[i];
    //DEBUG("UTEVT");
    //(s.U.transpose() * E * s.V.transpose()).debug();
    double diff = normInf(b2 - (s.U.transpose() * E * s.V.transpose()));
    //DEBUG(diff);
    assert(diff < highPrecEps);
    DEBUG("testSVDAuto passed");
}

void testSVD2Auto()
{
    DEBUG("testSVD2Auto");
    Matrix<double> b2(4, 2);//Ford example
    b2(0, 0) = 2;
    b2(0, 1) = 5;
    b2(1, 0) = 1;
    b2(1, 1) = 4;
    b2(2, 0) = -1;
    b2(2, 1) = 6;
    b2(3, 0) = 7;
    b2(3, 1) = 8;
    b2 = b2.transpose();
    //b2.debug();
    SVD s(b2);
    /*DEBUG("svd");
    s.svds.debug();
    DEBUG("U");
    s.U.debug();
    DEBUG("V");
    s.V.debug();*/

    Matrix<double> E(s.U.getRows(), s.V.getColumns());
    for(int i = 0; i < s.svds.getSize(); ++i) E(i, i) = s.svds[i];
    //DEBUG("UTEVT");
    //(s.U.transpose() * E * s.V.transpose()).debug();

    double diff = normInf(b2 - (s.U.transpose() * E * s.V.transpose()));
    //DEBUG(diff);
    assert(diff < highPrecEps);

    assert(s.rank() == 2);
    assert(int(s.norm2() * 10) == 132);
    DEBUG("testSVD2Auto passed");
}

void testTridiagonalAuto()
{
    DEBUG("testTridiagonalAuto");
    //From Burden - correct
    Matrix<double> b2(3, 3);
    b2(0, 0) = 3;
    b2(0, 1) = 1;
    b2(0, 2) = 0;
    b2(1, 0) = 1;
    b2(1, 1) = 3;
    b2(1, 2) = 1;
    b2(2, 0) = 0;
    b2(2, 1) = 1;
    b2(2, 2) = 3;
    //b2.debug();
    pair<Vector<double>, Matrix<double> > EQ = QREigenSymmetric(b2);
    Vector<double> eig = EQ.first;
    /*DEBUG("eigVa");
    eig.debug();
    DEBUG("eigVe");
    EQ.second.debug();*/
    Matrix<double> D(eig.getSize(), eig.getSize());
    for(int i = 0; i < eig.getSize(); ++i) D(i, i) = eig[i];
    //DEBUG("VTDV");
    //(EQ.second.transpose() * D * EQ.second).debug();
    double diff = normInf(b2 - (EQ.second.transpose() * D * EQ.second));
    assert(diff < highPrecEps);
    //DEBUG("VtV");
    //(EQ.second.transpose() * EQ.second).debug();
    diff = normInf(Matrix<double>::identity(3) - (EQ.second.transpose() * EQ.second));
    assert(diff < highPrecEps);
    //DEBUG("VVt");
    //(EQ.second * EQ.second.transpose()).debug();
    diff = normInf(Matrix<double>::identity(3) - (EQ.second * EQ.second.transpose()));
    assert(diff < highPrecEps);
    //test inverse
    Matrix<double> D2(eig.getSize(), eig.getSize());
    for(int i = 0; i < eig.getSize(); ++i) D2(i, i) = 1/eig[i];
    //DEBUG("VTDVVTD2V");
    //(EQ.second.transpose() * D * EQ.second * EQ.second.transpose() * D2 * EQ.second).debug();
    diff = normInf(Matrix<double>::identity(3) - (EQ.second.transpose() * D * EQ.second * EQ.second.transpose() * D2 * EQ.second));
    assert(diff < highPrecEps);

    for(int j = 0; j < eig.getSize(); ++j)
    {//test mult by matrix and eigenvalue
        Vector<double> eig1(eig.getSize());
        for(int i = 0; i < eig.getSize(); ++i) eig1[i] = EQ.second(j, i);
        //DEBUG(j);
        //eig1.debug();
        /*DEBUG("Aev");
        (b2 * eig1).debug();
        DEBUG("eev");
        (eig1 * eig[j]).debug();*/
        diff = normInf((b2 * eig1) - (eig1 * eig[j]));
        assert(diff < highPrecEps);
    }
    DEBUG("testTridiagonalAuto passed");
}

template<typename SOLVER>//LUP or QR
void testEQSolverAutoHelper(int n = 100)
{
    Matrix<double> A(n, n);
    for(int c = 0; c < n; ++c)
    {
        Vector<double> cc = GlobalRNG().randomUnitVector(n);
        for(int r = 0; r < n; ++r) A(r, c) = cc[r];
    }
    SOLVER s(A);
    Vector<double> solution = GlobalRNG().randomUnitVector(n), b = A * solution,
        x = s.solve(b), r = A * x - b;
    //DEBUG(norm(r));
    //expect backward stable always unless singular or ill-cond?
    assert(norm(r) < defaultPrecEps * max(1.0, norm(b)));
    pair<double, double> errors = estimateMatrixEquationError2Norm(A, b, x);
    /*DEBUG(errors.first);
    DEBUG(errors.second);
    DEBUG(norm(solution - x));
    DEBUG(norm(solution - x)/norm(solution));*/
    assert(norm(solution - x) <= errors.first);
    assert(norm(solution - x)/norm(solution) < errors.second);
    //inverse test
    double diff = normInf(Matrix<double>::identity(n) - A * inverse(s, n));
    assert(diff < defaultPrecEps);//high fails
}
void testEQSolverAuto()
{
    DEBUG("testEQSolverAuto");
    int n = 1;
    for(int i = 0; i < n; ++i)
    {
        testEQSolverAutoHelper<LUP<> >();
    }
    for(int i = 0; i < n; ++i)
    {
        testEQSolverAutoHelper<QRDecomposition>();
    }
    DEBUG("testEQSolverAuto passed");
}

void testTridiagSolveAuto()
{//Fausett example
    DEBUG("testTridiagSolveAuto");
    TridiagonalMatrix<double> T(4);
    T(0, 0) = 2;
    T(0, 1) = 2;
    T(1, 0) = 2;
    T(1, 1) = 4;
    T(1, 2) = 4;
    T(2, 1) = 1;
    T(2, 2) = 3;
    T(2, 3) = 3;
    T(3, 2) = 2;
    T(3, 3) = 5;
    //DEBUG("T");
    //T.debug();
    Vector<double> r(4);
    r[0] = 4;
    r[1] = 6;
    r[2] = 7;
    r[3] = 10;
    Vector<double> x = solveTridiag(T, r), answer;
    answer.append(1);
    answer.append(1);
    answer.append(0);
    answer.append(2);
    //DEBUG("x expect 1 1 0 2");
    //x.debug();
    assert(normInf(x - answer) < highPrecEps);
    DEBUG("testTridiagSolveAuto passed");
}

void testCholAuto()
{
    DEBUG("testCholAuto");
    Matrix<double> b2(3, 3);
    b2(0, 0) = 1;
    b2(0, 1) = 4;
    b2(0, 2) = 5;
    b2(1, 0) = 4;
    b2(1, 1) = 20;
    b2(1, 2) = 32;
    b2(2, 0) = 5;
    b2(2, 1) = 32;
    b2(2, 2) = 64;

    Cholesky<double> cho(b2);
    //this is LLT; note that LTL isn't the same
    double diff = normInf(b2 - cho.l * cho.l.transpose());
    assert(diff < highPrecEps);

    Vector<double> b;
    b.append(3);
    b.append(7);
    b.append(8);
    LUP<double> lup(b2);
    diff = normInf(cho.solve(b) - lup.solve(b));
    assert(diff < highPrecEps);
    DEBUG("testCholAuto passed");
}

void testDetAuto(int n = 100)
{
    DEBUG("testDetAuto");
    Matrix<double> A(n, n);
    for(int c = 0; c < n; ++c)
    {
        Vector<double> cc = GlobalRNG().randomUnitVector(n);
        for(int r = 0; r < n; ++r) A(r, c) = cc[r];
    }
    LUP<double> lup(A);
    QRDecomposition qr(A);
    double diff = abs(qr.logAbsDet() - lup.logAbsDet());
    assert(diff < defaultPrecEps);
    DEBUG("testDetAuto passed");
}

void testEigAsymetricAuto()
{
    DEBUG("testEigAsymetricAuto");
    //From Fausett
    Matrix<double> b2(4, 4);
    b2(0, 0) = 11;
    b2(0, 1) = -26;
    b2(0, 2) = 3;
    b2(0, 3) = -12;
    b2(1, 0) = 3;
    b2(1, 1) = -12;
    b2(1, 2) = 3;
    b2(1, 3) = -6;
    b2(2, 0) = 31;
    b2(2, 1) = -99;
    b2(2, 2) = 15;
    b2(2, 3) = -44;
    b2(3, 0) = 9;
    b2(3, 1) = -10;
    b2(3, 2) = -3;
    b2(3, 3) = -4;
    //b2.debug();
    pair<Vector<complex<double> >, Matrix<double> > result = QREigen(b2);
    Vector<complex<double> > eig = result.first;
    //DEBUG("eig");
    //eig.debug();
    //DEBUG("eigve");
    //result.second.debug();
    Matrix<double> D(eig.getSize(), eig.getSize());
    for(int i = 0; i < eig.getSize(); ++i) D(i, i) = eig[i].real();
    QRDecomposition qr(result.second);
    //DEBUG("V-1DV");
    //(result.second * D * inverse(qr, 4)).debug();
    double diff = normInf(b2 - (result.second * D * inverse(qr, 4)));
    assert(diff < defaultPrecEps);
    DEBUG("testEigAsymetricAuto passed");
}

void testCGAuto()
{
    DEBUG("testCGAuto");
    Matrix<double> Ad(3, 3);
    Ad(0, 0) = 5;
    Ad(0, 1) = 1;
    Ad(0, 2) = 1;
    Ad(1, 0) = 1;
    Ad(1, 1) = 4;
    Ad(1, 2) = 1;
    Ad(2, 0) = 1;
    Ad(2, 1) = 1;
    Ad(2, 2) = 6;
    //DEBUG("Ad");
    //Ad.debug();
    SparseMatrix<> A(3, 3);
    for(int r = 0; r < 3; ++r)
            for(int c = 0; c < 3; ++c) A.set(r, c, Ad(r, c));
    Vector<double> b(3);
    b[0] = 1;
    b[1] = 2;
    b[2] = 3;
    pair<Vector<double>, double> result = conjugateGradientSolve(A, b, findJacobiPreconditioner(A, false), false);
    LUP<double> lup(Ad);
    double diff = normInf(result.first - lup.solve(b));
    assert(diff < highPrecEps);
    DEBUG("testCGAuto passed");
}

void testSparseMatrixAuto()
{
    DEBUG("testSparseMatrixAuto");
    SparseMatrix<double> m(4, 4);
    Matrix<double> md(4, 4);
    //DEBUG("m");
    //m.debug();
    m.set(0, 0, 2);
    md(0, 0) = 2;
    assert(normInf(md - toDense<double>(m)) < highPrecEps);
    //DEBUG("m");
    //m.debug();
    m = SparseMatrix<double>::identity(4) * 10;
    m = m * m;
    md = Matrix<double>::identity(4) * 10;
    md = md * md;
    assert(normInf(md - toDense<double>(m)) < highPrecEps);
    for(int i = 0; i < 4; ++i) for(int j = 0; j < 4; ++j) if(i > j) m.set(i, j, 1);
    for(int i = 0; i < 4; ++i) for(int j = 0; j < 4; ++j) if(i > j) md(i, j) = 1;
    assert(normInf(md - toDense<double>(m)) < highPrecEps);
    DEBUG("testSparseMatrixAuto passed");
}

void matrixTestAllAuto()
{
    testQRRank1UpdateAuto();
    testSVDAuto();
    testSVD2Auto();
    testTridiagonalAuto();
    testEQSolverAuto();
    testTridiagSolveAuto();
    testCholAuto();
    testDetAuto();
    testEigAsymetricAuto();
    testCGAuto();
    testSparseMatrixAuto();
}

void testAllAutoNumericalMethods()
{
    DEBUG("testAllAutoNumericalMethods");
    testELessAuto();
    FFTTestAuto();
    matrixTestAllAuto();
}

}//end namespace
#endif

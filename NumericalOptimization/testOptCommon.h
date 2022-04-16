#ifndef IGMDK_TEST_OPT_COMMON_H
#define IGMDK_TEST_OPT_COMMON_H
#include "../Utils/DEBUG.h"
#include "../ExternalMemoryAlgorithms/CSV.h"
using namespace std;
namespace igmdk{

template<typename TEST_SET, typename FUNCTION> void debugResultHelperBatch(
    Vector<pair<Vector<double>, double> > const& results,
    FUNCTION const& f, Vector<Vector<string> > & matrix, int start)
{
    int k = results.getSize();
    double timediff = 1.0 * (clock() - start)/CLOCKS_PER_SEC/k;
    pair<Vector<double>, double> answer = f.getAnswer();
    double cap = -4;//"solved" if got few digits or more
    IncrementalStatistics xe, ye, gn, se;
    for(int i = 0; i < k; ++i)
    {//cap change at eps and digits at 0, digits negative for smaller is better
        pair<Vector<double>, double> const& result = results[i];
        double eps = numeric_limits<double>::epsilon(),
            relAbsXDigits = min(0.0, log10(max(eps, normInf(result.first - answer.first)/
            max(1.0, normInf(answer.first))))),
            relAbsYDigits = min(0.0, log10(max(eps, abs(result.second - answer.second)/
                max(1.0, abs(answer.second)))));
        xe.addValue(relAbsXDigits);
        ye.addValue(relAbsYDigits);
        se.addValue(relAbsYDigits <= cap);
        double gradNorm = norm(estimateGradientCD(result.first, f)),
            normalizedGradNorm = gradNorm/max(1.0, abs(result.second));
        gn.addValue(normalizedGradNorm);
    }
    double relAbsXDigits = xe.getMean(), relAbsYDigits = ye.getMean(),
        sePercentage = -se.getMean();
    DEBUG(relAbsXDigits);
    DEBUG(relAbsYDigits);
    DEBUG(sePercentage);
    DEBUG(TEST_SET::evalCount/k);//ok to round down
    matrix.lastItem().append(to_string(relAbsXDigits));
    matrix.lastItem().append(to_string(relAbsYDigits));
    matrix.lastItem().append(to_string(sePercentage));
    matrix.lastItem().append(to_string(TEST_SET::evalCount/k));
    TEST_SET::evalCount = 0;
    double normalizedGradNorm = gn.getMean();
    TEST_SET::evalCount = 0;
    DEBUG(normalizedGradNorm);
    matrix.lastItem().append(to_string(normalizedGradNorm));
    DEBUG(timediff);
    matrix.lastItem().append(to_string(timediff));
}

void createMinReport(string const& prefix,
    Vector<Vector<string> > const& matrix)
{
    int reportNumber = time(0);
    string filename = prefix + to_string(reportNumber) + ".csv";
    createCSV(matrix, filename.c_str());
    Vector<string> names;
    names.append("XError");
    names.append("YError");
    names.append("SEP");
    names.append("NEvals");
    names.append("ScaledGradNorm");
    names.append("TimeSeconds");
    createAugmentedCSVFiles(matrix, names, filename, 1);
}

}//end namespace igmdk
#endif

#ifndef IGMDK_CSV_H
#define IGMDK_CSV_H

#include "File.h"
#include "../RandomNumberGeneration/MultipleComparison.h"
using namespace std;
namespace igmdk{

void createCSV(Vector<Vector<string> > const& matrix, const char* filename)
{
    ofstream file(filename);
    assert(file);
    for(int i = 0; i < matrix.getSize(); ++i)
    {
        for(int j = 0; j < matrix[i].getSize(); ++j)
        {
            if(j > 0) file << ",";
            file << matrix[i][j];
        }
        file << endl;
    }
}

Vector<Vector<Vector<string> > > splitRegularMatrix(
    Vector<Vector<string> > const& matrix, int nMetrics)
{//must have proper number of columns in every row
    assert(nMetrics > 0);//calculate the number of row, columns,
    //and metrics
    for(int i = 0; i < matrix.getSize(); ++i)
        assert((matrix[i].getSize() - 1) % (nMetrics + 1) == 0);
    int nNewRows = matrix.getSize() + 1,
        nNewColumns = 1 + (matrix[0].getSize() - 1)/(nMetrics + 1);
    //do the splitting
    Vector<Vector<Vector<string> > > result(nMetrics,
        Vector<Vector<string> >(nNewRows, Vector<string>(nNewColumns)));
    for(int i = 0; i < nMetrics; ++i)
    {//copy over algorithms names from first row
        for(int c = 1; c < nNewColumns; ++c) result[i][0][c] =
            matrix[1][1 + (c - 1) * (nMetrics + 1)];
        for(int r = 1; r < nNewRows; ++r)
        {//copy over problem metricNames
            result[i][r][0] = matrix[r - 1][0];
            //copy over all relevant columns
            for(int c = 1; c < nNewColumns; ++c) result[i][r][c] =
                matrix[r - 1][2 + i + (c - 1) * (nMetrics + 1)];
        }
    }
    return result;
}


string cellValue(int r, int c)
{//convert cell and reference
    return "INDIRECT(ADDRESS(" + to_string(r + 1) + ";" +
        to_string(c + 1) + ";4))";
}
string fixNumber(int r, int c)
{//convert cell, reference, and convert to number if scientific notation
    return "VALUE(TRIM(" + cellValue(r, c) + "))";
}
string rankFormula(int r, int c, int c0, int cLast)
{//rank column c in range [c0, cLast] in given row
    return "RANK(" + cellValue(r, c) + ";" + cellValue(r, c0) +
        ":" + cellValue(r, cLast) + ";1)";
}
string minFormula(int r, int c0, int cLast)
    {return "MIN(" + cellValue(r, c0) + ":" + cellValue(r, cLast) + ")";}
string aveFormula(int r0, int c0, int rLast, int cLast)
{//average over a row to a column
    assert(r0 == rLast || c0 == cLast);
    return "AVERAGE(" + cellValue(r0, c0) + ":" +
        cellValue(rLast, cLast) + ")";
}

void augmentComparableMatrix(Vector<Vector<string> >& matrix, int nRepeats = 1)
{//assume first row is algorithm names + first column problem names
    int nDataRows = matrix.getSize() - 1, nColumns = matrix[0].getSize();
    //append empty row as separator
    matrix.append(Vector<string>(nColumns, ""));
    //extract numerical values in all columns
    int fixedStart = matrix.getSize();
    for(int r = 0; r < nDataRows; ++r)
    {
        Vector<string> newRow(nColumns, "");
        for(int c = 1; c < nColumns; c++)
            newRow[c] = string("=") + fixNumber(1 + r, c);
        matrix.append(newRow);
    }
    //append empty row as separator
    matrix.append(Vector<string>(nColumns, ""));
    //make rank formulas for all data points
    for(int r = 0; r < nDataRows; ++r)
    {
        Vector<string> newRow(nColumns, "");
        for(int c = 1; c < nColumns; c++) newRow[c] = string("=") +
            rankFormula(fixedStart + r, c, 1, nColumns - 1);
        matrix.append(newRow);
    }
    //average the ranks in each column
    Vector<string> newRow2(nColumns, "Ave Ranks");
    for(int c = 1; c < nColumns; c++) newRow2[c] = string("=") + aveFormula(
        matrix.getSize() - nDataRows, c, matrix.getSize() - 1, c);
    matrix.append(newRow2);
    //rank the averages
    Vector<string> newRow(nColumns, "Total Rank");
    for(int c = 1; c < nColumns; c++) newRow[c] = string("=") +
        rankFormula(matrix.getSize() - 1, c, 1, nColumns - 1);
    matrix.append(newRow);
    //find significant rankDifference
    double maxDiff = findNemenyiSignificantAveRankDiff(nColumns - 1,
        nDataRows/nRepeats);
    Vector<string> newRow3(nColumns, toStringDouble(maxDiff));
    newRow3[0] = "Significant Diff";
    matrix.append(newRow3);
    Vector<string> newRow4(nColumns, "Cutoff Rank");
    for(int c = 1; c < nColumns; c++) newRow4[c] = string("=") +
        minFormula(matrix.getSize() - 3, 1, nColumns - 1) + "+" +
        cellValue(matrix.getSize() - 1, c);
    matrix.append(newRow4);
    Vector<string> newRow5(nColumns, "Same as Best");
    for(int c = 1; c < nColumns; c++) newRow5[c] = string("=") + "IF(" +
        cellValue(matrix.getSize() - 1, c) + ">" +
        cellValue(matrix.getSize() - 4, c) + ";1;0)";
    matrix.append(newRow5);
}

void createAugmentedCSVFiles(Vector<Vector<string> > const& matrix,
    Vector<string> const& metricNames, string filename, int nRepeats = 1)
{
    Vector<Vector<Vector<string> > > pieces(
        splitRegularMatrix(matrix, metricNames.getSize()));
    for(int i = 0; i < pieces.getSize(); ++i)
    {
        augmentComparableMatrix(pieces[i], nRepeats);
        createCSV(pieces[i], (metricNames[i] + "_" + filename).c_str());
    }
}

}//end namespace
#endif

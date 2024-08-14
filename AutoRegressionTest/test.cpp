#include "../Utils/UtilsTestAuto.h"
#include "../Sorting/SortTestAuto.h"
#include "../RandomTreap/DynamicSortedSequenceTestAuto.h"
#include "../HashTable/HashTableTestAuto.h"
#include "../Heaps/HeapTestAuto.h"
#include "../Graphs/GraphsTestAuto.h"
#include "../ExternalMemoryAlgorithms/ExternalMemoryAlgorithmsTestAuto.h"
#include "../StringAlgorithms/StringAlgorithmsTestAuto.h"
#include "../Compression/CompressionTestAuto.h"
#include "../MiscAlgs/MiscAlgsTestAuto.h"
#include "../Optimization/OptTestAuto.h"
#include "../LargeNumbers/LargeNumberTestAuto.h"
#include "../ComputationalGeometry/ComputationalGeometryTestAuto.h"
#include "../ErrorCorrectingCodes/ErrorCorrectingCodesTestAuto.h"
#include "../Cryptography/CryptographyTestAuto.h"
#include "../NumericalMethods/NumericalMethodsTestAuto.h"
#include "../FinancialCalculations/FinancialCalculationsTestAuto.h"

using namespace igmdk;

int main()
{
    DEBUG("All Tests Auto");
    testAllAutoUtils();
    testAllAutoSort();
    testAllAutoDynamicSortedSequence();
    testAllAutoHashTable();
    testAllAutoHeaps();
    testAllAutoGraphs();
    testAllAutoExternalMemoryAlgorithms();
    testAllAutoStringAlgorithms();
    testAllAutoCompression();
    testAllAutoMiscAlgorithms();
    testAllAutoOpt();
    testAllAutoComputationalGeometry();
    testAllAutoErrorCorrectingCodes();
    testAllAutoCryptography();
    testAllAutoNumericalMethods();
    testAllAutoFinancialCalculations();
    DEBUG("All Tests Auto passed");

	return 0;
}

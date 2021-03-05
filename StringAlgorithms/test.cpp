#include "StringAlgorithmsTestAuto.h"

#include <string>
#include <iostream>
#include <cstdlib>
#include <ctime>
using namespace std;
using namespace igmdk;

void DDDLCS()
{
    Vector<unsigned char> x, y;
    x.append('s');
    x.append('i');
    x.append('n');
    x.append('k');

    y.append('t');
    y.append('h');
    y.append('i');
    y.append('n');
    y.append('k');
    typedef Diff<unsigned char> D;
    Vector<D::EditResult> SinkIntoThink = D::diff(x, y);

    cout << "breakpoint" << endl;

}

void DDDSuffixIndex()
{
    string test = "mississippi";
    Vector<char> temp;
    for(int i = 0; i < test.length(); ++i) temp.append(test[i]);
    SuffixIndex<char> Mississippi(temp);
    cout << "breakpoint" << endl;
}

int main()
{
	testAllAutoStringAlgorithms();

    DDDLCS();
	DDDSuffixIndex();

	return 0;
}

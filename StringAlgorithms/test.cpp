#include "StringAlgorithmsTestAuto.h"

#include <string>
#include <iostream>
#include <cstdlib>
#include <ctime>
using namespace std;
using namespace igmdk;

void testREAuto()
{
    DEBUG("testRE");
    string reS = "((A*B|AC)D)";
    DEBUG(reS);
    RegularExpressionMatcher re(reS);

    cout << "breakpoint" << endl;
    DEBUG(re.matches("ABCCBD"));
    //assert(!re.matches("ABCCBD"));
    DEBUG(re.matches("BCD"));
    assert(!re.matches("BCD"));
    DEBUG(re.matches("ABD"));
    assert(re.matches("ABD"));
    DEBUG(re.matches("ACD"));
    assert(re.matches("ACD"));
    DEBUG(re.matches("AABD"));
    assert(re.matches("AABD"));
    DEBUG("testREAuto Passed");
}



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

void testSuffixIndexAuto2()
{
    string s = "aaa";
    int n = s.length();
    Vector<char> w(n, 0);
    for(int i = 0; i < n; ++i){
        w[i] = s[i];
    }
    SuffixIndex<char> index(w);
    assert(index.sa[0] == 2);
    //assert(index.lcpa[0] == 1);
    assert(index.sa[1] == 1);
    //assert(index.lcpa[1] == 1);
    assert(index.sa[2] == 0);
    //assert(index.lcpa[2] == 2);
    for(int i = 0; i < n; ++i)
    {
        DEBUG(index.lcpa[i]);
    }
    string p = "a";
    pair<int, int> lr = index.interval((char*)p.c_str(), p.length());
    assert(lr.first == 0);
    assert(lr.second == 2);
}

int main()
{
	testSuffixIndexAuto2();
	return 0;
    testAllAutoStringAlgorithms();
    DDDLCS();
	testREAuto();
	DDDSuffixIndex();

	return 0;
}

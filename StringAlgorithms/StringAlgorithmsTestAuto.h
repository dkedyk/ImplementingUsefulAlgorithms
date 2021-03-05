#ifndef IGMDK_STRING_ALGORITHMS_TEST_AUTO_H
#define IGMDK_STRING_ALGORITHMS_TEST_AUTO_H
#include <string>
using namespace std;
#include "SubstringSearch.h"
#include "RegEx.h"
#include "Diff.h"
#include "SuffixArray.h"
#include "../Utils/Utils.h"

namespace igmdk{

void testHashQAuto()
{
    DEBUG("testHashQAuto");
	string pattern = "at";
	string text = "Four score and seven years ago our fathers brought forth"
        " on this continent a new nation ...";
	HashQ<string> matcher(pattern, pattern.length());
	for(int start = 0, i = 0; start != -1; ++i)
    {
        pair<int, int> result = matcher.findNext(text, text.length(), start);
        int match = result.first;
        start = result.second;
        if(i == 0)
        {
            assert(match == 36);
            for(int j = 0; j < pattern.length(); ++j)
                assert(text[match + j] == pattern[j]);
            assert(start = 38);
        }
        else if(i == 1)
        {
            assert(match == 82);
            for(int j = 0; j < pattern.length(); ++j)
                assert(text[match + j] == pattern[j]);
            assert(start = 84);
        }
        else if(i == 2)
        {
            assert(match == -1);
            assert(start = -1);
        }
    }
    DEBUG("testHashQAuto passed");
}

void testWuManberAuto()
{
    DEBUG("testWuManberAuto");
	string pattern = "at";
	string pattern2 = "forth";
	string text = "Four score and seven years ago our fathers brought forth"
        " on this continent a new nation ...";
	Vector<pair<string, int> > patterns;
    patterns.append(pair<string, int>(pattern, pattern.length()));
    patterns.append(pair<string, int>(pattern2, pattern2.length()));
    WuManber<string> matcher(patterns);
    for(int start = 0, j = 0; start != -1; ++j)
    {
        pair<Vector<int>, int> result = matcher.findNext(text, text.length(),
            start);
        start = result.second;
        assert(result.first.getSize() == 2);
        if(j == 0)
        {
            assert(result.first[0] == 36);
            for(int j = 0; j < pattern.length(); ++j)
                assert(text[result.first[0] + j] == pattern[j]);
            assert(result.first[1] == -1);
            assert(start = 38);
        }
        else if(j == 1)
        {
            assert(result.first[0] == -1);
            assert(result.first[1] == 51);
            for(int j = 0; j < pattern2.length(); ++j)
                assert(text[result.first[1] + j] == pattern2[j]);
            assert(start = 53);
        }
        else if(j == 2)
        {
            assert(result.first[0] == 82);
            for(int j = 0; j < pattern.length(); ++j)
                assert(text[result.first[0] + j] == pattern[j]);
            assert(result.first[1] == -1);
            assert(start = 84);
        }
        else if(j == 3)
        {
            for(int i = 0; i < patterns.getSize(); ++i)
                assert(result.first[i] == -1);
            assert(start = -1);
        }
    }
    DEBUG("testWuManberAuto passed");
}

void testREAuto()
{
    DEBUG("testRE");
    string reS = "((A*B|AC)D)";
    RegularExpressionMatcher re(reS);
    assert(!re.matches("ABCCBD"));
    assert(!re.matches("BCD"));
    assert(re.matches("ABD"));
    assert(re.matches("ACD"));
    assert(re.matches("AABD"));
    DEBUG("testREAuto Passed");
}

void testShiftAndExtendedAuto()
{
    DEBUG("testShiftAndAutoExtended");
	string pattern = "at";
	string text = "Four score and seven years ago our fathers brought forth"
        " on this continent a new nation ...";
	ShiftAndExtended matcher((unsigned char*)pattern.c_str(),
        pattern.length());
	for(int i = 0;; ++i)
    {
        int match = matcher.findNext((unsigned char*)text.c_str(),
            text.length());
        if(i == 0)
        {
            assert(match == 36);
            assert(text[match] == 'a');
            assert(text[match + 1] == 't');
        }
        else if(i == 1)
        {
            assert(match == 82);
            assert(text[match] == 'a');
            assert(text[match + 1] == 't');
        }
        else if(i == 2)
        {
            assert(match == -1);
        }
        if(match == -1) break;
    }
    DEBUG("testShiftAndExtendedAuto passed");
}

void testDiffAutoHelper(string const& a, string const& b)
{
    Vector<char> av, bv;
    for(int i = 0; i < a.length(); ++i) av.append(a[i]);
    for(int i = 0; i < b.length(); ++i) bv.append(b[i]);
    typedef Diff<char> D;
    assert(D::applyDiff(av, D::diff(av, bv)) == bv);
}
void testDiffAuto()
{
    DEBUG("testDiffAuto");
    testDiffAutoHelper("hey", "hi");
    testDiffAutoHelper("hi", "hey");
    testDiffAutoHelper("hey", "a");
    testDiffAutoHelper("a", "hey");
    testDiffAutoHelper("hi", "hi");
    DEBUG("testDiffAuto passed");
}

template<typename CHAR> struct SuffixTestComparator
{
    CHAR const*const s;
    int n;
    bool operator()(int a, int b)const
    {
        assert(a >= 0 && a < n && b >= 0 && b < n);
        Vector<CHAR> va, vb;
        for(; a < n; ++a) va.append(s[a]);
        for(; b < n; ++b) vb.append(s[b]);
        LexicographicComparator<Vector<CHAR> > c;
        return c(va, vb);
    }
};
template<typename CHAR>
int suffixLCP(int a, int b, CHAR const*const s, int n)
{
    assert(a >= 0 && a < n && b >= 0 && b < n);
    int result = 0;
    while(a < n && b < n && s[a++] == s[b++]) ++result;
    return result;
}
void testSuffixIndexAuto()
{
    DEBUG("testSuffixIndexAuto");
    int n = 10000;
    Vector<char> w(n, 0);
    for(int i = 0; i < n; ++i) w[i] = GlobalRNG().next();
    SuffixIndex<char> si(w);
    SuffixTestComparator<char> c = {w.getArray(), n};
    assert(isSorted(si.sa.getArray(), 0, n - 1, c));
    for(int i = 0; i < n; ++i)
    {
        int prev = (i == 0 ? n : i) - 1;
        assert(si.lcpa[i] ==
            suffixLCP(si.sa[prev], si.sa[i], w.getArray(), n));
    }
    DEBUG("testSuffixIndexAuto passed");
}

void testSuffixIndexAuto2()
{
    DEBUG("testSuffixIndexAuto2");
    string s = "aaa";
    int n = s.length();
    Vector<char> w(n, 0);
    for(int i = 0; i < n; ++i) w[i] = s[i];
    SuffixIndex<char> index(w);
    assert(index.sa[0] == 2);//"a"
    assert(index.lcpa[0] == 1);//"aaa" vs "a"
    assert(index.sa[1] == 1);//"aa"
    assert(index.lcpa[1] == 1);//"a" vs "aa"
    assert(index.sa[2] == 0);//"aaa"
    assert(index.lcpa[2] == 2);//"aa" vs "aaa"
    string p = "a";
    pair<int, int> lr = index.interval((char*)p.c_str(), p.length());
    assert(lr.first == 0);
    assert(lr.second == 2);
    DEBUG("testSuffixIndexAuto2 passed");
}

template<typename CHAR> struct BWTTestComparator
{
    CHAR const*const s;
    int n;
    bool operator()(int a, int b)const
    {
        assert(a >= 0 && a < n && b >= 0 && b < n);
        Vector<CHAR> va, vb;
        for(int i = 0; i < n; ++i) va.append(s[(a + i) % n]);
        for(int i = 0; i < n; ++i) vb.append(s[(b + i) % n]);
        LexicographicComparator<Vector<CHAR> > c;
        return c(va, vb);
    }
};
void testBWTAuto()
{
    DEBUG("timeBWTAuto");
    int n = 10000;
    Vector<char> w(n, 0);
    for(int i = 0; i < n; ++i) w[i] = GlobalRNG().next();
    Vector<int> sa = suffixArray<BWTRank>(w.getArray(), n);
    BWTTestComparator<char> c = {w.getArray(), n};
    assert(isSorted(sa.getArray(), 0, sa.getSize() - 1, c));

    DEBUG("testBWTAuto passed");
}

void testAllAutoStringAlgorithms()
{
    DEBUG("testAllAutoStringAlgorithms");
    testHashQAuto();
    testWuManberAuto();
    testREAuto();
    testShiftAndExtendedAuto();
    testDiffAuto();
    testSuffixIndexAuto();
    testSuffixIndexAuto2();
    testBWTAuto();
}

}//end namespace
#endif

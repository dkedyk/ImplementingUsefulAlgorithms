#include <iostream>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include "ScrapTrie.h"
#include "LeftLeaningRedBlackTree.h"
#include "../RandomNumberGeneration/Random.h"
using namespace std;
using namespace igmdk;

void timeRT()
{
	LeftLeaningRedBlackTree<int, int> trie;
	//PatriciaTrie<int, int> trie;
	int N = 1500000;
	for(int i = 0; i < N; ++i)
	{
	    //DEBUG("BEGIN");
	    //DEBUG(i);
		trie.insert(i, i);
		//DEBUG("END");
	}
	//system("PAUSE");
	for(int j = 0; j < 1; ++j)
	{
		for(int i = 0; i < N; ++i)
		{
			assert(trie.find(i));
			//DEBUG(*trie.find(i));
			//cout << *trie.find(i) << endl;
			assert(*trie.find(i) == i);
			trie.remove(i);
			/*for(int k = 0; k <=i; ++k)
			{
				assert(!trie.find(k));
			}
			for(int k = i+1; k < 1500; ++k)
			{
				assert(trie.find(k));
			}*/
		}
	}
}

void timeRT2()
{
    enum{N = 10000};
    int data[N];
    for(int i = 0; i < N; ++i) data[i] = GlobalRNG().next();
	LeftLeaningRedBlackTree<int, int> trie;
	//PatriciaTrie<int, int> trie;
	for(int i = 0; i < N; ++i)
	{
		trie.insert(data[i], i);

	}
	for(int j = 0; j < 10000; ++j)
	{
		for(int i = 0; i < N; ++i)
		{
			assert(trie.find(data[i]));
			//cout << *trie.find(i) << endl;
			assert(*trie.find(data[i]) == i);
			//trie.remove(i);
			/*for(int k = 0; k <=i; ++k)
			{
				assert(!trie.find(k));
			}
			for(int k = i+1; k < 1500; ++k)
			{
				assert(trie.find(k));
			}*/
		}
	}
}

int main()
{
	clock_t start = clock();

	timeRT();
	int tFL = (clock() - start);
    cout << "FL: "<<tFL << endl;
	//system("PAUSE");
	return 0;
}

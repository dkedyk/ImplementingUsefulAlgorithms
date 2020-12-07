#include <cassert>
#include <cstdlib>
#include <new>
#include <iostream>
#include <cmath>
#include <ctime>
#include "SortedLinkedList.h"
using namespace std;
using namespace igmdk;

void timeSRT()
{
	typedef pair<int, int> Item;
	LinkedList<Item> ll;
	for(int i = 100; i>=0; --i)
	{
		ll.prepend(Item(i, i));
	}
	ll.sort();
	for(LinkedList<Item>::Node* iter = ll.root;iter; iter = iter->next)
        cout << "key" <<iter->item.first<< endl;
}

int main()
{
	clock_t start = clock();


	timeSRT();
	int tFL = (clock() - start);
    cout << "FL: "<<tFL << endl;
	system("PAUSE");
	return 0;
}

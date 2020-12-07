#include <cassert>
#include <cstdlib>
#include <new>
#include <iostream>
#include <cmath>
#include <ctime>
#include "PairingHeap.h"
#include "ScrapHeap.h"
#include "SymmetricMinMaxHeap.h"
using namespace igmdk;

void timeSRT()
{
	IndexedPaHeap<int> heap;
	//IndexedPointerHeap<int> heap;
	//BucketQueue<int> heap(9);
	//SymmetricMinMaxHeap<int> heap;
	//PairingHeap<int> heap;
	int N = 15000000;
	//PairingHeap<int>::PairingHandle a[500000];
	for(int i = 0; i < N; ++i)
	{
		heap.insert(rand()%10, i);
		//heap.insert(i);
		//a[i] = rand();
		//heap.insert(a[i]);
	}
	//system("PAUSE");
	for(int i = 0; i < N; ++i)
	{
		//heap.decreaseKey(p[i], a[i]->element - i);
		//heap.changeKey(i, a[i] - i);
	}
	for(int i = 0; i < N; ++i)
	{
		//cout <<"&&&"<< heap.getMax() << endl;
		//heap.deleteMax();
		//if(heap.isEmpty())break;
		//cout <<"&&&"<< heap.getMin() << endl;
		heap.deleteMin();
		//heap.map.remove(i);
	}
	//cout << "next" << endl;
	/*for(int i = 0; i < 10; ++i)
	{
		heap.insert(rand());
	}
	for(int i = 0; i < 10; ++i)
	{
		cout << heap.getMin() << endl;
		heap.deleteMin();
	}*/
	/*int a[10];
	for(int i = 0; i < 10; ++i)
	{
		a[i] = rand();
	}
	Heap<int> heap(a,10,10);
	for(int i = 0; i < 10; ++i)
	{
		cout << a[i] << endl;
	}
	for(int i = 0; i < 10; ++i)
	{
		heap.changeKey(i, rand());
	}
	//heap.changeKey(5, -10);
	for(int i = 0; i < 10; ++i)
	{
		//cout << a[i] << endl;
	}
	//heap.changeKey(5, 100000000);
	for(int i = 0; i < 10; ++i)
	{
		cout << a[i] << endl;
	}
	for(int i = 0; i < 10; ++i)
	{
		cout << heap.getMin() << endl;
		heap.deleteMin();
	}*/
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

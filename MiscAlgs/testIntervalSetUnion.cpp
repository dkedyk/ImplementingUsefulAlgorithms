#include "IntervalSetUnion.h"
#include <iostream>
using namespace igmdk;

int main()
{
	IntervalSetUnion iu;

	iu.split(5);
	DEBUG(iu.find(2));
	iu.split(3);
	DEBUG(iu.find(2));
	iu.merge(3);
	DEBUG(iu.find(2));
	return 0;
}

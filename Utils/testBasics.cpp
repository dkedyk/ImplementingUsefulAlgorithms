#include "Stack.h"
#include "Queue.h"
#include "GCFreeList.h"
#include "UtilsTestAuto.h"
using namespace igmdk;

void DDDVector()
{
    Vector<int> Vector0to4;
    for(int i = 0; i < 5; ++i ) Vector0to4.append(i);

    cout << "breakpoint" << endl;
}

void DDDStack()
{
    Stack<int> Stack0to4;
    for(int i = 0; i < 5; ++i ) Stack0to4.push(i);

    cout << "breakpoint" << endl;
}

void DDDQueue()
{
    Queue<int> Queue1to4;
    for(int i = 0; i < 5; ++i ) Queue1to4.push(i);
    Queue1to4.pop();

    cout << "breakpoint" << endl;
}

void DDDList()
{
    SimpleDoublyLinkedList<int> List0to2;
    for(int i = 2; i >= 0; --i ) List0to2.prepend(i);

    cout << "breakpoint" << endl;
}

void DDDFreelist()
{
    Freelist<int> Freelist0to14R5(8);
    int* items[15];
    DEBUG("a");
    for(int i = 0; i < 15; ++i ) items[i] = new(Freelist0to14R5.allocate())int(i);
    DEBUG("r");
    for(int i = 0; i < 5; ++i ) Freelist0to14R5.remove(items[i]);

    cout << "breakpoint" << endl;
}



int main()
{
    testAllAutoUtils();
    //return 0;
    DDDVector();
    //DDDDeque();
    DDDStack();
    DDDQueue();
    DDDList();
    DDDFreelist();
    return 0;
}

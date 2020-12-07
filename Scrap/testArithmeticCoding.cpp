#include "ArithmeticCoding.h"
#include <iostream>
#include <ctime>
using namespace igmdk;

struct ugly
{
    unsigned char a:4;
    unsigned char b:4;
    unsigned char c;
};

Vector<unsigned char> compress(Vector<unsigned char> const& byteArray)
{
    return Arithmetic::AdaptiveCompress(byteArray);
}

Vector<unsigned char> uncompress(Vector<unsigned char> const& byteArray)
{
    return Arithmetic::AdaptiveUncompress(byteArray);
}

void timeRT()
{
    char text[] = "hello arithmetic coding";
    Vector<unsigned char> uncompressed0;
    for(int i = 0; i < sizeof(text)-1; ++i) uncompressed0.append(text[i]);
    Vector<unsigned char> uncompressed(uncompressed0);

    cout << "uncompressed.size" << uncompressed.getSize() << endl;
    Vector<unsigned char> compressed(compress(uncompressed));
    cout << "compressed.size" << compressed.getSize() << endl;
    Vector<unsigned char> uncompressed2(uncompress(compressed));
    cout << "uncompressed2.size" << uncompressed2.getSize() << endl;
    assert(uncompressed.getSize() ==  uncompressed2.getSize());
    for(int i = 0; i < uncompressed2.getSize(); ++i) assert(uncompressed2[i] == uncompressed[i]);
}

int main()
{
	clock_t start = clock();
	//for(;;)
	for(int i =0; i < 1; ++i)
	{
	    timeRT();
	}

	int tFL = (clock() - start);
    cout << "FL: "<<tFL << endl;
	return 0;
}

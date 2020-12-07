#ifndef IGMDK_REED_SOLOMON_H
#define IGMDK_REED_SOLOMON_H
#include "../Utils/Vector.h"
#include "../Utils/Bits.h"
namespace igmdk{

class GF2mArithmetic
{
    int n;
    Vector<int> el2p, p2el;
public:
    int one()const{return 1;}
    int alpha()const{return 2;}
    GF2mArithmetic(int primPoly)
    {//m is the highest set bit of polynomial
        int m = lgFloor(primPoly);
        assert(m <= 16);//avoid using too much memory
        n = twoPower(m);
        el2p = p2el = Vector<int>(n - 1);//0 has no corresponding power
        p2el[0] = 1;//a^0 = 1, a^(n - 1) also 1 so don't store it, and
        //implicitly el2p[1] = 0
        for(int p = 1; p < n - 1; ++p)
        {//calculate a^p from a^(p - 1)
            int e = p2el[p - 1] << 1;//multiply by x
            if(e >= n) e = sub(e, primPoly);//reduce if needed
            el2p[e - 1] = p;
            p2el[p] = e;
        }
    }
    int elementToPower(int x)const
    {
        assert(x > 0 && x < n);
        return el2p[x - 1];
    }
    int powerToElement(int x)const
    {
        assert(x >= 0 && x < n - 1);
        return p2el[x];
    }//both + and - just xor
    int add(int a, int b)const{return a ^ b;}
    int sub(int a, int b)const{return add(a, b);}
    int mult(int a, int b)const
    {//add in power basis and convert back
        return a == 0 || b == 0 ? 0 :
            powerToElement((elementToPower(a) + elementToPower(b)) % (n - 1));
    }
    int div(int a, int b)const
    {//subtract in power basis and convert back
        assert(b != 0);
        return a == 0 ? 0 : powerToElement((elementToPower(a) + (n - 1) -
            elementToPower(b)) % (n - 1));
    }
};

template<typename ITEM, typename ARITHMETIC>
struct Poly: public ArithmeticType<Poly<ITEM, ARITHMETIC> >
{
    Vector<ITEM> storage;
    ARITHMETIC ari;
public:
    int getSize()const{return storage.getSize();}
    int degree()const{return getSize() - 1;}
    Poly(ARITHMETIC const& theAri, Vector<ITEM> const& coefs =//default is 0
         Vector<ITEM>(1, 0)): ari(theAri), storage(coefs)
    {
        assert(getSize() > 0);
        trim();
    }
    static Poly zero(ARITHMETIC const& theAri){return Poly(theAri);}
    ITEM const& operator[](int i)const{return storage[i];}
    void trim()
        {while(getSize() > 1 && storage.lastItem() == 0) storage.removeLast();}
    Poly& operator+=(Poly const& rhs)
    {//add term-by-term, no carry
        while(getSize() < rhs.getSize()) storage.append(0);
        for(int i = 0; i < min(getSize(), rhs.getSize()); ++i)
            storage[i] = ari.add(storage[i], rhs[i]);
        trim();
        return *this;
    }
    Poly& operator-=(Poly const& rhs)
    {//subtract term-by-term, no carry
        while(getSize() < rhs.getSize()) storage.append(0);
        for(int i = 0; i < min(getSize(), rhs.getSize()); ++i)
            storage[i] = ari.sub(storage[i], rhs[i]);
        trim();
        return *this;
    }
    Poly& operator*=(ITEM const& a)
    {
        for(int i = 0; i < getSize(); ++i)
            storage[i] = ari.mult(storage[i], a);
        trim();
        return *this;
    }
    Poly operator*(ITEM const& a)const
    {
        Poly temp(*this);
        temp *= a;
        return temp;
    }
    Poly& operator<<=(int p)
    {
        assert(p >= 0);
        if(p > 0)
        {
            for(int i = 0; i < p; ++i) storage.append(0);
            for(int i = getSize() - 1; i >= p; --i)
            {
                storage[i] = storage[i - p];
                storage[i - p] = 0;
            }
        }
        return *this;
    }
    Poly& operator>>=(int p)
    {
        assert(p >= 0);
        if(p >= getSize()) storage = Vector<ITEM>(1);
        if(p > 0)
        {
            for(int i = 0; i < getSize() - p; ++i) storage[i] = storage[i + p];
            for(int i = 0; i < p; ++i) storage.removeLast();
        }
        return *this;
    }
    Poly& operator*=(Poly const& rhs)
    {//multiply each term of rhs and sum up
        Poly temp(*this);
        *this *= rhs[0];
        for(int i = 1; i < rhs.getSize(); ++i)
        {
            temp <<= 1;
            *this += temp * rhs[i];
        }
        return *this;
    }
    static Poly makeX(ARITHMETIC const& ari)
    {//x = 1 * x + 0 * 1
        Vector<ITEM> coefs(2);
        coefs[1] = ari.one();
        return Poly(ari, coefs);;
    }
    Poly& reduce(Poly const& rhs, Poly& q)
    {//quotient-remainder division, similar to numbers
        assert(rhs.storage.lastItem() != 0 && q == zero(ari));
        Poly one(ari, Vector<ITEM>(1, ari.one()));
        while(getSize() >= rhs.getSize())
        {//field guarantees exact division
            int diff = getSize() - rhs.getSize();
            ITEM temp2 = ari.div(storage.lastItem(), rhs.storage.lastItem());
            assert(storage.lastItem() ==
                ari.mult(temp2, rhs.storage.lastItem()));
            *this -= (rhs << diff) * temp2;
            q += (one << diff) * temp2;
        }
        return *this;
    }
    Poly& operator%=(Poly const& rhs)
    {
        Poly dummyQ(ari);
        return reduce(rhs, dummyQ);
    }
    bool operator==(Poly const& rhs)const
    {
        if(getSize() != rhs.getSize()) return false;
        for(int i = 0; i < getSize(); ++i)if(storage[i] != rhs[i])return false;
        return true;
    }
    ITEM eval(ITEM const& x)const
    {//Horner's algorithm
        ITEM result = storage[0], xpower = x;
        for(int i = 1; i < getSize(); ++i)
        {
            result = ari.add(result, ari.mult(xpower, storage[i]));
            xpower = ari.mult(xpower, x);
        }
        return result;
    }
    Poly formalDeriv()const
    {
        Vector<ITEM> coefs(getSize() - 1);
        for(int i = 0; i < coefs.getSize(); ++i)
            for(int j = 0; j < i + 1; ++j)
                coefs[i] = ari.add(coefs[i], storage[i + 1]);
        return Poly(ari, coefs);
    }

    void debug()
    {
        DEBUG(getSize());
        for(int i = 0; i < getSize(); ++i)
        {
            DEBUG(int(storage[i]));
        }
    }
};

class ReedSolomon
{
    int n, k;
    GF2mArithmetic ari;
    typedef Poly<unsigned char, GF2mArithmetic> P;
    typedef Vector<unsigned char> V;
    P generator;
    pair<P, P> findLocatorAndEvaluator(P const& syndromePoly, int t)const
    {
        P evPrev(ari, V(1, ari.one())), ev = syndromePoly,
            locPrev = P::zero(ari), loc = evPrev;
        evPrev <<= t;
        while(ev.degree() >= t/2)
        {
            P q(ari);
            evPrev.reduce(ev, q);
            swap(ev, evPrev);
            locPrev -= q * loc;
            swap(loc, locPrev);
        }//normalize them
        if(loc != P::zero(ari))
        {
            int normalizer = ari.div(ari.one(), loc[0]);
            loc *= normalizer;
            ev *= normalizer;
        }
        return make_pair(loc, ev);
    }
public:
    ReedSolomon(int theK = 223, int primPoly = 301): ari(primPoly),
        generator(ari, V(1, 1)), k(theK),
        n(twoPower(lgFloor(primPoly)) - 1)
    {
        assert(k < n && numeric_limits<unsigned char>::digits == 8);
        P x = P::makeX(ari);
        for(int i = 0, aPower = ari.alpha(); i < n - k; ++i)
        {
            generator *= (x - P(ari, V(1, aPower)));
            aPower = ari.mult(aPower, ari.alpha());
        }
        assert(generator.getSize() == n - k + 1);
    }
    V lengthPadBlock(V block)
    {
        assert(block.getSize() < k);
        block.append(block.getSize());
        while(block.getSize() < k) block.append(0);
        return block;
    }
    pair<V, bool> lengthUnpadBlock(V block)
    {
        assert(block.getSize() == k);
        while(block.getSize() >= 0 && block.lastItem() == 0)block.removeLast();
        bool correct = block.getSize() >= 0 &&
            block.lastItem() == block.getSize() - 1;
        assert(correct);
        if(correct) block.removeLast();
        return make_pair(block, correct);
    }
    V encodeBlock(V const& block)const
    {
        assert(block.getSize() == k);
        P c(ari, block);//init c
        c <<= (n - k);//make space for code
        c += c % generator;//add code
        //beware of poly trim if block is 0
        while(c.storage.getSize() < n) c.storage.append(0);
        return c.storage;
    }
    pair<V, bool> decodeBlock(V const& code)const
    {//calculate syndrome polynomial
        assert(code.getSize() == n);
        P c(ari, code);
        int t = n - k, aPower = ari.alpha();
        V syndromes(t);
        for(int i = 0; i < t; ++i)
        {
            syndromes[i] = c.eval(aPower);
            aPower = ari.mult(aPower, ari.alpha());
        }
        P s(ari, syndromes);
        if(s == P::zero(ari))//no error if yes
        {//take out check data and restore trimmed 0's
            c >>= t;
            while(c.storage.getSize() < k) c.storage.append(0);
            return make_pair(c.storage, true);
        }//find locator and evaluator polys
        pair<P, P> locEv = findLocatorAndEvaluator(s, t);
        if(locEv.first == P::zero(ari)) return make_pair(code, false);
        //find locator roots
        V roots;
        for(int i = 1; i < n + 1; ++i)
            if(locEv.first.eval(i) == 0) roots.append(i);
        if(roots.getSize() == 0) return make_pair(code, false);
        //find error values
        P fd = locEv.first.formalDeriv();
        V errors;
        for(int i = 0; i < roots.getSize(); ++i) errors.append(ari.sub(0,
            ari.div(locEv.second.eval(roots[i]), fd.eval(roots[i]))));
        //correct errors
        while(c.storage.getSize() < n) c.storage.append(0);
        for(int i = 0; i < roots.getSize(); ++i)
        {
            int location = ari.elementToPower(ari.div(ari.one(), roots[i]));
            assert(location < c.getSize());
            c.storage[location] = ari.add(c.storage[location], errors[i]);
        }
        if(c % generator != P::zero(ari)) return make_pair(code, false);
        c >>= t;
        return make_pair(c.storage, true);
    }
};

}//end namespace
#endif

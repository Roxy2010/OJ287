#ifndef _tvalarray_h
#define _tvalarray_h

#if defined(DEBUG) || defined(ARRAY_BOUNDS_CHECK)
#include <iostream>
#include <cstdlib> /* abort() */
#ifndef _COREDEBUG_STR
#define _COREDEBUG_STR "To debug the core, use this command:\n"\
"    gdb -c core /path/to/executable             (gdb) bt"
#endif
#endif


#if defined(DEBUG) || defined(ARRAY_BOUNDS_CHECK)
inline bool resize_warning()
{
    std::cerr<<"General information: do not use tvalarray< tvalarray<T> >::resize(unsigned, tvalarray<T>),\n"
             <<"because does not work!\n"
             <<"In case of hyper-arrays use resize(tvalarray< tvalarray<T> >&, unsigned, tvalarray<T>) instead!\n";
    return true;
}

inline bool display_resize_warning()
{
    static const bool resize_warning_buffer=resize_warning();
    return resize_warning_buffer;
}
#endif


/**
 * Array class with debugging support and a minimal set of
 * std::valarray-compatible methods.
 *
 * @version 0.5, 12/22/2008
 * @author Peter Csizmadia (2000-2003,2007)
 */
template <class T>
class tvalarray
{
private:
    T* x;
    unsigned n;
    const bool ownsData;

protected:
    /** Create an array that uses the specified low level array. */
    tvalarray(T* array, unsigned _n, bool realloc);

public:
    /** Create an empty array. */
    tvalarray() : ownsData(true) {
	x = 0;
       	n = 0;
    }

    /** Create an array with n elements. */
    tvalarray(unsigned _n) : ownsData(true) {
	if(_n > 0) {
	    n=_n;
	    x = new T[n];
	    *this = T(); // initialize elements
	} else {
	    n = 0;
	    x = 0;
	}
    }

    /** Create an array with n initialized elements. */
    tvalarray(const T& value, unsigned _n);

    /** Create an array with n initialized elements. */
    tvalarray(const T* values, unsigned _n);

    /** Copy constructor. */
    tvalarray(const tvalarray& a) : ownsData(true) {
	n = a.n;
	if(n > 0) {
	    x = new T[n];
	    const T* y = a.x;
	    for(unsigned i = 0; i < n; ++i) {
		x[i] = y[i];
	    }
	} else {
	    n = 0;
	    x = 0;
	}
    }

    ~tvalarray() {
	if(x != 0 && ownsData) {
	    delete[] x;
	}
    }

    /** Gets the number of elements. */
    unsigned size() const {
	return n;
    }

    /** Gets the underlying array. */
    operator T*() {
	return x;
    }

    /** Gets the underlying array. */
    operator const T*() const {
	return x;
    }

    /** Copy (but do not resize!!!). */
    tvalarray& operator =(const tvalarray& a) {
#if defined(DEBUG) || defined(ARRAY_BOUNDS_CHECK)
	if(n != a.n) {
	    std::cerr << "Fatal error: tvalarray a = b with different sizes ("
		 << n << " != " << a.n << ")\n" _COREDEBUG_STR << std::endl;
	    ::abort();
	}
#endif
	const T* y = a.x;
	for(unsigned i = 0; i < n; ++i) {
	    x[i] = y[i];
	}
       	return *this;
    }

    /**
     * Copy (also realloc if necessary). The instance for Tp!=tvalarray<>.
     */
    template <class Tp>
    friend void copy(const tvalarray<Tp>& source, tvalarray<Tp>& dest);

    /**
     * Copy (also realloc if necessary). The instance for Tp==tvalarray<>.
     */
    template <class Tp>
    friend void copy(const tvalarray< tvalarray<Tp> >& source, tvalarray< tvalarray<Tp> >& dest);

    /** Addition assignment. */
    tvalarray& operator +=(const tvalarray& a) {
#if defined(DEBUG) || defined(ARRAY_BOUNDS_CHECK)
	if(n != a.n) {
	    std::cerr << "Fatal error: tvalarray a += b with different sizes ("
		 << n << " != " << a.n << ")\n" _COREDEBUG_STR << std::endl;
	    ::abort();
	}
#endif
	const T* y = a.x;
	for(unsigned i = 0; i < n; ++i) {
	    x[i] += y[i];
	}
       	return *this;
    }

    /** Addition assignment. */
    tvalarray& operator +=(const T& c) {
	for(unsigned i = 0; i < n; ++i) {
	    x[i] += c;
	}
	return *this;
    }

    /** Bitwise exclusive OR assignment. */
    tvalarray& operator ^=(const tvalarray& a) {
#if defined(DEBUG) || defined(ARRAY_BOUNDS_CHECK)
	if(n != a.n) {
	    std::cerr << "Fatal error: tvalarray a ^= b with different sizes ("
		 << n << " != " << a.n << ")\n" _COREDEBUG_STR << std::endl;
	    ::abort();
	}
#endif
	const T* y = a.x;
	for(unsigned i = 0; i < n; ++i) {
	    x[i] ^= y[i];
	}
       	return *this;
    }

    /** Bitwise exclusive OR assignment. */
    tvalarray& operator ^=(const T& c) {
	for(unsigned i = 0; i < n; ++i) {
	    x[i] ^= c;
	}
	return *this;
    }

    /** Bitwise AND assignment. */
    tvalarray& operator &=(const tvalarray& a) {
#if defined(DEBUG) || defined(ARRAY_BOUNDS_CHECK)
	if(n != a.n) {
	    std::cerr << "Fatal error: tvalarray a &= b with different sizes ("
		 << n << " != " << a.n << ")\n" _COREDEBUG_STR << std::endl;
	    ::abort();
	}
#endif
	const T* y = a.x;
	for(unsigned i = 0; i < n; ++i) {
	    x[i] &= y[i];
	}
       	return *this;
    }

    /** Bitwise AND assignment. */
    tvalarray& operator &=(const T& c) {
	for(unsigned i = 0; i < n; ++i) {
	    x[i] &= c;
	}
	return *this;
    }

    /** Bitwise inclusive OR assignment. */
    tvalarray& operator |=(const tvalarray& a) {
#if defined(DEBUG) || defined(ARRAY_BOUNDS_CHECK)
	if(n != a.n) {
	    std::cerr << "Fatal error: tvalarray a |= b with different sizes ("
		 << n << " != " << a.n << ")\n" _COREDEBUG_STR << std::endl;
	    ::abort();
	}
#endif
	const T* y = a.x;
	for(unsigned i = 0; i < n; ++i) {
	    x[i] |= y[i];
	}
       	return *this;
    }

    /** Bitwise inclusive OR assignment. */
    tvalarray& operator |=(const T& c) {
	for(unsigned i = 0; i < n; ++i) {
	    x[i] |= c;
	}
	return *this;
    }

    /** Bitwise left shift assignment. */
    tvalarray& operator <<=(const tvalarray& a) {
#if defined(DEBUG) || defined(ARRAY_BOUNDS_CHECK)
	if(n != a.n) {
	    std::cerr << "Fatal error: tvalarray a <<= b with different sizes ("
		 << n << " != " << a.n << ")\n" _COREDEBUG_STR << std::endl;
	    ::abort();
	}
#endif
	const T* y = a.x;
	for(unsigned i = 0; i < n; ++i) {
	    x[i] <<= y[i];
	}
       	return *this;
    }

    /** Bitwise left shift assignment. */
    tvalarray& operator <<=(const T& c) {
	for(unsigned i = 0; i < n; ++i) {
	    x[i] <<= c;
	}
	return *this;
    }

    /** Bitwise right shift assignment. */
    tvalarray& operator >>=(const tvalarray& a) {
#if defined(DEBUG) || defined(ARRAY_BOUNDS_CHECK)
	if(n != a.n) {
	    std::cerr << "Fatal error: tvalarray a >>= b with different sizes ("
		 << n << " != " << a.n << ")\n" _COREDEBUG_STR << std::endl;
	    ::abort();
	}
#endif
	const T* y = a.x;
	for(unsigned i = 0; i < n; ++i) {
	    x[i] >>= y[i];
	}
       	return *this;
    }

    /** Bitwise right shift assignment. */
    tvalarray& operator >>=(const T& c) {
	for(unsigned i = 0; i < n; ++i) {
	    x[i] >>= c;
	}
	return *this;
    }

    /** Multiplication assignment. */
    tvalarray& operator *=(const tvalarray& a) {
#if defined(DEBUG) || defined(ARRAY_BOUNDS_CHECK)
	if(n != a.n) {
	    std::cerr << "Fatal error: tvalarray a *= b with different sizes ("
		 << n << " != " << a.n << ")\n" _COREDEBUG_STR << std::endl;
	    ::abort();
	}
#endif
	const T* y = a.x;
	for(unsigned i = 0; i < n; ++i) {
	    x[i] *= y[i];
	}
       	return *this;
    }

    /** Multiplication assignment. */
    tvalarray& operator *=(const T& c) {
	for(unsigned i = 0; i < n; ++i) {
	    x[i] *= c;
	}
	return *this;
    }

    /** Division assignment. */
    tvalarray& operator /=(const tvalarray& a) {
#if defined(DEBUG) || defined(ARRAY_BOUNDS_CHECK)
	if(n != a.n) {
	    std::cerr << "Fatal error: tvalarray a /= b with different sizes ("
		 << n << " != " << a.n << ")\n" _COREDEBUG_STR << std::endl;
	    ::abort();
	}
#endif
	const T* y = a.x;
	for(unsigned i = 0; i < n; ++i) {
	    x[i] /= y[i];
	}
       	return *this;
    }

    /** Division assignment. */
    tvalarray& operator /=(const T& c) {
	for(unsigned i = 0; i < n; ++i) {
	    x[i] /= c;
	}
	return *this;
    }

    /** Remainder assignment. */
    tvalarray& operator %=(const tvalarray& a) {
#if defined(DEBUG) || defined(ARRAY_BOUNDS_CHECK)
	if(n != a.n) {
	    std::cerr << "Fatal error: tvalarray a %= b with different sizes ("
		 << n << " != " << a.n << ")\n" _COREDEBUG_STR << std::endl;
	    ::abort();
	}
#endif
	const T* y = a.x;
	for(unsigned i = 0; i < n; ++i) {
	    x[i] %= y[i];
	}
       	return *this;
    }

    /** Remainder assignment. */
    tvalarray& operator %=(const T& c) {
	for(unsigned i = 0; i < n; ++i) {
	    x[i] %= c;
	}
	return *this;
    }

    /** Assignment operation. */
    tvalarray& operator =(const T& c) {
	for(unsigned i = 0; i < n; ++i) {
	    x[i] = c;
	}
	return *this;
    }

    /** Subtraction assignment. */
    tvalarray& operator -=(const T& c) {
	for(unsigned i = 0; i < n; ++i) {
	    x[i] -= c;
	}
	return *this;
    }

    /** Subtraction assignment. */
    tvalarray& operator -=(const tvalarray& a) {
#if defined(DEBUG) || defined(ARRAY_BOUNDS_CHECK)
	if(n != a.n) {
	    std::cerr << "Fatal error: tvalarray a -= b with different sizes ("
		 << n << " != " << a.n << ")\n" _COREDEBUG_STR << std::endl;
	    ::abort();
	}
#endif
	const T* y = a.x;
	for(unsigned i = 0; i < n; ++i) {
	    x[i] -= y[i];
	}
       	return *this;
    }

    /** Resize array and initialize all elements. */
    void resize(unsigned _n, T value = T()) {
#if defined(DEBUG) || defined(ARRAY_BOUNDS_CHECK)
	display_resize_warning();
#endif
	if(_n > n) {
	    realloc(_n);
	}
	if((n = _n) > 0) {
	    *this = value; // initialize elements (does not work for T==tvalarray<>)
	}
    }

    /** Resize array and initialize all elements. */
    template <class Tp>
    friend void resize(tvalarray<Tp>& arr, unsigned _n, Tp value /*= Tp()*/);

    /** Resize array and initialize all elements in hyper-arrays. */
    // template <class Tp>
    // friend void resize(tvalarray<tvalarray<Tp>>& arr, unsigned _n, tvalarray<Tp> value = tvalarray<Tp>{});

    /** Gets the kth element. */
    T operator [](int k) const {
#if defined(DEBUG) || defined(ARRAY_BOUNDS_CHECK)
	if(k < 0 || k >= (int)n) {
	    std::cerr << "Fatal error: tvalarray index " << k
		 << " out of range [0, " << n << "[\n" _COREDEBUG_STR << std::endl;
	    ::abort();
	}
#endif
	return x[k];
    }

    /** Gets the kth element. */
    T& operator [](int k) {
#if defined(DEBUG) || defined(ARRAY_BOUNDS_CHECK)
	if(k < 0 || k >= (int)n) {
	    std::cerr << "Fatal error: tvalarray index " << k
		 << " out of range [0, " << n << "[\n" _COREDEBUG_STR << std::endl;
	    ::abort();
	}
#endif
	return x[k];
    }

    /** Returns a shifted array. */
    tvalarray shift(int k) const;

    /** Returns a rotated array. */
    tvalarray cshift(int k) const;

public:
    /** (Re)allocate array. */
    void realloc(int _n) {
	if(ownsData) {
	    if ( x!=0 ) delete [] x;
	    if(_n > 0) { n = _n; x = new T[n]; }
	    else { n=0; x=0; }
	} else {
#if defined(DEBUG) || defined(ARRAY_BOUNDS_CHECK)
	    std::cerr << "Fatal error: reallocating tvalarray not owning data\n" 
			 << _COREDEBUG_STR << std::endl;
	    ::abort();
#endif
	    x = 0;
	    n = 0;
	}
    }
};

/**
 * Array wrapper class that uses the specified low level array instead of
 * allocating data.
 *
 * @version 0.5, 02/11/2007
 * @author Peter Csizmadia (2007)
 */
template <class T>
class tvalarray_wrapper: public tvalarray<T>
{
public:
    /** Create an array using the specified low level array. */
    tvalarray_wrapper(T* array, unsigned n): tvalarray<T>(array, n, false) {
    }
};

template <class T> tvalarray<T> tvalarray<T>::shift(int k) const
{
    tvalarray<T> result(n);
    T* q = result.x;
    if((k > 0 && (unsigned)k >= n) || (k < 0 && (unsigned)(-k) >= n)) {
	T* end = result.x + n;
	while(q != end) {
	    *(q++) = 0;
	}
    } else if(k >= 0) {
	T* q = result.x;
	T* end = x + n;
	for(T* p = x + k; p != end; ++p) {
	    *(q++) = *p;
	}
	for(int i = 0; i < k; ++i) {
	    *(q++) = 0;
	}
    } else {
	k = -k;
	for(int i = 0; i < k; ++i) {
	    *(q++) = 0;
	}
	T* end = x + n - k;
	for(T* p = x; p != end; ++p) {
	    *(q++) = *p;
	}
    }
    return result;
}

template <class T> tvalarray<T> tvalarray<T>::cshift(int k) const
{
    if(k < 0) {
	k = n - (-k) % n;
    } else if((unsigned)k >= n) {
	k = k % n;
    }
    tvalarray<T> result(n);
    T* end = x + n;
    T* q = result.x;
    for(T* p = x + k; p != end; ++p) {
	*(q++) = *p;
    }
    end = x + k;
    for(T* p = x; p != end; ++p) {
	*(q++) = *p;
    }
    return result;
}

template <class T> tvalarray<T>::tvalarray(const T& value, unsigned _n)
 : ownsData(true)
{
    if(_n > 0) {
	n = _n;
	x = new T[n];
	*this = value; // initialize elements
    } else {
	n=0;
	x = 0;
    }
}

template <class T> tvalarray<T>::tvalarray(T* array, unsigned _n, bool realloc)
 : ownsData(realloc)
{
    if(_n > 0) {
	n=_n;
	if(realloc) {
	    x = new T[n];
	    for(unsigned i = 0; i < n; ++i) {
		x[i] = array[i]; // initialize elements
	    }
	} else {
	    x = array;
	}
    } else {
	n=0;
	x = 0;
    }
}

template <class T> tvalarray<T>::tvalarray(const T* values, unsigned _n)
 : ownsData(true)
{
    if(_n > 0) {
	n=_n;
	x = new T[n];
	for(unsigned i = 0; i < n; ++i) {
	    x[i] = values[i]; // initialize elements
	}
    } else {
	n=0;
	x = 0;
    }
}

/**
 * Copy (also realloc if necessary). The instance for T!=tvalarray<>.
 */
template <class T>
void copy(const tvalarray<T>& source, tvalarray<T>& dest)
{
    if ( dest.n!=source.n )
    {
	dest.realloc(source.n);
    }
    const T* y = source.x;
    for(unsigned i = 0; i < dest.n; ++i)
    {
	dest.x[i] = y[i]; // <- Here T==tvalarray<something> would not work.
    }
}

/**
 * Copy (also realloc if necessary). The instance for T==tvalarray<>.
 */
template <class T>
void copy(const tvalarray< tvalarray<T> >& source, tvalarray< tvalarray<T> >& dest)
{
    if ( dest.n!=source.n )
    {
	dest.realloc(source.n);
    }
    const tvalarray<T>* y = source.x;
    for(unsigned i = 0; i < dest.n; ++i) {
	copy(y[i], dest.x[i]); // <- Here T==tvalarray<something> would not work.
    }
}

/** Resize array and initialize all elements in arrays. */
template <class T>
void resize(tvalarray< T >& arr, unsigned _n, T value/*=T()*/)
{
    if(_n > arr.n) {
	arr.realloc(_n);
    }
    if(arr.n > 0) {
	arr=value;
    }
}

/** Resize array and initialize all elements in hyper-arrays. */
template <class T>
void resize(tvalarray< tvalarray<T> >& arr, unsigned _n, tvalarray<T> value=tvalarray<T>())
{
    if(_n > arr.n) {
	arr.realloc(_n);
    }
    if(arr.n > 0) {
	for ( unsigned int i=0 ; i<arr.n ; ++i ) copy(value, arr[i]);
    }
}

#ifdef DEBUG
template <class T>
inline std::ostream& operator <<(std::ostream& out, const tvalarray<T>& x)
{
    for(unsigned i = 0; i < x.size(); ++i) {
	out << x[i];
	if(i+1 < x.size()) {
	    out << ", ";
	}
    }
    return out;
}
#endif /* DEBUG */

#endif /* _tvalarray_h */

#ifndef Spline3_h
#define Spline3_h

/**
 * Cubic spline interpolation and integration.
 * @version 09/16/2009
 * @since 05/13/1992
 * @author Peter Csizmadia
 */
template <class T>
class Spline3 {

private:
    int N;
    const T* xArray;
    const T* yArray;
    T* aArray;
    T* bArray;
    T* cArray;
    T* dArray;

public:
    Spline3(const T* xarr, const T* yarr, int n) {
	N = n;
	xArray = xarr;
	yArray = yarr;
	init();
    }

    ~Spline3() {
	delete[] aArray;
	delete[] bArray;
	delete[] cArray;
	delete[] dArray;
    }

    T getX(int i) const {
	return xArray[i];
    }

    T interpolate(T x, int i = -1) const {
	if(i < 0) {
	    i = findIndex(x);
	}
	T h = x - xArray[i];
	T a = aArray[i];
	T b = bArray[i];
	T c = cArray[i];
	T d = dArray[i];
	return a + (b + (c + d*h)*h)*h;
    }

    T integrate() const {
	T tot = 0;
	for(unsigned i = 0; i < N - 1; ++i) {
	    T h = xArray[i+1] - xArray[i];
	    T a = aArray[i];
	    T b = bArray[i];
	    T c = cArray[i];
	    T d = dArray[i];
	    tot += (a + (b/2 + (c/3 + d*h/4)*h)*h)*h;
	}
	return tot;
    }

    T x_resampled(int n, int i) const {
	return xArray[0] + (xArray[N - 1] - xArray[0])*i/(n - 1);
    }

    void resample(T* outdata, int n) const {
	int j = 0;
	for(int i = 0; i < N - 1; ++i) {
	    T a = aArray[i];
	    T b = bArray[i];
	    T c = cArray[i];
	    T d = dArray[i];
	    T x;
	    while((x = x_resampled(n, j)) < xArray[i + 1] && j < n) {
		T h = x - xArray[i];
		outdata[j++] = a + (b + (c + d*h)*h)*h;
	    }
	}
	T a = aArray[N - 2];
	T b = bArray[N - 2];
	T c = cArray[N - 2];
	T d = dArray[N - 2];
	while(j < n) {
	    T x = x_resampled(n, j);
	    T h = x - xArray[N - 2];
	    outdata[j++] = a + (b + (c + d*h)*h)*h;
	}
    }

private:
    int findIndex(T x) const {
	int i;
	if(x < xArray[0]) {
	    return -1;
	} else {
	    for(i = 0; i < N - 1 && (x < xArray[i] || x >= xArray[i+1]); ++i);
	    if(i == N - 1) {
		if(x > xArray[N-1]) {
		    return -1;
		}
		--i;
	    }
	}
	return i;
    }

    void init() {
	aArray = new T[N];
	bArray = new T[N];
	cArray = new T[N];
	dArray = new T[N];
	T *alpha = new T[N-1], *beta = new T[N-1];
	T *h = new T[N-1],     *H = new T[N-1];
	T a,b,c,d;

	alpha[0] = 0;
	beta[0] = 0;
	h[0] = xArray[1] - xArray[0];
	H[0] = yArray[1] - yArray[0];
	for(int i = 0; i < N-2; ++i) {
	    H[i+1] = yArray[i+2] - yArray[i+1];
	    a = h[i];
	    b = h[i+1] = xArray[i+2] - xArray[i+1];
	    c =-2*(a+b);
	    d = c-a*alpha[i];
	    alpha[i+1] = b/d;
	    beta[i+1] = (a*beta[i] - 3*(H[i+1]/b-H[i]/a)) / d;
	}
	for(int i = N-2; i>=0; --i) {
	    aArray[i] = yArray[i];
	}
	cArray[0] = 0;
	cArray[N-2] = beta[N-2];
	for(int i = N-3; i > 0; --i) {
	    cArray[i] = alpha[i]*cArray[i+1] + beta[i];
	}
	for(int i = 0; i < N-2; ++i) {
	    dArray[i] = (cArray[i+1]-cArray[i]) / (3*h[i]);
	}
	dArray[N-2] = -cArray[N-2] / (3*h[N-2]);
	for(int i=0; i<N-2; ++i) {
	    bArray[i] = H[i]/h[i] - h[i]*(cArray[i+1]+2*cArray[i])/3;
	}
	bArray[N-2] = H[N-2]/h[N-2] - 2/3.0*h[N-2]*cArray[N-2];
	delete[] alpha;  delete[] h;
	delete[] beta;   delete[] H;
    }
};

template <class T>
inline void spline3resample(const T* xarr1, const T* yarr1, int n1,
			    T* yarr2, int n2)
{
    Spline3<T> spline(xarr1, yarr1, n1);
    spline.resample(yarr2, n2);
}

#endif /* Spline3_H */

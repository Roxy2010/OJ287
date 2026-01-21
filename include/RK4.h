#ifndef RK4_h
#define RK4_h

#include "tvalarray.h"
#include <string>
#include <iostream>

/**
 * Ordinary differential equation interface.
 *
 * @version 0.5, 08/31/2009
 * @since   GridRipper 0.4, 08/02/2006
 * @author  Peter Csizmadia
 */
template <class T>
class ODE
{
public:
    enum PrintWhat { SUMS, COMPONENTS };
private:
    int numComponents;

    const std::string* componentNames;

protected:
    /**
     * Constructs the ODE.
     * @param nc         number of field components
     * @param compnames  array of component names
     */
    ODE(int nc, const std::string* compnames) {
	numComponents = nc;
	componentNames = compnames;
    }

public: 
    virtual ~ODE() { }

    /**
     * Evaluates the ODE.
     * @param f       vector of field values
     * @param offset  offset in field vector
     * @param x       the x coordinate
     * @param df      output vector for the derivatives
     */
    virtual void eval(const T* f, int offset, T x, T* eval_rx) =0;

    /**
     * Gets the number of field components.
     * @return the number of components
     */
    int getNumComponents() const {
	return numComponents;
    }

    /**
     * Gets the name of a component.
     * @param k    the component index
     * @return the component name
     */
    std::string getComponentName(int k) const {
	return componentNames[k];
    }

    /**
     * Sets a parameter or a variable.
     * @param name   parameter name
     * @param value  parameter value
     * @param vars   array of variables.
     * @return true in case of successful setting,
     *         false if no parameter exists with the specified name
     */
    virtual bool set(const std::string& name, T value, tvalarray<T>& vars);
};

/**
 * Fourth order Runge-Kutta method for ODE integration.
 * @version 0.5, 08/28/2009
 * @since  GridRipper 0.4, 08/02/2006
 * @author Peter Csizmadia
 */
template <class T>
class RK4
{
    tvalarray<T> tmpf1;
    tvalarray<T> tmpf2;
    tvalarray<T> tmpf3;
    tvalarray<T> tmpf4;
    ODE<T>* ode;

public:
    RK4(ODE<T>& ode);

    void integrate(tvalarray<T>& f, T x, T dx);

private:
    void integrate(T x, T dx, T* f0, T* f1, T* f2, T* f3, T* f4);

    void step(T* f0, T* fin, T x, T* fout, T dx);

    void step4(T* f0, T* f1, T* f2, T* f3, T x, T* f4, T dx);
};

template <class T>
bool ODE<T>::set(const std::string& name, T value, tvalarray<T>& vars)
{
    for(int i = 0; i < getNumComponents(); ++i) {
	if(getComponentName(i) == name) {
	    vars[i] = value;
	    return true;
	}
    }
    return false;
}

template <class T>
RK4<T>::RK4(ODE<T>& ode): ode(&ode)
{
    int nc = ode.getNumComponents();
    tmpf1.resize(nc);
    tmpf2.resize(nc);
    tmpf3.resize(nc);
    tmpf4.resize(nc);
}

template <class T>
void RK4<T>::integrate(tvalarray<T>& f, T x, T dx)
{
    integrate(x, dx, f, tmpf1, tmpf2, tmpf3, tmpf4);
    for(unsigned i = 0; i < f.size(); ++i) {
	f[i] = tmpf4[i];
    }
}

template <class T>
void RK4<T>::integrate(T x, T dx, T* f0, T* f1, T* f2, T* f3, T* f4)
{
    // Phi(l+1/4) = Phi(l) + Phi'(l)*dx/2
    step(f0, f0, x, f1, 0.5*dx);

    // Phi(l+2/4) = Phi(l) + Phi'(l+1/4)*dx/2
    step(f0, f1, x + 0.5*dx, f2, 0.5*dx);

    // Phi(l+3/4) = Phi(l) + Phi'(l+2/4)*dx
    step(f0, f2, x + 0.5*dx, f3, dx);

    // Phi(l+1) = 1/6*(2*Phi(l+1/4) + 4*Phi(l+2/4) + 2*Phi(l+3/4)
    //                 - 2*Phi(l) + Phi'(l+3/4)*dx)
    step4(f0, f1, f2, f3, x + dx, f4, dx);
}

template <class T>
void RK4<T>::step(T* f0, T* fin, T x, T* fout, T dx)
{
    // fout = f0 + fin'*dt;
    ode->eval(fin, 0, x, fout);
    for(unsigned i = 0; i < tmpf1.size(); ++i) {
	fout[i] = f0[i] + fout[i]*dx;
    }
}

template <class T>
void RK4<T>::step4(T* f0, T* f1, T* f2, T* f3, T x, T* f4, T dx)
{
    tvalarray<T> w(f4, tmpf1.size());
    ode->eval(f3, 0, x, w);

    // f4 = (2*f1 + 4*f2 + 2*f3 + w*dt - 2*f0)/6;
    for(unsigned i = 0; i < w.size(); ++i) {
	f4[i] = (-2*f0[i] + 2*f1[i] + 4*f2[i] + 2*f3[i] + w[i]*dx)/6;
    }
}

#endif /* RK4_h */

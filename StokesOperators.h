/**
* @file StokesOperators.h
* 
* @brief Defines the boundary integral operators used in solving the Stokes Mobility problem.
* 
* The goal here is to first implement the single layer operator. It is defined by the integral
* @f{equation}{
S_j[\rho](\boldsymbol{x}) = \sum_j\int_\Gamma G_{i,j}(\boldsymbol{x} - \boldsymbol{y}) \rho_j(\boldsymbol{y}) d\boldsymbol{y},
 @f}
 where @f$ G_{i,j}(\boldsymbol{r}) @f$ is the Stokeslet tensor:
 @f{equation}{
G_{i , j}(\boldsymbol{r}) = \frac{1}{8 \pi}\left(\frac{\delta_{i , j}}{|\boldsymbol{r}|} + \frac{r_i r_j}{|\boldsymbol{r}|^3}\right).
 @f}

* Given a SphericalVectorField, one can approximate the boundary integral over the sphere's surface
* via the quadrature mentioned SphericalCalc.h . This method is implemented by the SingleLayerDirect class. 
* However it is not accurate for points near the sphere's surface, so we provide an alternate method of evaluation, the SingleLayerHarmonic class which
* evaluates the operator acting on a VSHSeries for the function. 
* 
* We indroduce the StokesPressure operator as an intermediary to the StokesTractionHarmonic operator. The pressure is computed the same way as
* SingleLayerHarmonic is.
* 
* Finally we have the traction operator. It is given by the integral
* @f{equation}{
T_i[\boldsymbol{\rho}](\boldsymbol{x}) = \int_\Gamma T_{i,j,k}(\boldsymbol{x} - \boldsymbol{y})n_k(\boldsymbol{x})\rho_j(\boldsymbol{y}) dy,
 @f}
 *Where @f$ T(\boldsymbol{r}) @f$ is the third rank tensor
 * @f{equation}{
\boldsymbol{T}(\boldsymbol{r}) = -\frac{3}{4 \pi} \frac{r_i r_j r_k}{|\boldsymbol{r}|^5}.
 @f}
 * As with the single layer operator, we provide two methods of evaluating the traction. The first, StokesTractionDirect, 
 * evaluates using the usual quadrature and is valid for points away from the surface of the source sphere. For surfaces near or on the sphere,
 * we Evaluate using StokesTractionHarmonic, which takes in SingleLayerHarmonic and StokesPressure objects and computes the matrix product
 * @f{equation}{
\boldsymbol{t} = \boldsymbol{\sigma} \cdot \boldsymbol{n} = \left[ - pI + \nabla \boldsymbol{u} + \nabla\boldsymbol{u}^T\right] \cdot \boldsymbol{n}  
 @f}
 * exactly, up to the VSHSeries truncation and fourier coefficient accuracy.
*/


#pragma once
#include "Matrix.h"
#include "SphericalHarmonics.h"
#include <cmath>










/**
 * @brief Functor for a term in the Single Layer integral operator evaluated using a VSHSeries representation.
*/
class SingleLayerHarmonicTerm : public SphericalVectorField
{
private:

    int m_;
    int n_;

    Legendre* poly;

    double rhohatVr;
    double rhohatVi;
    double rhohatWr;
    double rhohatWi;
    double rhohatXr;
    double rhohatXi;

    VReal vr;
    VImag vi;
    WReal wr;
    WImag wi;
    XReal xr;
    XImag xi;

    /**
 * @brief Coefficient of VReal or VImag in the SingleLayerHarmonic series.
 * @param r Radius from a SphereCoord.
 * @param n Degree of the VReal or VImag term.
 * @param rhohatV Normalized Fourier coefficient of the associated VReal or VImag function.
 * @param rhohatW Normalized Fourier coefficient of the associated WReal or WImag function.
 * @return @f{equation}{
f_{n,m}^V(r) = \frac{n}{(2n + 1)(2n+3)}r^{-(n+2)} \hat{\rho}_{n,m}^V + \frac{n+1}{4n + 2}\left(r^{-n+2} - r^{-n}\right)\hat{\rho}_{n,m}^W
 @f}
*/
    inline double fV(double r, int n, double rhohatV, double rhohatW) const
    {

        if (n == 0)
            return 0;

        double m = double(n);

        double temp1 = m / ((2.0 * m + 1.0) * (2.0 * m + 3.0) * pow(r, n + 2));
        temp1 *= rhohatV;

        double temp2 = m / (4.0 * m + 2.0) * (1.0 / pow(r, n + 2) - 1.0 / pow(r, n));
        temp2 *= rhohatW;

        return temp1 + temp2;
    }
    /**
     * @brief Computes the derivative of fV() with respect to r. Used in helper functions for SingleLayerHarmonic.
     * @param r Radius from a SphereCoord.
     * @param n Degree of the VReal or VImag term.
     * @param rhohatV Normalized Fourier coefficient of the associated VReal or VImag function.
     * @param rhohatW Normalized Fourier coefficient of the associated WReal or WImag function.
     * @return @f{equation}{
    f_{n,m}^{V'}(r) = \frac{-n(n+2)}{(2n + 1)(2n+3)}r^{-(n+3)} \hat{\rho}_{n,m}^V + \frac{n+1}{4n + 2}\left(n r^{-(n+1) - (n+2)r^{-(n+3)}}\right)\hat{\rho}_{n,m}^W
     @f}
    */
    inline double fVprime(double r, int n, double rhohatV, double rhohatW) const
    {

        if (n == 0)
            return 0;

        double m = double(n);

        double temp1 = -m * (m + 2.0) / ((2.0 * m + 1.0) * (2.0 * m + 3.0) * pow(r, n + 3));
        temp1 *= rhohatV;

        double temp2 = m / (4.0 * m + 2.0) * (-(m + 2.0) / pow(r, n + 3) + m / pow(r, n + 1));
        temp2 *= rhohatW;

        return temp1 + temp2;
    }
    /**
     * @brief Coefficient of WReal or WImag in the SingleLayerHarmonic series.
     * @param r Radius from a SphereCoord.
     * @param n Degree of the WReal or WImag term.
     * @param rhohatW Normalized Fourier coefficient of the associated WReal or WImag function.
     * @return @f{equation}{
     f_{n,m}^W(r) =\frac{n}{(2n + 1)(2n - 1)} r^{-n}\hat{\rho}_{n,m}^W
     @f}
    */
    inline double fW(double r, int n, double rhohatW) const
    {
        if (n == 0)
            return 0;

        double m = double(n);

        double temp1 = (m + 1.0) / ((2.0 * m + 1.0) * (2.0 * m - 1.0) * pow(r, n));
        temp1 *= rhohatW;





        return temp1;

    }
    /**
     * @brief Coefficient of WReal or WImag in the SingleLayerHarmonic series.
     * @param r Radius from a SphereCoord.
     * @param n Degree of the WReal or WImag term.
     * @param rhohatW Normlized Fourier coefficient of the associated WReal or WImag function.
     * @return @f{equation}{
     f_{n,m}^{W'}(r) =\frac{-n^2}{(2n + 1)(2n - 1)} r^{-(n+1)}\hat{\rho}_{n,m}^W
     @f}
    */
    inline double fWprime(double r, int n, double rhohatW) const
    {
        if (n == 0)
            return 0;
        double m = double(n);

        double temp1 = -(m + 1.0) * m / ((2.0 * m + 1.0) * (2.0 * m - 1.0) * pow(r, n + 1));
        temp1 *= rhohatW;





        return temp1;

    }
    /**
     * @brief Coefficient of XReal or XImag in the SingleLayerHarmonic series.
     * @param r Radius from a SphereCoord.
     * @param n Degree of the XReal or XImag term.
     * @param rhohatW Normalized Fourier coefficient of the associated XReal or XImag function.
     * @return @f{equation}{
     f_{n,m}^X(r) =\frac{1}{2n + 1} r^{-(n+1)}\hat{\rho}_{n,m}^X
     @f}
    */
    inline double fX(double r, int n, double rhohatX) const
    {
        if (n == 0)
            return 0;

        double m = double(n);

        double temp = 1.0 / ((2.0 * m + 1.0) * pow(r, n + 1));
        temp *= rhohatX;

        return temp;
    }
    /**
     * @brief Coefficient of XReal or XImag in the SingleLayerHarmonic series.
     * @param r Radius from a SphereCoord.
     * @param n Degree of the XReal or XImag term.
     * @param rhohatW Normalized Fourier coefficient of the associated XReal or XImag function.
     * @return @f{equation}{
     f_{n,m}^{X'}(r) =\frac{-(n+1)}{2n + 1} r^{-(n+2)}\hat{\rho}_{n,m}^X
     @f}
    */
    inline double fXprime(double r, int n, double rhohatX) const
    {

        if (n == 0)
            return 0;

        double m = double(n);

        double temp = -(m + 1.0) / ((2.0 * m + 1.0) * pow(r, n + 2));
        temp *= rhohatX;

        return temp;
    }
    


public:

    SingleLayerHarmonicTerm(Legendre* p = nullptr): poly(p) , vr(p) , vi(p) ,
                                           wr(p) , wi(p),
                                           xr(p) , xi(p){}
    /**
     * @brief Constructs the term in the series.
     * @param series VSHSeries corresponding to data to evaluate. 
     * @param m The order of the term.
     * @param n The degree of the term. 
    */
    SingleLayerHarmonicTerm(const VSHSeries& series,  int m, int n):
        poly(series.P),
        vr(series.P ,m,n),
        vi(series.P,m,n),
        wr(series.P,m,n),
        wi(series.P,m,n),
        xr(series.P,m,n),
        xi(series.P,m,n)
    {
        m_ = m;
        n_ = n;

        rhohatVr = series.terms[n][m].rhohatVr;
        rhohatVi = series.terms[n][m].rhohatVi;
        rhohatWr = series.terms[n][m].rhohatWr;
        rhohatWi = series.terms[n][m].rhohatWi;
        rhohatXr = series.terms[n][m].rhohatXr;
        rhohatXi = series.terms[n][m].rhohatXi;

    }


    SingleLayerHarmonicTerm(const VSHSeries& series, int m, int n , Legendre* p) :
        poly(p),
        vr(p, m, n),
        vi(p, m, n),
        wr(p, m, n),
        wi(p, m, n),
        xr(p, m, n),
        xi(p, m, n)
    {
        m_ = m;
        n_ = n;

        rhohatVr = series.terms[n][m].rhohatVr;
        rhohatVi = series.terms[n][m].rhohatVi;
        rhohatWr = series.terms[n][m].rhohatWr;
        rhohatWi = series.terms[n][m].rhohatWi;
        rhohatXr = series.terms[n][m].rhohatXr;
        rhohatXi = series.terms[n][m].rhohatXi;

    }

    SingleLayerHarmonicTerm(const SingleLayerHarmonicTerm & term):
        vr( term.vr),
        vi( term.vi),
        wr( term.wr),
        wi( term.wi),
        xr( term.xr),
        xi( term.xi)
    {
        
        int m = m_ = term.m_;
        int n = n_ = term.n_;

        rhohatVr = term.rhohatVr;
        rhohatVi = term.rhohatVi;
        rhohatWr = term.rhohatWr;
        rhohatWi = term.rhohatWi;
        rhohatXr = term.rhohatXr;
        rhohatXi = term.rhohatXi;

    }
    SingleLayerHarmonicTerm(const SingleLayerHarmonicTerm& term , Legendre * p) :
        vr(p , term.m_ , term.n_),
        vi(p , term.m_ , term.n_),
        wr(p , term.m_, term.n_),
        wi(p , term.m_ , term.n_),
        xr(p, term.m_, term.n_),
        xi(p , term.m_ , term.n_)
    {

        int m = m_ = term.m_;
        int n = n_ = term.n_;

        rhohatVr = term.rhohatVr;
        rhohatVi = term.rhohatVi;
        rhohatWr = term.rhohatWr;
        rhohatWi = term.rhohatWi;
        rhohatXr = term.rhohatXr;
        rhohatXi = term.rhohatXi;

    }



    
    SingleLayerHarmonicTerm& operator =(const SingleLayerHarmonicTerm& term)
    {

        poly = term.poly;
        int m = m_ = term.m_;
        int n = n_ = term.n_;
        vr = VReal(poly ,  m, n);
        vi = VImag(poly, m, n);
        wr = WReal(poly, m, n);
        wi = WImag(poly, m, n);
        xr = XReal(poly, m, n);
        xi = XImag(poly, m, n);

        rhohatVr = term.rhohatVr;
        rhohatVi = term.rhohatVi;
        rhohatWr = term.rhohatWr;
        rhohatWi = term.rhohatWi;
        rhohatXr = term.rhohatXr;
        rhohatXi = term.rhohatXi;

        return *this;

    }

    

    /**
     * @brief Evaluates the term.
     * @return @f{equation}{
     f_{n,m}^{V}(r)\boldsymbol{V}_n^m(\theta , \phi)
     +  f_{n,m}^{W}(r) \boldsymbol{W}_n^m(\theta , \phi) 
    +  f_{n,m}^{X}(r) \boldsymbol{X}_n^m(\theta , \phi) 
    @f}

     see fV , fW, fX , VReal , VImag, WReal, WImag, XReal, XImag.
    */
    RectCoord  operator()(const SphereCoord & s) const
    {
        RectCoord temp = fV(s.r, n_, rhohatVr, rhohatWr) * vr(s) + fV(s.r, n_, rhohatVi, rhohatWi) * vi(s);
        temp = temp + fW(s.r, n_, rhohatWr) * wr(s) + fW(s.r, n_, rhohatWi) * wi(s);
        temp = temp + fX(s.r, n_, rhohatXr) * xr(s) + fX(s.r, n_, rhohatXi) * xi(s);

        return temp;
    }
    /**
    * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of the Single Layer operator.
    * @returns The expression in operator()(), with @f$ \boldsymbol{e}_r @f$ distributed to all terms.
    * See e.g. VReal::eR().
    */
    double eR(const SphereCoord & s) const
    {
        double temp = fV(s.r, n_, rhohatVr, rhohatWr) * vr.eR(s) + fV(s.r, n_, rhohatVi, rhohatWi) * vi.eR(s);
        temp = temp + fW(s.r, n_, rhohatWr) * wr.eR(s) + fW(s.r, n_, rhohatWi) * wi.eR(s);
        temp = temp + fX(s.r, n_, rhohatXr) * xr.eR(s) + fX(s.r, n_, rhohatXi) * xi.eR(s);

        return temp;
    }
    /**
    * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of the Single Layer operator.
    * @returns The expression in operator()(), with @f$ \boldsymbol{e}_\theta @f$ distributed to all terms.
    * See e.g. VReal::eTheta().
    */
    double eTheta(const SphereCoord & s) const
    {
        double temp = fV(s.r, n_, rhohatVr, rhohatWr) * vr.eTheta(s) + fV(s.r, n_, rhohatVi, rhohatWi) * vi.eTheta(s);
        temp = temp + fW(s.r, n_, rhohatWr) * wr.eTheta(s) + fW(s.r, n_, rhohatWi) * wi.eTheta(s);
        temp = temp + fX(s.r, n_, rhohatXr) * xr.eTheta(s) + fX(s.r, n_, rhohatXi) * xi.eTheta(s);

        return temp;
    }

    /**
    * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of the Single Layer operator.
    * @returns The expression in operator()(), with @f$ \boldsymbol{e}_\phi @f$ distributed to all terms.
    * See e.g. VReal::ePhi().
    */
    double ePhi(const SphereCoord & s) const
    {
        double temp = fV(s.r, n_, rhohatVr, rhohatWr) * vr.ePhi(s) + fV(s.r, n_, rhohatVi, rhohatWi) * vi.ePhi(s);
        temp = temp + fW(s.r, n_, rhohatWr) * wr.ePhi(s) + fW(s.r, n_, rhohatWi) * wi.ePhi(s);
        temp = temp + fX(s.r, n_, rhohatXr) * xr.ePhi(s) + fX(s.r, n_, rhohatXi) * xi.ePhi(s);

        return temp;
    }

    /**
   * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of the Single Layer integral operator,then differentiates with respect to @f$ r @f$.
   * @return The expression in operator()(), with @f$ \boldsymbol{e}_r @f$ distributed to all terms, and @f$f @f$s replaced by @f$f' @f$s.
   * See e.g. fVprime() , VReal::eR().
   */
    double eR_dR(const SphereCoord & s) const
    {
        double temp = fVprime(s.r, n_, rhohatVr, rhohatWr) * vr.eR(s) + fVprime(s.r, n_, rhohatVi, rhohatWi) * vi.eR(s);
        temp = temp + fWprime(s.r, n_, rhohatWr) * wr.eR(s) + fWprime(s.r, n_, rhohatWi) * wi.eR(s);
        temp = temp + fXprime(s.r, n_, rhohatXr) * xr.eR(s) + fXprime(s.r, n_, rhohatXi) * xi.eR(s);

        return temp;
    }

    /**
   * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of the Single Layer integral operator,then differentiates with respect to @f$ r @f$.
   * @return The expression in operator()(), with @f$ \boldsymbol{e}_\theta @f$ distributed to all terms, and @f$f @f$s replaced by @f$f' @f$s.
   * See e.g. fVprime() , VReal::eTheta().
   */
    double eTheta_dR(const SphereCoord & s) const
    {
        double temp = fVprime(s.r, n_, rhohatVr, rhohatWr) * vr.eTheta(s) + fVprime(s.r, n_, rhohatVi, rhohatWi) * vi.eTheta(s);
        temp = temp + fWprime(s.r, n_, rhohatWr) * wr.eTheta(s) + fWprime(s.r, n_, rhohatWi) * wi.eTheta(s);
        temp = temp + fXprime(s.r, n_, rhohatXr) * xr.eTheta(s) + fXprime(s.r, n_, rhohatXi) * xi.eTheta(s);

        return temp;
    }

    /**
   * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of the Single Layer integral operator,then differentiates with respect to @f$ r @f$.
   * @return The expression in operator()(), with @f$ \boldsymbol{e}_\phi @f$ distributed to all terms, and @f$f @f$s replaced by @f$f' @f$s.
   * See e.g. fVprime() , VReal::ePhi().
   */
    double ePhi_dR(const SphereCoord & s) const
    {
        double temp = fVprime(s.r, n_, rhohatVr, rhohatWr) * vr.ePhi(s) + fVprime(s.r, n_, rhohatVi, rhohatWi) * vi.ePhi(s);
        temp = temp + fWprime(s.r, n_, rhohatWr) * wr.ePhi(s) + fWprime(s.r, n_, rhohatWi) * wi.ePhi(s);
        temp = temp + fXprime(s.r, n_, rhohatXr) * xr.ePhi(s) + fXprime(s.r, n_, rhohatXi) * xi.ePhi(s);

        return temp;
    }

    /**
   * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of the Single Layer integral operator,then differentiates with respect to @f$ \theta @f$.
   * @return The expression in operator()(), e.g VReal::eR_dTheta() distributed.
   */
    double eR_dTheta(const SphereCoord & s) const
    {
        double temp = fV(s.r, n_, rhohatVr, rhohatWr) * vr.eR_dTheta(s) + fV(s.r, n_, rhohatVi, rhohatWi) * vi.eR_dTheta(s);
        temp = temp + fW(s.r, n_, rhohatWr) * wr.eR_dTheta(s) + fW(s.r, n_, rhohatWi) * wi.eR_dTheta(s);
        temp = temp + fX(s.r, n_, rhohatXr) * xr.eR_dTheta(s) + fX(s.r, n_, rhohatXi) * xi.eR_dTheta(s);

        return temp;
    }

    /**
   * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of the Single Layer integral operator,then differentiates with respect to @f$ \theta @f$.
   * @return The expression in operator()(), e.g VReal::eTheta_dTheta() distributed.
   */
    double eTheta_dTheta(const SphereCoord & s) const
    {
        double temp = fV(s.r, n_, rhohatVr, rhohatWr) * vr.eTheta_dTheta(s) + fV(s.r, n_, rhohatVi, rhohatWi) * vi.eTheta_dTheta(s);
        temp = temp + fW(s.r, n_, rhohatWr) * wr.eTheta_dTheta(s) + fW(s.r, n_, rhohatWi) * wi.eTheta_dTheta(s);
        temp = temp + fX(s.r, n_, rhohatXr) * xr.eTheta_dTheta(s) + fX(s.r, n_, rhohatXi) * xi.eTheta_dTheta(s);

        return temp;
    }

    /**
   * @brief Extracts the @f$ \boldsymbol{e}_\Phi @f$ component of the Single Layer integral operator,then differentiates with respect to @f$ \theta @f$.
   * @return The expression in operator()(), e.g VReal::ePhi_dTheta() distributed.
   */
    double ePhi_dTheta(const SphereCoord & s) const
    {
        double temp = fV(s.r, n_, rhohatVr, rhohatWr) * vr.ePhi_dTheta(s) + fV(s.r, n_, rhohatVi, rhohatWi) * vi.ePhi_dTheta(s);
        temp = temp + fW(s.r, n_, rhohatWr) * wr.ePhi_dTheta(s) + fW(s.r, n_, rhohatWi) * wi.ePhi_dTheta(s);
        temp = temp + fX(s.r, n_, rhohatXr) * xr.ePhi_dTheta(s) + fX(s.r, n_, rhohatXi) * xi.ePhi_dTheta(s);

        return temp;
    }

    /**
   * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of the Single Layer integral operator,then differentiates with respect to @f$ \phi @f$.
   * @return The expression in operator()(), e.g VReal::eR_dPhi() distributed.
   */
    double eR_dPhi(const SphereCoord & s) const
    {
        double temp = fV(s.r, n_, rhohatVr, rhohatWr) * vr.eR_dPhi(s) + fV(s.r, n_, rhohatVi, rhohatWi) * vi.eR_dPhi(s);
        temp = temp + fW(s.r, n_, rhohatWr) * wr.eR_dPhi(s) + fW(s.r, n_, rhohatWi) * wi.eR_dPhi(s);
        temp = temp + fX(s.r, n_, rhohatXr) * xr.eR_dPhi(s) + fX(s.r, n_, rhohatXi) * xi.eR_dPhi(s);

        return temp;
    }

    /**
  * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of the Single Layer integral operator,then differentiates with respect to @f$ \phi @f$.
  * @return The expression in operator()(), e.g VReal::eTheta_dPhi() distributed.
  */
    double eTheta_dPhi(const SphereCoord & s) const
    {
        double temp = fV(s.r, n_, rhohatVr, rhohatWr) * vr.eTheta_dPhi(s) + fV(s.r, n_, rhohatVi, rhohatWi) * vi.eTheta_dPhi(s);
        temp = temp + fW(s.r, n_, rhohatWr) * wr.eTheta_dPhi(s) + fW(s.r, n_, rhohatWi) * wi.eTheta_dPhi(s);
        temp = temp + fX(s.r, n_, rhohatXr) * xr.eTheta_dPhi(s) + fX(s.r, n_, rhohatXi) * xi.eTheta_dPhi(s);

        return temp;
    }

    /**
  * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of the Single Layer integral operator,then differentiates with respect to @f$ \phi @f$.
  * @return The expression in operator()(), e.g VReal::ePhi_dPhi() distributed.
  */
    double ePhi_dPhi(const SphereCoord & s) const
    {
        double temp = fV(s.r, n_, rhohatVr, rhohatWr) * vr.ePhi_dPhi(s) + fV(s.r, n_, rhohatVi, rhohatWi) * vi.ePhi_dPhi(s);
        temp = temp + fW(s.r, n_, rhohatWr) * wr.ePhi_dPhi(s) + fW(s.r, n_, rhohatWi) * wi.ePhi_dPhi(s);
        temp = temp + fX(s.r, n_, rhohatXr) * xr.ePhi_dPhi(s) + fX(s.r, n_, rhohatXi) * xi.ePhi_dPhi(s);

        return temp;
    }

    


};

/**
 * @brief Functor for Single Layer operator. Sums SingleLayerHarmonicTerm objects.
*/
class SingleLayerHarmonic : public VectorFieldSum
{
private:

    Legendre* P;
    
    std::vector<std::vector<SingleLayerHarmonicTerm>> terms;

public:
    
   SphereCoordSystem* system; /**< The center of the Spherical Coordinate system of the boundary integral. */

    /**
     * @brief Constructs the Functor.
     * @param n number of the terms in the series.
     * @param c center of the SphereCoord in the integral.
    */
    SingleLayerHarmonic(int n ,Legendre* p, SphereCoordSystem*sys) : P(p)
    {

        system = sys;
        terms.resize(n + 1);

        for (size_t i = 0; i < terms.size(); i++)
        {
            terms[i].resize(i + 1);
            for (int j = 0; j <= i; j++)
                append(&terms[i][j]);
        }
    }

    /**
     * @brief Generates the Functor from a Fourier Series representation.
     * @param series the Source boundary data.
    */
    SingleLayerHarmonic(const VSHSeries& series) : P(series.P)
    {
        terms.resize(series.N + 1);
        for (int i = 0; i <= series.N; i++)
            terms[i].resize( i + 1);
        for (int n = 0; n <= series.N; n++)
            for (int m = 0; m <= n; m++)
            {
                terms[n][m] = SingleLayerHarmonicTerm(series, m, n);
                append(&terms[n][m]);
            }

        system = series.system;

    }
    /**
     * @brief Copy Consructor
     * @param slh Single Layer to copy.
    */
    SingleLayerHarmonic(const SingleLayerHarmonic& slh): P(slh.P)
    {
        terms = slh.terms;
        for (int n = 0; n < terms.size(); n++)
            for (int m = 0; m <= n; m++)
            {
                append(&terms[n][m]);
            }

        system= slh.system;
    }

    SingleLayerHarmonic(const SingleLayerHarmonic& slh , Legendre* p) : P(p)
    {
        terms = slh.terms;
        for (int n = 0; n < terms.size(); n++)
            for (int m = 0; m <= n; m++)
            {
                terms[n][m] = SingleLayerHarmonicTerm(slh.terms[n][m], p);
                append(&terms[n][m]);
            }

        system = slh.system;
    }

   

    /**
     * @brief Assignment operator.
     * @param slh Single Layer to assign.
    */
    SingleLayerHarmonic & operator =(const SingleLayerHarmonic& slh)
    {
        terms = slh.terms;
        system = slh.system;
        clear();

        for (int n = 0; n < terms.size(); n++)
            for (int m = 0; m <= n; m++)
                append(&terms[n][m]);

        return *this;
    }

   


    
    /**
     * @brief Generates the Functor from a Fourier Series representation.
     * @param series the Source boundary data.
    */
    void solve(VSHSeries& series)
    {
        clear();
        terms.clear();
        terms.resize(series.N + 1);
        for (int i = 0; i <= series.N; i++)
            terms[i].resize(i + 1);
        for (int n = 0; n <= series.N; n++)
            for (int m = 0; m <= n; m++)
            {
                terms[n][m] = SingleLayerHarmonicTerm(series, m, n);
                append(&terms[n][m]);
            }

        system = series.system;
    }

    /**
     * @brief Generates the Functor from boundary data. The Series is constructed as an intermediate step.
     * @param series the Source boundary data.
    */
    void solve(SphericalVectorField* rho, int N ,QuadConstants* consts, Legendre* p,SphereCoordSystem* sys)
    {
        VSHSeries series( N,consts,p, sys);
        series.approximate(rho);

        solve(series);
    }

    /**
     * @brief Evaluates the Single Layer operator at a point.
     * @return The sum of terms. See SingleLayerHarmonicTerm::operator()().
    */
    RectCoord  operator()(const SphereCoord & s) const
    {

        SphereCoord temp = recast(s, system);
        int N = terms.size();
        P->populate(cos(temp.s.theta), N, N);
        return VectorFieldSum::operator()(temp);

    }


    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of the Single Layer operator.
     * @return the sum of values of SingleLayerHarmonicTerm::eR().
    */
    double eR(const SphereCoord & s) const
    {
        double temp = 0.0;

        SphereCoord sc = recast(s, system);

        P->populate(cos(temp.s.theta), N + 1, N  + 1);
     
        for(int i = 0; i < terms.size();i++)
            for (int j = 0; j < terms[i].size();j++)
            {
                temp += terms[i][j].eR(sc);
            }

        return temp;
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of the Single Layer operator.
     * @return the sum of values of SingleLayerHarmonicTerm::eTheta().
    */
    double eTheta(const SphereCoord & s) const
    {
        double temp = 0.0;

       P->populate(cos(temp.s.theta), N + 1, N  + 1);
     
        SphereCoord sc = recast(s, system);

        for (int i = 0; i < terms.size(); i++)
            for (int j = 0; j < terms[i].size(); j++)
            {
                temp += terms[i][j].eTheta(sc);
            }

        return temp;
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of the Single Layer operator.
     * @return the sum of values of SingleLayerHarmonicTerm::ePhi().
    */
    double ePhi(const SphereCoord & s) const
    {
        double temp = 0.0;

        SphereCoord sc = recast(s, system);

     P->populate(cos(temp.s.theta), N + 1, N  + 1);

        for (int i = 0; i < terms.size(); i++)
            for (int j = 0; j < terms[i].size(); j++)
            {
                temp += terms[i][j].ePhi(sc);
            }

        return temp;
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of the Single Layer operator, then differentiates with respect to @f$ r @f$.
     * @return the sum of values of SingleLayerHarmonicTerm::eR_dR().
    */
    double eR_dR(const SphereCoord & s) const
    {
        double temp = 0.0;

        SphereCoord sc = recast(s, system);

     P->populate(cos(temp.s.theta), N + 1, N  + 1);

        for (int i = 0; i < terms.size(); i++)
            for (int j = 0; j < terms[i].size(); j++)
            {
                temp += terms[i][j].eR_dR(sc);
            }

        return temp;
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of the Single Layer operator, then differentiates with respect to @f$ r @f$.
     * @return the sum of values of SingleLayerHarmonicTerm::eTheta_dR().
    */
    double eTheta_dR(const SphereCoord & s) const
    {
        double temp = 0.0;

        SphereCoord sc = recast(s, system);

     P->populate(cos(temp.s.theta), N + 1, N  + 1);

        for (int i = 0; i < terms.size(); i++)
            for (int j = 0; j < terms[i].size(); j++)
            {
                temp += terms[i][j].eTheta_dR(sc);
            }

        return temp;
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of the Single Layer operator, then differentiates with respect to @f$ r @f$.
     * @return the sum of values of SingleLayerHarmonicTerm::ePhi_dR().
    */
    double ePhi_dR(const SphereCoord & s) const
    {
        double temp = 0.0;

     P->populate(cos(temp.s.theta), N + 1, N  + 1);

        SphereCoord sc = recast(s, system);

        for (int i = 0; i < terms.size(); i++)
            for (int j = 0; j < terms[i].size(); j++)
            {
                temp += terms[i][j].ePhi_dR(sc);
            }

        return temp;
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of the Single Layer operator, then differentiates with respect to @f$ \theta @f$.
     * @return the sum of values of SingleLayerHarmonicTerm::eR_dTheta().
    */
    double eR_dTheta(const SphereCoord & s) const
    {
        double temp = 0.0;

     P->populate(cos(temp.s.theta), N + 1, N  + 1);

        SphereCoord sc = recast(s, system);

        for (int i = 0; i < terms.size(); i++)
            for (int j = 0; j < terms[i].size(); j++)
            {
                temp += terms[i][j].eR_dTheta(sc);
            }

        return temp;
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of the Single Layer operator, then differentiates with respect to @f$ \theta @f$.
     * @return the sum of values of SingleLayerHarmonicTerm::eTheta_dTheta().
    */
    double eTheta_dTheta(const SphereCoord & s) const
    {
        double temp = 0.0;

     P->populate(cos(temp.s.theta), N + 1, N  + 1);

        SphereCoord sc = recast(s, system);

        for (int i = 0; i < terms.size(); i++)
            for (int j = 0; j < terms[i].size(); j++)
            {
                temp += terms[i][j].eTheta_dTheta(sc);
            }

        return temp;
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of the Single Layer operator, then differentiates with respect to @f$ r @f$.
     * @return the sum of values of SingleLayerHarmonicTerm::eR_dR().
    */
    double ePhi_dTheta(const SphereCoord & s) const
    {
        double temp = 0.0;

     P->populate(cos(temp.s.theta), N + 1, N  + 1);

        SphereCoord sc = recast(s, system);

        for (int i = 0; i < terms.size(); i++)
            for (int j = 0; j < terms[i].size(); j++)
            {
                temp += terms[i][j].ePhi_dTheta(sc);
            }

        return temp;
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of the Single Layer operator, then differentiates with respect to @f$ \phi @f$.
     * @return the sum of values of SingleLayerHarmonicTerm::eR_dPhi().
    */
    double eR_dPhi(const SphereCoord & s) const
    {
        double temp = 0.0;

        SphereCoord sc = recast(s, system);

     P->populate(cos(temp.s.theta), N + 1, N  + 1);

        for (int i = 0; i < terms.size(); i++)
            for (int j = 0; j < terms[i].size(); j++)
            {
                temp += terms[i][j].eR_dPhi(sc);
            }

        return temp;
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of the Single Layer operator, then differentiates with respect to @f$ \phi @f$.
     * @return the sum of values of SingleLayerHarmonicTerm::eTheta_dPhi().
    */
    double eTheta_dPhi(const SphereCoord & s) const
    {
        double temp = 0.0;

     P->populate(cos(temp.s.theta), N + 1, N  + 1);

        SphereCoord sc = recast(s, system);

        for (int i = 0; i < terms.size(); i++)
            for (int j = 0; j < terms[i].size(); j++)
            {
                temp += terms[i][j].eTheta_dPhi(sc);
            }

        return temp;
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of the Single Layer operator, then differentiates with respect to @f$ \phi @f$.
     * @return the sum of values of SingleLayerHarmonicTerm::ePhi_dPhi().
    */
    double ePhi_dPhi(const SphereCoord & s) const
    {
        double temp = 0.0;

     P->populate(cos(temp.s.theta), N + 1, N  + 1);

        SphereCoord sc = recast(s, system);

        for (int i = 0; i < terms.size(); i++)
            for (int j = 0; j < terms[i].size(); j++)
            {
                temp += terms[i][j].ePhi_dPhi(sc);
            }

        return temp;
    }

    
};


/**
 * @brief Term in the Evaluation of StokesPressure.
*/
class StokesPressureTerm : public SphericalScalarFunction
{
private:


    int m_;
    int n_;

    Legendre* poly;
    double rhohatWr;
    double rhohatWi;

    YReal Yr;
    YImag Yi;

    /**
     * @brief Coefficient of WReal or WImag in the StokesPressure series.
     * @param r Radius from a SphereCoord.
     * @param n Degree of the WReal or WImag term.
     * @param rhohatW Normalized Fourier coefficient of the associated WReal or WImag function.
     * @return @f{equation}{
    g_{n,m}(r) = n r^{-(n+1) }\hat{\rho}_{n,m}^W
     @f}
    */
    inline double g(double r, int n, double rhohatW) const
    {
        double m = double(n);

        if (isnan(rhohatW))
            return 0.0;

        double temp = m / pow(r, n + 1);
        temp *= rhohatW;



        return temp;
    }

public:

    StokesPressureTerm(Legendre*p = nullptr) : poly(p) , Yr(p) , Yi(p) {}
    /**
     * @brief Constructs the term.
     * @param series the series to extract information from.
     * @param m The order of the term.
     * @param n The Degree of the term.
    */
    StokesPressureTerm(const VSHSeries& series, int m, int n): poly(series.P),Yr(series.P, m , n) , Yi(series.P, m , n)
    {
        m_ = m;
        n_ = n;
        rhohatWr = series.terms[n][m].rhohatWr;
        rhohatWi = series.terms[n][m].rhohatWi;

    }

    StokesPressureTerm(const VSHSeries& series, int m, int n , Legendre*p) :poly(p), Yr(p, m, n), Yi(p, m, n)
    {
        m_ = m;
        n_ = n;
        rhohatWr = series.terms[n][m].rhohatWr;
        rhohatWi = series.terms[n][m].rhohatWi;

    }

    StokesPressureTerm(const StokesPressureTerm& term) : Yr(term.Yr), Yi(term.Yi)
    {
        m_ = term.m_;
        n_ = term.n_;

        rhohatWr = term.rhohatWr;
        rhohatWi = term.rhohatWi;
    }

    StokesPressureTerm(const StokesPressureTerm& term , Legendre*p) :poly(p), Yr(p, term.m_,term.n_) , Yi(p, term.m_,term.n_)
    {
        m_ = term.m_;
        n_ = term.n_;

        rhohatWr = term.rhohatWr;
        rhohatWi = term.rhohatWi;
    }

    StokesPressureTerm & operator =(const StokesPressureTerm& term)
    {
        m_ = term.m_;
        n_ = term.n_;
        poly = term.poly;

        Yr = YReal(poly, m_, n_);
        Yi = YImag(poly, m_, n_);

        rhohatWr = term.rhohatWr;
        rhohatWi = term.rhohatWi;

        return *this;
    }
    /**
     * @brief Evaluates the term.
     * @return @f$ g_{n,m}(r) \mathfrak{R}(Y_n^m)(\theta , \phi) + g_{n,m}(r) \mathfrak{I}(Y_n^m)(\theta , \phi).@f$
    */
    double operator()(const SphereCoord & s) const
    {
        double temp = g(s.r, n_, rhohatWr) * Yr(s) + g(s.r, n_, rhohatWi) * Yi(s);

        //std::cout << temp << "\n";

        return temp;
    }

};

/**
 * @brief Evaluates the pressure operator by summing StokesPressureTerm objects.
*/
class StokesPressure : public SphereScalFunctionSum
{
private:
    
    Legendre* P;
    std::vector<std::vector<StokesPressureTerm>> terms;

public:

    SphereCoordSystem* system;

    /**
     * @brief Constructs a pressure term.
     * @param n The number of terms in the sum.
     * @param c The center of the source sphere.
    */
    StokesPressure(int n, Legendre* p , SphereCoordSystem* sys) : P(p)
    {

        system = sys;
        terms.resize(n+1);

        for (size_t i = 0; i < terms.size(); i++)
        {
            terms[i].resize(i + 1);
            for (int j = 0; j <= i; j++)
                append(&terms[i][j]);
        }
    }

    /**
     * @brief Constructs pressure from fourier data.
     * @param series Series to extract data from.
     * @param n The name.
    */
    StokesPressure(const VSHSeries& series, std::string n = "") : P(series.P)
    {
        name = n;
        system = series.system;
        terms.resize(series.N + 1);
        for (int i = 0; i <= series.N; i++)
            terms[i].resize(i + 1);

        for (int n = 0; n <= series.N; n++)
            for (int m = 0; m <= n; m++)
            {
                terms[n][m] = StokesPressureTerm(series, m, n);
                append(&terms[n][m]);
            }

    }

    /**
     * @brief Copy Constructor.
     * @param p the pressure to copy.
    */
    StokesPressure(const StokesPressure& p): P(p.P)
    {
        system = p.system;
        

        terms = p.terms;

        for (int n = 0; n < terms.size(); n++)
            for (int m = 0; m <= n; m++)
            {
                append(&terms[n][m]);
            }
    }


    StokesPressure(const StokesPressure& p , Legendre* poly) : P(poly)
    {
        system = p.system;


        terms = p.terms;

        for (int n = 0; n < terms.size(); n++)
            for (int m = 0; m <= n; m++)
            {
                terms[n][m] = StokesPressureTerm(p.terms[n][m], poly);
                append(&terms[n][m]);
            }
    }

    

    /**
     * @brief Assignment operator
     * @param p the pressure to copy
     * @return the 
    */
    const StokesPressure & operator =(const StokesPressure& p)
    {
        terms = p.terms;
        system = p.system;
        clear();

        for (int n = 0; n < terms.size(); n++)
            for (int m = 0; m <= n; m++)
                append(&terms[n][m]);



        return *this;
    }

    

    /**
     * @brief solves for pressure from series data.
    */
    void solve(const VSHSeries& series)
    {
        clear();
        terms.clear();
        terms.resize(series.N + 1);
        for (int i = 0; i <= series.N; i++)
            terms[i].resize(i + 1);

        for (int n = 0; n <= series.N; n++)
            for (int m = 0; m <= n; m++)
            {
                terms[n][m] = StokesPressureTerm(series, m, n , P);
                append(&terms[n][m]);
            }

        system = series.system;
    }

    /**
     * @brief solves for pressure from boundary data. Computes a VSHSeries as intermediary.
     * @param rho the density to solve.
     * @param N maximum degree.
     * @param center center of the spherical coordinate system.
    */
    void solve(SphericalVectorField* rho, QuadConstants* consts, int N,Legendre*p, SphereCoordSystem* sys)
    {
        VSHSeries series(N, consts ,p, sys);
        series.approximate(rho);

        solve(series);
    }

    /**
     * @brief Evaluates the pressure. 
     * @return the sum of StokesPressureTerm::operator()().
    */
    double operator()(const SphereCoord & s) const
    {

        SphereCoord temp = recast(s, system);
        int N = terms.size();
        P->populate(cos(temp.s.theta), N, N);
        return SphereScalFunctionSum::operator()(temp);

    }

};

vec<3> RecttoVec(RectCoord x)
{
    vec<3> temp;
    temp[0] = x.x;
    temp[1] = x.y;
    temp[2] = x.z;
    return temp;
}

RectCoord VectoRect(vec<3> temp)
{
    RectCoord x;
    x.x = temp[0];
    x.y = temp[1];
    x.z = temp[2];
    return x;
}

/**
 * @brief Evaluates the Single Layer operator for points away from the surface, acting on a SphericalVectorField.
*/
class SingleLayerDirect : public SphericalVectorField
{
private:

    SphereCoordSystem* system;
    SphericalVectorField* u;
    QuadConstants* consts;

   Matrix<3, 3> Kernel(const SphereCoord & x, const SphereCoord &y) const
   {
       Matrix<3, 3> temp;

       //convert to the argument of the stokeslet.
       vec<3> r = RecttoVec(x - y);
       double R = sqrt(dot(r,r));

       //iterate over entries over the tensor.
       for(int i = 0; i < 3; i++)
           for (int j = 0; j < 3; j++)
           {
               //formula for the i,jth entry of the stokeslet
               temp[i][j] = (double)(i ==j)/ R + r[i] * r[j] / (R * R * R);
               temp[i][j] /= 8.0 * consts->PI;
           }

       return temp;
   } 
public:

    /**
     * @brief Constructs the functor, assigned to a SphericalVectorField.
     * @param rho The density the opeator is acting on.
     * @param c The center of the sphere to integrate over.
     * @param na the name.
    */
    SingleLayerDirect(SphericalVectorField *rho,QuadConstants* cts, SphereCoordSystem* sys, std::string na = "Direct Single Layer")
    {
        u = rho;
        name = na;
        system = sys;
        consts = cts;
    }

    /**
     * @brief Evaluates the operator.
     * @return 
    */
    RectCoord  operator()(const SphereCoord & x) const
    {

        

        vec<3> total = 0.0;
        //loop over trapezoidal(outside) and legendre(inside) nodes.
        for (int p = 0; p < consts->NUMTRAPNODES; p++)
            for (int i = 0; i < consts->NUMGLNODES; i++)
            {
                //generate point on unit sphere, note theta is the first argument, phi the second.
                SurfaceCoord s(consts->GLnodes[i] , 2.0 * consts->PI * (double)p / (double)consts->NUMTRAPNODES);
                //convert to spherical coordinate.
                SphereCoord y(s , system);
                //add contribution at y to integral.
                total = total +  dot(Kernel(x , y), RecttoVec((*u)(y)))  * consts->GLweights[i] * (2.0 * consts->PI / consts->NUMTRAPNODES);
            }

        return VectoRect(total);
    }
};

/**
 * @brief Evaluates the Single Layer operator for points away from the surface, acting on SphereData.
*/
class SingleLayerDirectDiscrete : public SphericalVectorField
{
private:

    SphereCoordSystem* system;
    SphereData* data;
    QuadConstants* consts;

    Matrix<3, 3> Kernel(const SphereCoord& x, const SphereCoord& y) const
    {
        Matrix<3, 3> temp;

        //convert to the argument of the stokeslet.
        vec<3> r = RecttoVec(x - y);
        double R = sqrt(dot(r, r));

        //iterate over entries over the tensor.
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                //formula for the i,jth entry of the stokeslet
                temp[i][j] = (double)(i == j) / R + r[i] * r[j] / (R * R * R);
                temp[i][j] /= 8.0 * consts->PI;
            }

        return temp;
    }
public:

    /**
     * @brief Constructs the functor, assigned to a SphereData.
     * @param d The density the opeator is acting on.
     * @param c The center of the sphere to integrate over.
     * @param na the name.
    */
    SingleLayerDirectDiscrete(SphereData* d, QuadConstants* cts, SphereCoordSystem* sys, std::string na = "Direct Single Layer Discrete")
    {
        data = d;
        name = na;
        system = sys;
        consts = cts;
    }

    /**
     * @brief Evaluates the operator.
     * @return
    */
    RectCoord  operator()(const SphereCoord& x) const
    {



        vec<3> total = 0.0;
        //loop over trapezoidal(outside) and legendre(inside) nodes.
        for (int p = 0; p < consts->NUMTRAPNODES; p++)
            for (int i = 0; i < consts->NUMGLNODES; i++)
            {
                //generate point on unit sphere, note theta is the first argument, phi the second.
                SurfaceCoord s(consts->GLnodes[i], 2.0 * consts->PI * (double)p / (double)consts->NUMTRAPNODES);
                //convert to spherical coordinate.
                SphereCoord y(s, system);
                //add contribution at y to integral.
                total = total + dot(Kernel(x, y), RecttoVec((*data)[i][p])) * consts->GLweights[i] * (2.0 * consts->PI / consts->NUMTRAPNODES);
            }

        return VectoRect(total);
    }
};


/**
 * @brief Evaluates the traction using quadrature. Valid for points far from the source surface. Acts on a SphericalVectorField.
*/
class StokesTractionDirect : public SphericalVectorField
{
private:

    SphericalVectorField* u;
    QuadConstants* consts;
    SphereCoordSystem* system;

    Matrix<3, 3> Kernel(const SphereCoord & x,const SphereCoord &y) const
    {
        Matrix<3, 3> temp;


        //convert to the argument of the stokeslet.
        vec<3> r = RecttoVec(x - y);
        double R = sqrt(dot(r, r));
        double rn = dot(r,RecttoVec(e_r(x)));
        //iterate over entries over the tensor.
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                //formula for the i,jth entry of the stokeslet
                temp[i][j] = r[i] * r[j] * rn / pow(R , 5);
                temp[i][j] *= - 0.75 / consts->PI;
            }

        return temp;
    }
public:

    /**
     * @brief Contruct the Functor from a density.
     * @param rho the function the density will act on.
     * @param c The center of the spherical coordinate system of the boundary integral.
     * @param na the name.
    */
    StokesTractionDirect(SphericalVectorField* rho,QuadConstants* cts,  SphereCoordSystem*sys, std::string na = "Traction Direct")
    {
        u = rho;
        name = na;
        system = sys;
        consts = cts;
    }

    /**
     * @brief Evaluates the traction at the target point.
    */
    RectCoord operator()(const SphereCoord & x) const
    {



        vec<3> total = 0.0;
        //loop over trapezoidal(outside) and legendre(inside) nodes.
        for (int p = 0; p < consts->NUMTRAPNODES; p++)
            for (int i = 0; i < consts->NUMGLNODES; i++)
            {
                //generate point on unit sphere, note theta is the first argument, phi the second.
                SurfaceCoord s(consts->GLnodes[i], 2.0 * consts->PI * (double)p / (double)consts->NUMTRAPNODES);
                //convert to spherical coordinate.
                SphereCoord y(s , system);
                //add contribution at y to integral.
                total = total + dot(Kernel(x, y), RecttoVec((*u)(y))) * consts->GLweights[i] * (2.0 * consts->PI / consts->NUMTRAPNODES);
                if (isnan(total[0]) || isnan(total[1]) || isnan(total[2]))
                {
                    std::cout << "nan detected from function " << name << " at " << y <<"\n";
                }
            }

        return VectoRect(total);
    }
};


/**
 * @brief Evaluates the traction using quadrature. Valid for points far from the source surface. Acts on SphereData.
*/
class StokesTractionDirectDiscrete : public SphericalVectorField
{
private:

    SphereData* data;
    QuadConstants* consts;
    SphereCoordSystem* system;

    Matrix<3, 3> Kernel(const SphereCoord& x, const SphereCoord& y) const
    {
        Matrix<3, 3> temp;


        //convert to the argument of the stokeslet.
        vec<3> r = RecttoVec(x - y);
        double R = sqrt(dot(r, r));
        double rn = dot(r, RecttoVec(e_r(x)));
        //iterate over entries over the tensor.
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                //formula for the i,jth entry of the stokeslet
                temp[i][j] = r[i] * r[j] * rn / pow(R, 5);
                temp[i][j] *= -0.75 / consts->PI;
            }

        return temp;
    }
public:

    /**
     * @brief Contruct the Functor from a density.
     * @param rho the function the density will act on.
     * @param c The center of the spherical coordinate system of the boundary integral.
     * @param na the name.
    */
    StokesTractionDirectDiscrete(SphereData* d,QuadConstants* cts, SphereCoordSystem* sys, std::string na = "Traction Direct Discrete")
    {
        data = d;
        name = na;
        system = sys;
        consts = cts;
    }

    /**
     * @brief Evaluates the traction at the target point.
    */
    RectCoord operator()(const SphereCoord& x) const
    {



        vec<3> total = 0.0;
        //loop over trapezoidal(outside) and legendre(inside) nodes.
        for (int p = 0; p < consts->NUMTRAPNODES; p++)
            for (int i = 0; i < consts->NUMGLNODES; i++)
            {
                //generate point on unit sphere, note theta is the first argument, phi the second.
                SurfaceCoord s(consts->GLnodes[i], 2.0 * consts->PI * (double)p / (double)consts->NUMTRAPNODES);
                //convert to spherical coordinate.
                SphereCoord y(s, system);
                //add contribution at y to integral.
                total = total + dot(Kernel(x, y), RecttoVec((*data)[i][p])) * consts->GLweights[i] * (2.0 * consts->PI / consts->NUMTRAPNODES);
                if (isnan(total[0]) || isnan(total[1]) || isnan(total[2]))
                {
                    std::cout << "nan detected from function " << name << " at " << y << "\n";
                }
            }

        return VectoRect(total);
    }
};

/**
 * @brief The Evaluation of the traction via Spherical Hermonics.
 * 
 * We evaluate the expression 
 * @f{equation}{
\boldsymbol{t} = \boldsymbol{\sigma} \cdot \boldsymbol{n} = \left[ - pI + \nabla \boldsymbol{u} + \nabla\boldsymbol{u}^T\right] \cdot \boldsymbol{n}  
 @f}
* where @f$ \boldsymbol{u} @f$ is given by a SingleLayerHarmonic object, and @f$ p @f$ is given by a StokesPressure object.
* The entries of the symmetric matrix @f$ \boldsymbol{\sigma} @f$ are below, see  
@link Traction.pdf Traction.pdf @endlink for a derivation.
@f{eqnarray}{
\sigma_{rr} &=&  2\frac{\partial \boldsymbol{u}_r}{\partial r} \\
\sigma_{r\theta} &=& \frac{\partial \boldsymbol{u}_\theta}{\partial r} + \frac{1}{r}\frac{\partial \boldsymbol{u}_r}{\partial \theta} 
                    - \frac{1}{r}\boldsymbol{u}_\theta \\
\sigma_{r\phi} &=& \frac{\partial \boldsymbol{u}_\phi}{\partial r}  + \frac{1}{r \sin\theta}\frac{\partial \boldsymbol{u}_r}{\partial \phi}
                   - \frac{1}{r}\boldsymbol{u}_\phi \\
\sigma_{\theta\theta} &=& \frac{2}{r}\boldsymbol{u}_r + \frac{2}{r}\frac{\partial \boldsymbol{u}_\theta}{\partial \theta} \\
\sigma_{\theta\phi} &=& \frac{1}{r} \frac{\partial \boldsymbol{u}_\phi}{\partial \theta} + \frac{1}{r\sin\theta}\frac{\partial \boldsymbol{u}_\theta}{\partial \phi} 
                           - \frac{\cot\theta}{r} \boldsymbol{u}_\phi  \\
\sigma_{\phi\phi} &=&  \frac{2}{r}\boldsymbol{u}_r + \frac{2\cot\theta}{r}\boldsymbol{u}_\theta 
                                + \frac{2}{r\sin\theta} \frac{\partial \boldsymbol{u}_\phi}{\partial \phi}
 @f}

 Once the entries of the matrix are computed the normal is recast into the spherical basis of the source particle at the target point,
 then the matrix-vector product is computed, returning a vector in the same coordintes.
*/
class StokesTractionHarmonic : public SphericalVectorField
{
private:
    SingleLayerHarmonic* slh;
    StokesPressure* sp;

    StokesTractionHarmonic(const StokesTractionHarmonic & t): slh(nullptr) , sp(nullptr) {}
    StokesTractionHarmonic(StokesTractionHarmonic && t) noexcept : slh(nullptr), sp(nullptr) {}
    StokesTractionHarmonic & operator =(const StokesTractionHarmonic & t){}
    StokesTractionHarmonic & operator =(StokesTractionHarmonic && t) noexcept {}

public:

    /**
     * @brief constructs the traction from a SingleLayerHarmonic and StokesPressure opbect.
    */
    StokesTractionHarmonic(SingleLayerHarmonic* u, StokesPressure* p ,std::string na = "Traction Harmonic"): 
        slh(u), 
        sp(p)
    {
        name = na;
    }


    

    /**
     * @brief Evaluates the traction at the target point.
    */
    RectCoord  operator()(const SphereCoord & s) const
    {
        
        SphereCoord temp = recast(s , slh->system);

       
        double sigmaRR = 2.0 * slh->eR_dR(temp);
        double sigmaRTheta = slh->eTheta_dR(temp) + slh->eR_dTheta(temp) / temp.r - slh->eTheta(temp) / temp.r;
        double sigmaRPhi = slh->ePhi_dR(temp) +  slh->eR_dPhi(temp) / (temp.r * sin(temp.s.theta)) - slh->ePhi(temp) / temp.r;

        
        double sigmaThetaTheta =  2.0 / temp.r * (slh->eR(temp) + slh->eTheta_dTheta(temp));
        double sigmaThetaPhi = 1.0 / temp.r * (slh->ePhi_dTheta(temp) + 1.0 / sin(temp.s.theta) * (slh->eTheta_dPhi(temp) - cos(temp.s.theta) * slh->ePhi(temp)));
        
        
        double sigmaPhiPhi =  2.0 / temp.r * (slh->eR(temp) + 1.0 / sin(temp.s.theta) * (cos(temp.s.theta) *slh->eTheta(temp) + slh->ePhi_dPhi(temp)));

        vec<3> n;

        n[0] = dot(e_r(s), e_r(temp));
        n[1] = dot(e_r(s), e_theta(temp));
        n[2] = dot(e_r(s), e_phi(temp));

        Matrix<3, 3> sigma;
        sigma[0][0] = sigmaRR - (*sp)(temp);
        sigma[1][0] = sigma[0][1] = sigmaRTheta;
        sigma[2][0] = sigma[0][2] = sigmaRPhi;
        sigma[1][1] = sigmaThetaTheta - (*sp)(temp);
        sigma[2][1] = sigma[1][2] = sigmaThetaPhi;
        sigma[2][2] = sigmaPhiPhi - (*sp)(temp);

        vec<3> t = dot(sigma, n);

        return t[0] * e_r(temp) + t[1] * e_theta(temp) + t[2] * e_phi(temp);
    }


};

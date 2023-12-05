#pragma once
#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <random>
#include <string>
#include <fstream>

#include "Legendre.h"
#include "SphericalCalc.h"

/**
* \file "SphericalHarmonics.h"
* \brief defninition of the Spherical Harmonic basis functions, and series representation.
* 
* Here we intruduce the Spherical Harmonic basis functions, The eigenfunctions for the Laplacian on a sphere.
* The complex valued scalar spherical harmonics are given in spherical coordinates by 
* @f{equation}{
   Y_n^m(\theta , \phi) = \alpha_{n,m} P_n^m(\cos(\theta)) e^{i m \phi}. \quad |m| \leq n.
 @f}
* Throughout this work, we do not deal in complex numbers, and define the real and imaginary parts of @f$ Y_n^m @f$ separately.
* From the above, we define the Vector Spherical Harmonics to be 
* @f{eqnarray}{
\boldsymbol{V}_n^m &=& -(n+1) Y_n^m \boldsymbol{e}_r + \nabla_\Gamma Y_n^m \\
\boldsymbol{W}_n^m &=& n Y_n^m 
\boldsymbol{e}_r + \nabla_\Gamma Y_n^m \\
\boldsymbol{X}_n^m &=& \boldsymbol{e}_r \times \nabla_\Gamma Y_n^m.
 @f}
 Given a density @f$ \boldsymbol{\rho}: \Gamma \rightarrow R^3 @f$ periodic in the azimuth angle, we can write it in the form  
 @f{equation}{
   \boldsymbol{\rho}(\theta , \phi) =\sum_{n = 0}^\infty\sum_{-n \leq m \leq n} \frac{\hat{\rho}_{n,m}^V}{\|\boldsymbol{V}_n^m\|_2^2}\boldsymbol{V}_n^m(\theta, \phi) +  \frac{\hat{\rho}_{n,m}^W}{\|\boldsymbol{W}_n^m\|_2^2}\boldsymbol{W}_n^m(\theta, \phi)
    +\frac{\hat{\rho}_{n,m}^X}{\|\boldsymbol{X}_n^m\|_2^2}\boldsymbol{X}_n^m(\theta, \phi).
 @f}
 With the @f$ \hat{\rho} @f$'s are found by forming the inner product
  @f{equation}{
    \hat{\rho}_{n,m}^U = \int_{\Gamma} \boldsymbol{\rho}(\theta , \phi) \cdot \bar{\boldsymbol{U}}_n^m(\theta , \phi) dS_x
 @f}
*
*/



/**
* @brief minimum macro.
*/
#define min(a , b) a < b? a : b

/**
* @brief maximum macro.
*/
#define max(a , b) a > b? a : b



/**
 * @brief The Real part of @f$ Y_n^m @f$.
 * 
 * Computes the real part of the scalar spherical harmonic and derivatives up to second order.
*/
class YReal : public SphericalScalarFunction
{
private:
    int m;
    int n;
    Legendre* poly;
    
    double PI = 4.0 * atan(1.0);

public:
    /**
     * @brief Constructor, sets the values of m and n.
     * @param a the order(m).
     * @param b the degree (n).
     * @note we require @f$ 0 \leq m \leq n @f$.
    */
    YReal(Legendre* P,const int a = 0, const int b = 0 ) : poly(P)
    {
        m = a;
        n = b;
        name = "YReal";

        if (m < 0)
            std::cout << "Index error in YReal, m < 0" << std::endl;

        if(m > n)
            std::cout << "Index error in YReal, m > n" << std::endl;
    }

    YReal(const YReal& yr , Legendre* P) : poly(P)
    {
        m = yr.m;
        n = yr.n;

        name = "YReal";
    }

    void reset(int a, int b)
    {
        m = a;
        n = b;
    }

    const YReal & operator =(const YReal& yr)
    {
        m = yr.m;
        n = yr.n;

        poly = yr.poly;

        name = "YReal";

        return *this;
    }

    /**
     * @brief Evaluates @f$ \mathfrak{R}(Y_n^m) @f$ at s.
     * @return @f{equation}{
      \alpha_{n,m} P_n^m(\cos(\theta)) \cos{m \phi}
       @f}
    */
    inline double operator()(const SphereCoord & x) const 
    {



        double coef = sqrt((2.0 * (double)n + 1.0) * (double)factorial(n - m) / (4.0 * PI * (double)factorial(n + m)));
        //std::cout << coef << std::endl;

        coef *= (*poly)(cos(x.s.theta), m, n);

        //std::cout << coef << std::endl;

        return coef * cos((double)m * x.s.phi);


    }

    /**
     * @brief Computes the partial derivative of @f$ Y_n^m @f$ with respect to @f$ \phi @f$.
     * @return  @f{equation}{
      \frac{\partial \mathfrak{R}(Y_n^m)}{\partial \phi} = -m \alpha_{n,m} P_n^m(\cos(\theta)) \sin{m \phi}
       @f}
    */
    inline double dPhi(const SphereCoord & x) const
    {
        double coef = sqrt((2.0 * (double)n + 1.0) * (double)factorial(n - m) / (4.0 * PI * (double)factorial(n + m)));
        //std::cout << coef << std::endl;

        coef *= (*poly)(cos(x.s.theta), m, n);

        //std::cout << coef << std::endl;

        return - coef *m * sin((double)m * x.s.phi);
    }

    /**
     * @brief Computes the second order partial derivative of @f$ Y_n^m @f$ with respect to @f$ \phi @f$ twice.
     * @return  @f{equation}{
      \frac{\partial^2 \mathfrak{R}(Y_n^m)}{\partial \phi^2} = -m^2 \mathfrak{R}(Y_n^m)
       @f}
    */
    inline double dPhidPhi(const SphereCoord & x) const
    {
        double m_(m);
        return -m_ * m_ * operator()(x);
    }

    /**
     * @brief Computes the partial derivative of @f$ Y_n^m @f$ with respect to @f$ \theta @f$.
     * @return @f{equation}{
      \alpha_{n,m} \frac{\partial P_n^m(\cos(\theta))}{\partial \theta} \cos{m \phi}
       @f}
       See also: @ref  Legendre::dTheta()
    */
    inline double dTheta(const SphereCoord & x) const
    {
        double coef = sqrt((2.0 * (double)n + 1.0) * (double)factorial(n - m) / (4.0 * PI * (double)factorial(n + m)));
        //std::cout << coef << std::endl;

        poly->populate(cos(x.s.theta), m + 1, n + 1);

        coef *= poly->dTheta(m, n, x.s.theta);


        return coef * cos((double)m * x.s.phi);
    }

    /**
     * @brief Computes the second partial derivative of @f$ Y_n^m @f$ with respect to @f$ \theta @f$ twice.
     * @return @f{equation}{
      \alpha_{n,m} \frac{\partial^2 P_n^m(\cos(\theta))}{\partial \theta^2} \cos{m \phi}
       @f}
       See also: @ref Legendre::dTheta()
    */
    inline double dThetadTheta(const SphereCoord & x) const
    {
        double coef = sqrt((2.0 * (double)n + 1.0) * (double)factorial(n - m) / (4.0 * PI * (double)factorial(n + m)));
        //std::cout << coef << std::endl;

        poly->populate(cos(x.s.theta), m + 2, n + 2);

        coef *= poly->dTheta(m, n, x.s.theta , 2);


        return coef * cos((double)m * x.s.phi);
    }

    /**
     * @brief Computes the second mixed partial derivative of @f$ Y_n^m @f$ with respect to @f$ \theta @f$ and @f$ \phi @f$.
     * @return @f{equation}{
      -m \alpha_{n,m} \frac{\partial P_n^m(\cos(\theta))}{\partial \theta} \sin{m \phi}
       @f}
       See also: @ref Legendre::dTheta() @ref dPhi()
    */
    inline double dThetadPhi(const SphereCoord & x) const
    {
        double coef = sqrt((2.0 * (double)n + 1.0) * (double)factorial(n - m) / (4.0 * PI * (double)factorial(n + m)));
        //std::cout << coef << std::endl;

        poly->populate(cos(x.s.theta), m + 1, n + 1);

        coef *= poly->dTheta(m, n, x.s.theta);


        return -   coef * m* sin((double)m * x.s.phi);
    }

    

};





/**
 * @brief The Imaginary part of @f$ Y_n^m @f$.
 *
 * Computes the imaginary part of the scalar spherical harmonic and derivatives up to second order.
*/
class YImag : public SphericalScalarFunction
{
private:
    int m;
    int n;
    Legendre* poly;
    double PI = 4.0 * atan(1.0);

public:

    /**
     * @brief Constructor, sets the values of m and n.
     * @param a the order(m).
     * @param b the degree (n).
     * @note we require @f$ 0 \leq m \leq n @f$.
    */
    YImag(Legendre* P ,int a = 0, int b = 0) : poly(P)
    {
        m = a;
        n = b;
        name = "YImag";

        if (m < 0)
            std::cout << "Index error in YImag, m < 0" << std::endl;

        if (m > n)
            std::cout << "Index error in YImag, m > n" << std::endl;

    }

    YImag(const YImag& yi , Legendre * P) : poly(P)
    {
        m = yi.m;
        n = yi.n;

        name = "YReal";
    }

    void reset(int a, int b)
    {
        m = a;
        n = b;
    }

    const YImag& operator =(const YImag& yi)
    {
        m = yi.m;
        n = yi.n;

        poly = yi.poly;

        name = "YReal";

        return *this;
    }


    /**
     * @brief Evaluates @f$ \mathfrak{I}(Y_n^m) @f$ at s.
     * @return @f{equation}{
      \alpha_{n,m} P_n^m(\cos(\theta)) \sin{m \phi}
       @f}
    */
    inline double operator()(const SphereCoord & x) const 
    {


        double coef = sqrt((2.0 * (double)n + 1.0) * (double)factorial(n - m) / (4.0 * PI * (double)factorial(n + m)));
        //std::cout << coef << std::endl;

        coef *= (*poly)(cos(x.s.theta), m, n);

        //std::cout << coef << std::endl;

        return coef * sin((double)m * x.s.phi);


    }

    /**
     * @brief Computes the partial derivative of @f$ Y_n^m @f$ with respect to @f$ \phi @f$.
     * @return  @f{equation}{
      \frac{\partial \mathfrak{R}(Y_n^m)}{\partial \phi} = m \alpha_{n,m} P_n^m(\cos(\theta)) \cos{m \phi}
       @f}
    */
    inline double dPhi(const SphereCoord & x) const
    {
        double coef = sqrt((2.0 * (double)n + 1.0) * (double)factorial(n - m) / (4.0 * PI * (double)factorial(n + m)));
        //std::cout << coef << std::endl;

        coef *= (*poly)(cos(x.s.theta), m, n);

        //std::cout << coef << std::endl;

        return coef * m * cos((double)m * x.s.phi);
    }

    /**
     * @brief Computes the second order partial derivative of @f$ Y_n^m @f$ with respect to @f$ \phi @f$ twice.
     * @return  @f{equation}{
      \frac{\partial^2 \mathfrak{I}(Y_n^m)}{\partial \phi^2} = -m^2 \mathfrak{I}(Y_n^m)
       @f}
    */
    inline double dPhidPhi(const SphereCoord & x) const
    {
        double m_(m);
        return -m_ * m_ * operator()(x);
    }

    /**
     * @brief Computes the partial derivative of @f$ Y_n^m @f$ with respect to @f$ \theta @f$.
     * @return @f{equation}{
      \alpha_{n,m} \frac{\partial P_n^m(\cos(\theta))}{\partial \theta} \cos{m \phi}
       @f}
       See also: @ref  Legendre::dTheta()
    */
    inline double dTheta(const SphereCoord & x) const
    {
        double coef = sqrt((2.0 * (double)n + 1.0) * (double)factorial(n - m) / (4.0 * PI * (double)factorial(n + m)));
        //std::cout << coef << std::endl;

        poly->populate(cos(x.s.theta), m + 1, n + 1);

        coef *= poly->dTheta(m, n, x.s.theta);


        return coef * sin((double)m * x.s.phi);
    }

    /**
     * @brief Computes the second partial derivative of @f$ Y_n^m @f$ with respect to @f$ \theta @f$ twice.
     * @return @f{equation}{
      \alpha_{n,m} \frac{\partial^2 P_n^m(\cos(\theta))}{\partial \theta^2} \cos{m \phi}
       @f}
       See also: @ref Legendre::dTheta()
    */
    inline double dThetadTheta(const SphereCoord & x) const
    {
        double coef = sqrt((2.0 * (double)n + 1.0) * (double)factorial(n - m) / (4.0 * PI * (double)factorial(n + m)));
        //std::cout << coef << std::endl;

        poly->populate(cos(x.s.theta), m + 1, n + 1);

        coef *= poly->dTheta(m, n, x.s.theta , 2);


        return coef * sin((double)m * x.s.phi);
    }

    /**
     * @brief Computes the second mixed partial derivative of @f$ Y_n^m @f$ with respect to @f$ \theta @f$ and @f$ \phi @f$.
     * @return @f{equation}{
      m \alpha_{n,m} \frac{\partial P_n^m(\cos(\theta))}{\partial \theta} \cos{m \phi}
       @f}
       See also: @ref Legendre::dTheta() @ref dPhi()
    */
    inline double dThetadPhi(const SphereCoord & x) const
    {
        double coef = sqrt((2.0 * (double)n + 1.0) * (double)factorial(n - m) / (4.0 * PI * (double)factorial(n + m)));
        //std::cout << coef << std::endl;

        poly->populate(cos(x.s.theta), m + 1, n + 1);

        coef *= poly->dTheta(m, n, x.s.theta);


        return coef * m * cos((double)m * x.s.phi);
    }

};

/**
 * @brief Functor for surface gradient of @f$ \mathfrak{R}(Y_n^m) @f$. Includes derivatives and spherical basis component extractions. 
*/
class GReal : public SphericalVectorField
{
private:
    int m;
    int n;
    YReal yr;

public:
    GReal(Legendre* P): yr(P){}
    /**
     * @brief Assigns the order and degree to the object.
    */
    GReal(Legendre*P ,  int m0, int n0) : yr(P ,m0, n0)
    {
        m = m0;
        n = n0;
    }

    void reset(int a, int b)
    {
        m = a;
        n = b;
        yr.reset(a, b);
    }

    /**
     * @brief The surface gradient of @f$ \mathfrak{R}(Y_n^m) @f$.
     * @return @f{equation}{
     * \nabla_\Gamma \mathfrak{R}(Y_n^m) = \frac{\partial \mathfrak{R}(Y_n^m)}{\partial \theta} \boldsymbol{e}_\theta + \frac{1}{\sin\theta} \frac{\partial \mathfrak{R}(Y_n^m)}{\partial\phi}\boldsymbol{e}_\phi
     * @f} evaluated at x. The partial derivatives are YReal::dPhi() and YReal::dTheta().
    */
    inline RectCoord  operator()(const SphereCoord & x) const
    {
        double dTheta = yr.dTheta(x);

        double dPhi = yr.dPhi(x);

        return dTheta * e_theta(x) + dPhi / sin(x.s.theta) * e_phi(x);


    }

    /**
     * @brief Evaluates the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \nabla_\Gamma \mathfrak{R}(Y_n^m) @f$.
     * @return @f$ \frac{\partial \mathfrak{R}(Y_n^m)}{\partial \theta} @f$ evaluated at x.
     * 
     * See also: YReal::dTheta().
    */
    inline double eTheta(const SphereCoord & x) const
    {
        return yr.dTheta(x);
    }

    /**
     * @brief Evaluates the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \nabla_\Gamma \mathfrak{R}(Y_n^m) @f$.
     * @return @f$ \frac{1}{\sin\theta} \frac{\partial \mathfrak{R}(Y_n^m)}{\partial \phi} @f$ evaluated at x.
     *
     * See also: YReal::dPhi().
    */
    inline double ePhi(const SphereCoord & x) const
    {
        return yr.dPhi(x) / sin(x.s.theta);
    }

    /**
     * @brief Evaluates the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \nabla_\Gamma \mathfrak{R}(Y_n^m)@f$, differentiated with respect to @f$\theta @f$.
     * @return @f$ \frac{\partial^2 \mathfrak{R}(Y_n^m)}{\partial \theta^2} @f$ evaluated at x.
     *
     * See also: YReal::dThetadTheta().
    */
    inline double eTheta_dTheta(const SphereCoord & x) const
    {
        return yr.dThetadTheta(x);
    }

    /**
     * @brief Evaluates the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \nabla_\Gamma \mathfrak{R}(Y_n^m)@f$, differentiated with respect to @f$\phi @f$.
     * @return @f$ \frac{\partial^2 \mathfrak{R}(Y_n^m)}{\partial \theta \partial \phi} @f$ evaluated at x.
     *
     * See also: YReal::dThetadPhi().
    */
    inline double eTheta_dPhi(const SphereCoord & x) const
    {
        return yr.dThetadPhi(x);
    }

    /**
     * @brief Evaluates the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \nabla_\Gamma \mathfrak{R}(Y_n^m)@f$, differentiated with respect to @f$\phi @f$.
     * @return @f$ \frac{1}{\sin \theta}\frac{\partial^2 \mathfrak{R}(Y_n^m)}{\partial \phi^2} @f$ evaluated at x.
     *
     * See also: YReal::dPhidPhi().
    */
    inline double ePhi_dPhi(const SphereCoord & x) const
    {
        return yr.dPhidPhi(x) / sin(x.s.theta);
    }

    /**
     * @brief Evaluates the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \nabla_\Gamma \mathfrak{R}(Y_n^m)@f$, differentiated with respect to @f$\theta @f$.
     * @return @f{equation}{ 
     \frac{\sin \theta \ \frac{\partial^2 \mathfrak{R}(Y_n^m)}{\partial \theta \partial \phi} 
            - \cos \theta \frac{\partial \mathfrak{R}(Y_n^m)}{\partial \phi} }{\sin^2 \theta}
     @f} evaluated at x.
     *
     * See also: YReal::dThetadPhi() , YReal::dPhi().
    */
    inline double ePhi_dTheta(const SphereCoord & x) const
    {
        return (yr.dThetadPhi(x) * sin(x.s.theta) - cos(x.s.theta) * yr.dPhi(x)) / (sin(x.s.theta) * sin(x.s.theta));
    }

    
};



/**
 * @brief Functor for surface gradient of @f$ \mathfrak{I}(Y_n^m) @f$. Includes derivatives and spherical basis component extractions.
*/
class GImag : public SphericalVectorField
{
private:
    int m;
    int n;
    YImag yi;

public:
    GImag(Legendre* P):yi(P) {}


    /**
     * @brief Assigns the order and degree to the object.
    */
    GImag(Legendre* P, int m0, int n0) : yi(P,m0,n0)
    {
        m = m0;
        n = n0;
    }

    void reset(int a, int b)
    {
        m = a;
        n = b;

        yi.reset(a, b);
    }

    
    /**
     * @brief The surface gradient of @f$ \mathfrak{I}(Y_n^m) @f$.
     * @return @f{equation}{
     * \nabla_\Gamma \mathfrak{I}(Y_n^m) = \frac{\partial \mathfrak{I}(Y_n^m)}{\partial \theta} \boldsymbol{e}_\theta + \frac{1}{\sin\theta} \frac{\partial \mathfrak{I}(Y_n^m)}{\partial\phi}\boldsymbol{e}_\phi
     * @f} evaluated at x. The partial derivatives are YImag::dPhi() and YImag::dTheta().
    */
    inline RectCoord  operator()(const SphereCoord & x) const
    {


        double dTheta = yi.dTheta(x);

        double dPhi = yi.dPhi(x);

        return dTheta * e_theta(x) + dPhi / sin(x.s.theta) * e_phi(x);


    }

    /**
     * @brief Evaluates the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \nabla_\Gamma \mathfrak{I}(Y_n^m) @f$.
     * @return @f$ \frac{\partial \mathfrak{I}(Y_n^m)}{\partial \theta} @f$ evaluated at x.
     *
     * See also: YImag::dTheta().
    */
    inline double eTheta(const SphereCoord & x) const
    {
        return yi.dTheta(x);
    }

    /**
     * @brief Evaluates the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \nabla_\Gamma \mathfrak{I}(Y_n^m) @f$.
     * @return @f$ \frac{1}{\sin\theta} \frac{\partial \mathfrak{I}(Y_n^m)}{\partial \phi} @f$ evaluated at x.
     *
     * See also: YImag::dPhi().
    */
    inline double ePhi(const SphereCoord & x) const
    {
        return yi.dPhi(x) / sin(x.s.theta);
    }

    /**
     * @brief Evaluates the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \nabla_\Gamma \mathfrak{I}(Y_n^m)@f$, differentiated with respect to @f$\theta @f$.
     * @return @f$ \frac{\partial^2 \mathfrak{I}(Y_n^m)}{\partial \theta^2} @f$ evaluated at x.
     *
     * See also: YImag::dThetadTheta().
    */
    inline double eTheta_dTheta(const SphereCoord & x) const
    {
        return yi.dThetadTheta(x);
    }

    /**
     * @brief Evaluates the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \nabla_\Gamma \mathfrak{I}(Y_n^m)@f$, differentiated with respect to @f$\phi @f$.
     * @return @f$ \frac{\partial^2 \mathfrak{I}(Y_n^m)}{\partial \theta \partial \phi} @f$ evaluated at x.
     *
     * See also: YImag::dThetadPhi().
    */
    inline double eTheta_dPhi(const SphereCoord & x) const
    {
        return yi.dThetadPhi(x);
    }

    /**
     * @brief Evaluates the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \nabla_\Gamma \mathfrak{I}(Y_n^m)@f$, differentiated with respect to @f$\phi @f$.
     * @return @f$ \frac{1}{\sin \theta}\frac{\partial^2 \mathfrak{I}(Y_n^m)}{\partial \phi^2} @f$ evaluated at x.
     *
     * See also: YImag::dPhidPhi().
    */
    inline double ePhi_dPhi(const SphereCoord & x) const
    {
        return yi.dPhidPhi(x) / sin(x.s.theta);
    }

    /**
     * @brief Evaluates the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \nabla_\Gamma \mathfrak{I}(Y_n^m)@f$, differentiated with respect to @f$\theta @f$.
     * @return @f{equation}{
     \frac{\sin \theta \ \frac{\partial^2 \mathfrak{I}(Y_n^m)}{\partial \theta \partial \phi}
            - \cos \theta \frac{\partial \mathfrak{I}(Y_n^m)}{\partial \phi} }{\sin^2 \theta}
     @f} evaluated at x.
     *
     * See also: YImag::dThetadPhi() , YImag::dPhi().
    */
    inline double ePhi_dTheta(const SphereCoord & x) const
    {
        return (yi.dThetadPhi(x) * sin(x.s.theta) - cos(x.s.theta) * yi.dPhi(x)) / (sin(x.s.theta) * sin(x.s.theta));
    }

    

};

/**
 * @brief Functor for @f$ \mathfrak{R}(\boldsymbol{V}_n^m) @f$. Also has functions for component extraction and differentiation.
*/
class VReal : public SphericalVectorField
{
    int m;
    int n;
    YReal yr;
    GReal gr;

public:


    /**
     * @brief Initializes the order and degree.
    */
    VReal(Legendre*P ,  int a = 0, int b = 0) : yr(P , a , b) , gr(P ,a,b)
    {
        m = a;
        n = b;
    }

    /**
     * @brief Evaluates  @f$ \mathfrak{R}(\boldsymbol{V}_n^m) @f$.
     * @return @f$ -(n+1) \mathfrak{R}(\boldsymbol{Y}_n^m)(x) \boldsymbol{e}_rs + \nabla_\Gamma \mathfrak{R}(\boldsymbol{Y}_n^m)(x) @f$.
     * 
     * See YReal::operator()() , GReal::operator()().
    */
    inline RectCoord operator()(const SphereCoord & x) const 
    {
        return -((double)n + 1.0) * yr(x) * e_r(x) + gr(x);
    }

    /**
     * @brief resets the indices.
    */
    void reset(int a, int b)
    {
        m = a;
        n = b;
        yr.reset(a, b);
        gr.reset(a, b);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of @f$ \mathfrak{R}(\boldsymbol{V}_n^m) @f$.
     * @return  @f$ -(n+1) \mathfrak{R}(\boldsymbol{Y}_n^m)(x) @f$.
     * 
     * See: YReal.
    */
    inline double eR(const SphereCoord & x) const
    {
        return -((double)n + 1.0) * yr(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \mathfrak{R}(\boldsymbol{V}_n^m) @f$.
     * @return  @f$ GReal::eTheta() @f$.
    */
    inline double eTheta(const SphereCoord & x) const
    {
        return gr.eTheta(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \mathfrak{R}(\boldsymbol{V}_n^m) @f$.
     * @return GReal::ePhi().
    */
    inline double ePhi(const SphereCoord & x) const
    {
        return gr.ePhi(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \mathfrak{R}(\boldsymbol{V}_n^m) @f$,then differentiate with respect to @f$\theta@f$.
     * @return GReal::eTheta_dTheta().
    */
    inline double eTheta_dTheta(const SphereCoord & x) const
    {
        return gr.eTheta_dTheta(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \mathfrak{R}(\boldsymbol{V}_n^m) @f$,then differentiate with respect to @f$\phi@f$.
     * @return GReal::eTheta_dPhi().
    */
    inline double eTheta_dPhi(const SphereCoord & x) const
    {
        return gr.eTheta_dPhi(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \mathfrak{R}(\boldsymbol{V}_n^m) @f$,then differentiate with respect to @f$\phi@f$.
     * @return GReal::ePhi_dPhi().
    */
    inline double ePhi_dPhi(const SphereCoord & x) const
    {
        return gr.ePhi_dPhi(x);
    }
    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \mathfrak{R}(\boldsymbol{V}_n^m) @f$,then differentiate with respect to @f$\theta@f$.
     * @return GReal::ePhi_dTheta().
    */
    inline double ePhi_dTheta(const SphereCoord & x) const
    {
        return gr.ePhi_dTheta(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of @f$ \mathfrak{R}(\boldsymbol{V}_n^m) @f$,then differentiate with respect to @f$\theta@f$.
     * @return  @f$ -(n+1)\frac{\partial \mathfrak{R}(\boldsymbol{Y}_n^m)}{\partial \theta}  @f$.
     * 
     * See: YReal::dTheta()
    */
    inline double eR_dTheta(const SphereCoord & x) const
    {
        return -((double) n+1.0)*yr.dTheta(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of @f$ \mathfrak{R}(\boldsymbol{V}_n^m) @f$,then differentiate with respect to @f$\phi@f$.
     * @return  @f$ -(n+1)\frac{\partial \mathfrak{R}(\boldsymbol{Y}_n^m)}{\partial \phi}  @f$.
     *
     * See: YReal::dPhi()
    */
    inline double eR_dPhi(const SphereCoord & x) const
    {
        return -((double)n + 1.0) * yr.dPhi(x);
    }

    

};

/**
 * @brief Functor for @f$ \mathfrak{I}(\boldsymbol{V}_n^m) @f$. Also has functions for component extraction and differentiation.
*/
class VImag : public SphericalVectorField
{
    int m;
    int n;
    YImag yi;
    GImag gi;

public:

    /**
     * @brief Initializes the order and degree.
    */
    VImag(Legendre*P , int a = 0, int b = 0) : yi(P , a, b), gi(P , a, b)
    {
        m = a;
        n = b;
    }

    /**
     * @brief Evaluates  @f$ \mathfrak{I}(\boldsymbol{V}_n^m) @f$.
     * @return @f$ -(n+1) \mathfrak{I}(\boldsymbol{Y}_n^m)(x) \boldsymbol{e}_rs + \nabla_\Gamma \mathfrak{I}(\boldsymbol{Y}_n^m)(x) @f$.
     *
     * See YImag::operator()() , GImag::operator()().
    */
    inline RectCoord operator()(const SphereCoord & x) const 
    {
        return -((double)n + 1.0) * yi(x) * e_r(x) + gi(x);
    }


    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of @f$ \mathfrak{I}(\boldsymbol{V}_n^m) @f$.
     * @return  @f$ -(n+1) \mathfrak{I}(\boldsymbol{Y}_n^m)(x) @f$.
     *
     * See: YImag::operator()().
    */
    inline double eR(const SphereCoord & x) const
    {
        return -((double)n + 1.0) * yi(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \mathfrak{I}(\boldsymbol{V}_n^m) @f$.
     * @return GImag::eTheta().
    */
    inline double eTheta(const SphereCoord & x) const
    {
        return gi.eTheta(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \mathfrak{I}(\boldsymbol{V}_n^m) @f$.
     * @return GImag::ePhi().
    */
    inline double ePhi(const SphereCoord & x) const
    {
        return gi.ePhi(x);
    }

    /**
    * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \mathfrak{I}(\boldsymbol{V}_n^m) @f$,then differentiate with respect to @f$\theta@f$.
    * @return GImag::eTheta_dTheta().
    */
    inline double eTheta_dTheta(const SphereCoord & x) const
    {
        return gi.eTheta_dTheta(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \mathfrak{I}(\boldsymbol{V}_n^m) @f$,then differentiate with respect to @f$\phi@f$.
     * @return GImag::eTheta_dPhi().
    */
    inline double eTheta_dPhi(const SphereCoord & x) const
    {
        return gi.eTheta_dPhi(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \mathfrak{I}(\boldsymbol{V}_n^m) @f$,then differentiate with respect to @f$\phi@f$.
     * @return GImag::ePhi_dPhi().
    */
    inline double ePhi_dPhi(const SphereCoord & x) const
    {
        return gi.ePhi_dPhi(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \mathfrak{I}(\boldsymbol{V}_n^m) @f$,then differentiate with respect to @f$\theta@f$.
     * @return GImag::ePhi_dTheta().
    */
    inline double ePhi_dTheta(const SphereCoord& x) const
    {
        return gi.ePhi_dTheta(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of @f$ \mathfrak{I}(\boldsymbol{V}_n^m) @f$,then differentiate with respect to @f$\theta@f$.
     * @return  @f$ -(n+1)\frac{\partial \mathfrak{R}(\boldsymbol{Y}_n^m)}{\partial \theta}  @f$.
     *
     * See: YImag::dTheta()
    */
    inline double eR_dTheta(const SphereCoord& x) const
    {
        return -((double)n + 1.0) * yi.dTheta(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of @f$ \mathfrak{I}(\boldsymbol{V}_n^m) @f$,then differentiate with respect to @f$\phi@f$.
     * @return  @f$ -(n+1)\frac{\partial \mathfrak{I}(\boldsymbol{Y}_n^m)}{\partial \phi}  @f$.
     *
     * See: YImag::dPhi()
    */
    inline double eR_dPhi(const SphereCoord & x) const
    {
        return -(double)(n + 1) * yi.dPhi(x);
    }

    /**
     * @brief resets the indices.
    */
    void reset(int a, int b)
    {
        m = a;
        n = b;
        yi.reset(a, b);
        gi.reset(a, b);
    }

    

};


/**
 * @brief Functor for @f$ \mathfrak{R}(\boldsymbol{W}_n^m) @f$. Also has functions for component extraction and differentiation.
*/
class WReal : public SphericalVectorField
{
    int m;
    int n;
    YReal yr;
    GReal gr;

public:

    /**
     * @brief Initializes the order and degree.
    */
    WReal(Legendre*P, int a = 0, int b = 0) : yr(P, a, b), gr(P, a, b)
    {
        m = a;
        n = b;
    }

    /**
     * @brief Evaluates  @f$ \mathfrak{R}(\boldsymbol{W}_n^m) @f$.
     * @return @f$ n \mathfrak{R}(\boldsymbol{Y}_n^m)(x) \boldsymbol{e}_rs + \nabla_\Gamma \mathfrak{R}(\boldsymbol{Y}_n^m)(x) @f$.
     *
     * See YReal::operator()() , GReal::operator()().
    */
    inline RectCoord  operator()(const SphereCoord & x) const
    {
        return (double)n * yr(x) * e_r(x) + gr(x);
    }

    /**
     * @brief resets the indices.
    */
    void reset(int a, int b)
    {
        m = a;
        n = b;
        yr.reset(m, n);
        gr.reset(m, n);
    }


    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of @f$ \mathfrak{R}(\boldsymbol{W}_n^m) @f$.
     * @return  @f$ n \mathfrak{R}(\boldsymbol{Y}_n^m)(x) @f$.
     *
     * See: YReal::operator()().
    */
    inline double eR(const SphereCoord & x) const
    {
        return n * yr(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \mathfrak{R}(\boldsymbol{W}_n^m) @f$.
     * @return GReal::eTheta().
    */
    inline double eTheta(const SphereCoord & x) const
    {
        return gr.eTheta(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \mathfrak{R}(\boldsymbol{W}_n^m) @f$.
     * @return GReal::ePhi().
    */
    inline double ePhi(const SphereCoord & x) const
    {
        return gr.ePhi(x);
    }

    /**
    * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \mathfrak{R}(\boldsymbol{W}_n^m) @f$,then differentiate with respect to @f$\theta@f$.
    * @return GReal::eTheta_dTheta().
    */
    inline double eTheta_dTheta(const SphereCoord & x) const
    {
        return gr.eTheta_dTheta(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \mathfrak{R}(\boldsymbol{W}_n^m) @f$,then differentiate with respect to @f$\phi@f$.
     * @return GReal::eTheta_dPhi().
    */
    inline double eTheta_dPhi(const SphereCoord & x) const
    {
        return gr.eTheta_dPhi(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \mathfrak{R}(\boldsymbol{W}_n^m) @f$,then differentiate with respect to @f$\phi@f$.
     * @return GReal::ePhi_dPhi().
    */
    inline double ePhi_dPhi(const SphereCoord & x) const
    {
        return gr.ePhi_dPhi(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \mathfrak{R}(\boldsymbol{W}_n^m) @f$,then differentiate with respect to @f$\theta@f$.
     * @return GReal::ePhi_dTheta().
    */
    inline double ePhi_dTheta(const SphereCoord& x) const
    {
        return gr.ePhi_dTheta(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of @f$ \mathfrak{R}(\boldsymbol{W}_n^m) @f$,then differentiate with respect to @f$\theta@f$.
     * @return  @f$ n \frac{\partial \mathfrak{R}(\boldsymbol{Y}_n^m)}{\partial \theta}  @f$.
     *
     * See: YReal::dTheta()
    */
    inline double eR_dTheta(const SphereCoord& x) const
    {
        return (double)n * yr.dTheta(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of @f$ \mathfrak{R}(\boldsymbol{W}_n^m) @f$,then differentiate with respect to @f$\phi@f$.
     * @return  @f$ n \frac{\partial \mathfrak{R}(\boldsymbol{Y}_n^m)}{\partial \phi}  @f$.
     *
     * See: YReal::dPhi()
    */
    inline double eR_dPhi(const SphereCoord & x) const
    {
        return (double)n * yr.dPhi(x);
    }

    

};

/**
 * @brief Functor for @f$ \mathfrak{I}(\boldsymbol{W}_n^m) @f$. Also has functions for component extraction and differentiation.
*/
class WImag : public SphericalVectorField
{
    int m;
    int n;
    YImag yi;
    GImag gi;

public:


    /**
     * @brief Initializes the order and degree.
    */
    WImag(Legendre*P, int a = 0, int b = 0) : yi(P, a, b), gi(P, a, b)
    {
        m = a;
        n = b;
    }

    /**
     * @brief Evaluates  @f$ \mathfrak{I}(\boldsymbol{W}_n^m) @f$.
     * @return @f$ n \mathfrak{I}(\boldsymbol{Y}_n^m)(x) \boldsymbol{e}_rs + \nabla_\Gamma \mathfrak{I}(\boldsymbol{Y}_n^m)(x) @f$.
     *
     * See YImag::operator()() , GImag::operator()().
    */
    inline RectCoord operator()(const SphereCoord & x) const
    {
        return (double)n * yi(x) * e_r(x) + gi(x);
    }

    /**
     * @brief resets the indices.
    */
    void reset(int a, int b)
    {
        m = a;
        n = b;
        yi.reset(m, n);
        gi.reset(m, n);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of @f$ \mathfrak{I}(\boldsymbol{W}_n^m) @f$.
     * @return  @f$ n \mathfrak{I}(\boldsymbol{Y}_n^m)(x) @f$.
     *
     * See: YImag::operator()().
    */
    inline double eR(const SphereCoord & x) const
    {
        return n * yi(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \mathfrak{I}(\boldsymbol{W}_n^m) @f$.
     * @return GImag::eTheta().
    */
    inline double eTheta(const SphereCoord & x) const
    {
        return gi.eTheta(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \mathfrak{I}(\boldsymbol{W}_n^m) @f$.
     * @return GImag::ePhi().
    */
    inline double ePhi(const SphereCoord & x) const
    {
        return gi.ePhi(x);
    }

    /**
    * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \mathfrak{I}(\boldsymbol{W}_n^m) @f$,then differentiate with respect to @f$\theta@f$.
    * @return GImag::eTheta_dTheta().
    */
    inline double eTheta_dTheta(const SphereCoord & x) const
    {
        return gi.eTheta_dTheta(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \mathfrak{I}(\boldsymbol{W}_n^m) @f$,then differentiate with respect to @f$\phi@f$.
     * @return GImag::eTheta_dPhi().
    */
    inline double eTheta_dPhi(const SphereCoord & x) const
    {
        return gi.eTheta_dPhi(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \mathfrak{I}(\boldsymbol{W}_n^m) @f$,then differentiate with respect to @f$\phi@f$.
     * @return GImag::ePhi_dPhi().
    */
    inline double ePhi_dPhi(const SphereCoord & x) const
    {
        return gi.ePhi_dPhi(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \mathfrak{I}(\boldsymbol{W}_n^m) @f$,then differentiate with respect to @f$\theta@f$.
     * @return GImag::ePhi_dTheta().
    */
    inline double ePhi_dTheta(const SphereCoord& x) const
    {
        return gi.ePhi_dTheta(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of @f$ \mathfrak{I}(\boldsymbol{W}_n^m) @f$,then differentiate with respect to @f$\theta@f$.
     * @return  @f$ n \frac{\partial \mathfrak{R}(\boldsymbol{Y}_n^m)}{\partial \theta}  @f$.
     *
     * See: YImag::dTheta()
    */
    inline double eR_dTheta(const SphereCoord& x) const
    {
        return (double)n * yi.dTheta(x);
    }


    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of @f$ \mathfrak{I}(\boldsymbol{W}_n^m) @f$,then differentiate with respect to @f$\phi@f$.
     * @return  @f$ n \frac{\partial \mathfrak{I}(\boldsymbol{Y}_n^m)}{\partial \phi}  @f$.
     *
     * See: YImag::dPhi()
    */
    inline double eR_dPhi(const SphereCoord& x) const
    {
        return (double)n * yi.dPhi(x);
    }

    
};



/**
 * @brief Functor for @f$ \mathfrak{R}(\boldsymbol{X}_n^m) @f$. Also has functions for component extraction and differentiation.
*/
class XReal : public SphericalVectorField
{
    int m;
    int n;
    GReal gr;

public:

    /**
     * @brief Initializes the order and degree.
    */
    XReal(Legendre* P, int a = 0, int b = 0) : gr(P, a, b)
    {
        m = a;
        n = b;
    }

    /**
     * @brief Evaluates  @f$ \mathfrak{R}(\boldsymbol{X}_n^m) @f$.
     * @return @f$  \boldsymbol{e}_rs \times \nabla_\Gamma \mathfrak{R}(\boldsymbol{Y}_n^m)(x) @f$.
     *
     * See GReal::operator()().
    */
    inline RectCoord operator()(const SphereCoord & x) const
    {
        return cross(e_r(x), gr(x));
    }

    /**
     * @brief resets the indices.
    */
    void reset(int a, int b)
    {
        m = a;
        n = b;
        gr.reset(m, n);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of @f$ \mathfrak{R}(\boldsymbol{X}_n^m) @f$.
     * @return 0.
    */
    inline double eR(const SphereCoord & x) const
    {
        return 0;
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \mathfrak{R}(\boldsymbol{X}_n^m) @f$.
     * @return the negative of GReal::ePhi().
    */
    inline double eTheta(const SphereCoord & x) const
    {
        return -gr.ePhi(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \mathfrak{R}(\boldsymbol{X}_n^m) @f$.
     * @return GReal::eTheta().
    */
    inline double ePhi(const SphereCoord & x) const
    {
        return gr.eTheta(x);
    }

    /**
    * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \mathfrak{R}(\boldsymbol{X}_n^m) @f$,then differentiate with respect to @f$\theta@f$.
    * @return The negative of GReal::ePhi_dTheta().
    */
    inline double eTheta_dTheta(const SphereCoord & x) const
    {
        return -gr.ePhi_dTheta(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \mathfrak{R}(\boldsymbol{X}_n^m) @f$,then differentiate with respect to @f$\phi@f$.
     * @return The negative of GReal::ePhi_dPhi().
    */
    inline double eTheta_dPhi(const SphereCoord & x) const
    {
        return -gr.ePhi_dPhi(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \mathfrak{R}(\boldsymbol{X}_n^m) @f$,then differentiate with respect to @f$\phi@f$.
     * @return GReal::eTheta_dPhi().
    */
    inline double ePhi_dPhi(const SphereCoord & x) const
    {
        return gr.eTheta_dPhi(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \mathfrak{R}(\boldsymbol{X}_n^m) @f$,then differentiate with respect to @f$\theta@f$.
     * @return GReal::eTheta_dTheta().
    */
    inline double ePhi_dTheta(const SphereCoord& x) const
    {
        return gr.eTheta_dTheta(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of @f$ \mathfrak{R}(\boldsymbol{X}_n^m) @f$,then differentiate with respect to @f$\theta@f$.
     * @return 0.
    */
    inline double eR_dTheta(const SphereCoord& x) const
    {
        return 0;
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of @f$ \mathfrak{R}(\boldsymbol{X}_n^m) @f$,then differentiate with respect to @f$\phi@f$.
     * @return 0.
    */
    inline double eR_dPhi(SphereCoord x) const
    {
        return 0;
    }

    
};

/**
 * @brief Functor for @f$ \mathfrak{I}(\boldsymbol{X}_n^m) @f$. Also has functions for component extraction and differentiation.
*/
class XImag : public SphericalVectorField
{
    int m;
    int n;
    GImag gi;

public:

    /**
     * @brief Initializes the order and degree.
    */
    XImag(Legendre*P,int a = 0, int b = 0) : gi(P, a, b)
    {
        m = a;
        n = b;
    }

    /**
     * @brief Evaluates  @f$ \mathfrak{I}(\boldsymbol{X}_n^m) @f$.
     * @return @f$  \boldsymbol{e}_rs \times \nabla_\Gamma \mathfrak{I}(\boldsymbol{Y}_n^m)(x) @f$.
     *
     * See GImag::operator()().
    */
    inline RectCoord operator()(const SphereCoord & x) const
    {
        return cross(e_r(x), gi(x));
    }

    /**
     * @brief resets the indices.
    */
    void reset(int a, int b)
    {
        m = a;
        n = b;
        gi.reset(m, n);
    }
    
    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of @f$ \mathfrak{I}(\boldsymbol{X}_n^m) @f$.
     * @return 0.
    */
    inline double eR(const SphereCoord & x) const
    {
        return 0;
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \mathfrak{I}(\boldsymbol{X}_n^m) @f$.
     * @return The negative GImag::ePhi().
    */
    inline double eTheta(const SphereCoord & x) const
    {
        return -gi.ePhi(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \mathfrak{I}(\boldsymbol{X}_n^m) @f$.
     * @return GImag::eTheta().
    */
    inline double ePhi(const SphereCoord & x) const
    {
        return gi.eTheta(x);
    }

    /**
    * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \mathfrak{I}(\boldsymbol{X}_n^m) @f$,then differentiate with respect to @f$\theta@f$.
    * @return The negative of GImag::ePhi_dTheta().
    */
    inline double eTheta_dTheta(const SphereCoord & x) const
    {
        return -gi.ePhi_dTheta(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\theta @f$ component of @f$ \mathfrak{I}(\boldsymbol{X}_n^m) @f$,then differentiate with respect to @f$\phi@f$.
     * @return The negative of GImag::ePhi_dPhi().
    */
    inline double eTheta_dPhi(const SphereCoord & x) const
    {
        return -gi.ePhi_dPhi(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \mathfrak{I}(\boldsymbol{X}_n^m) @f$,then differentiate with respect to @f$\phi@f$.
     * @return GImag::eTheta_dPhi().
    */
    inline double ePhi_dPhi(const SphereCoord & x) const
    {
        return gi.eTheta_dPhi(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_\phi @f$ component of @f$ \mathfrak{I}(\boldsymbol{X}_n^m) @f$,then differentiate with respect to @f$\theta@f$.
     * @return GImag::eTheta_dTheta().
    */
    inline double ePhi_dTheta(const SphereCoord& x) const
    {
        return gi.eTheta_dTheta(x);
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of @f$ \mathfrak{I}(\boldsymbol{X}_n^m) @f$,then differentiate with respect to @f$\theta@f$.
     * @return 0.
    */
    inline double eR_dTheta(const SphereCoord& x) const
    {
        return 0;
    }

    /**
     * @brief Extracts the @f$ \boldsymbol{e}_r @f$ component of @f$ \mathfrak{I}(\boldsymbol{X}_n^m) @f$,then differentiate with respect to @f$\phi@f$.
     * @return 0.
    */
    inline double eR_dPhi(SphereCoord x) const
    {
        return 0;
    }

    

};

/**
 * @brief Single term in the VSHSeries class. Contains each Vector Spherical Harmonic Basis function for as given index. 
*/
class VSHTerm : public SphericalVectorField
{
private:

    int m;
    int n;
    VReal vr;
    VImag vi;
    WReal wr;
    WImag wi;
    XReal xr;
    XImag xi;
public:
    double rhohatVr;
    double rhohatVi;
    double rhohatWr;
    double rhohatWi;
    double rhohatXr;
    double rhohatXi;



    VSHTerm(Legendre*P = nullptr):
    vr(P),
    vi(P),
    wr(P),
    wi(P),
    xr(P),
    xi(P){}

    /**
     * @brief Constructs a term in the series given the normalized coefficients.
     * @param rhohats in order, the VReal, VImag, WReal, WImag, XReal, XImag coefficients to be assigned.
    */
    VSHTerm(Legendre* P, int m0, int n0, std::array<double, 6> rhohats): 
        vr(P ,m0,n0), 
        vi(P,m0,n0),
        wr(P,m0,n0),
        wi(P,m0,n0),
        xr(P,m0,n0),
        xi(P,m0,n0)
    {
        m = m0;
        n = n0;

        rhohatVr = rhohats[0];
        rhohatVi = rhohats[1];
        rhohatWr = rhohats[2];
        rhohatWi = rhohats[3];
        rhohatXr = rhohats[4];
        rhohatXi = rhohats[5];

    }

    VSHTerm(const VSHTerm& term , Legendre*P):
        vr(P, term.m, term.n),
        vi(P, term.m, term.n),
        wr(P, term.m, term.n),
        wi(P, term.m, term.n),
        xr(P, term.m, term.n),
        xi(P, term.m, term.n)
    {
        m = term.m;
        n = term.n;



        rhohatVr = term.rhohatVr;
        rhohatVi = term.rhohatVi;
        rhohatWr = term.rhohatWr;
        rhohatWi = term.rhohatWi;
        rhohatXr = term.rhohatXr;
        rhohatXi = term.rhohatXi;
    }

    VSHTerm operator =(const VSHTerm& term)
    {
        m = term.m;
        n = term.n;


        vr.reset(m, n);
        vi.reset(m, n);
        wr.reset(m, n);
        wi.reset(m, n);
        xr.reset(m, n);
        xi.reset(m, n);

        rhohatVr = term.rhohatVr;
        rhohatVi = term.rhohatVi;
        rhohatWr = term.rhohatWr;
        rhohatWi = term.rhohatWi;
        rhohatXr = term.rhohatXr;
        rhohatXi = term.rhohatXi;

        return *this;
    }

    /**
     * @brief Addition operator between terms.
     * @param t The term to add.
     * @return A term whose coefficients are sums of the two terms.
    */
    VSHTerm plus(const VSHTerm & t , Legendre*P) const
    {
        std::array<double, 6> temp;

        temp[0] = rhohatVr + t.rhohatVr;
        temp[1] = rhohatVi + t.rhohatVi;
        temp[2] = rhohatWr + t.rhohatWr;
        temp[3] = rhohatWi + t.rhohatWi;
        temp[4] = rhohatXr + t.rhohatXr;
        temp[5] = rhohatXi + t.rhohatXi;

        return VSHTerm(P,m, n, temp);
    }

    /**
     * @brief Subtraction operator between terms.
     * @param t The term to add.
     * @return A term whose coefficients are the difference of the two terms.
    */
    VSHTerm minus(const VSHTerm & t , Legendre*P) const
    {
        std::array<double, 6> temp;

        temp[0] = rhohatVr - t.rhohatVr;
        temp[1] = rhohatVi - t.rhohatVi;
        temp[2] = rhohatWr - t.rhohatWr;
        temp[3] = rhohatWi - t.rhohatWi;
        temp[4] = rhohatXr - t.rhohatXr;
        temp[5] = rhohatXi - t.rhohatXi;

        return VSHTerm(P,m, n, temp);
    }

    /**
     * @brief Unry minus operator.
     * @return A term whose coefficients are the negative of the calling term.
    */
    VSHTerm uminus(Legendre*P) const
    {
        std::array<double, 6> temp;

        temp[0] = -rhohatVr;
        temp[1] = -rhohatVi;
        temp[2] = -rhohatWr;
        temp[3] = -rhohatWi;
        temp[4] = -rhohatXr;
        temp[5] = -rhohatXi;

        return VSHTerm(P,m, n, temp);
    }


    VSHTerm times(const double& a , Legendre*P) const
    {
        std::array<double, 6> temp;

        temp[0] = rhohatVr *a;
        temp[1] = rhohatVi * a;
        temp[2] = rhohatWr * a;
        temp[3] = rhohatWi * a;
        temp[4] = rhohatXr * a;
        temp[5] = rhohatXi * a;

        return VSHTerm(P,m, n, temp);
    }


    /**
     * @brief Computes the inner product between terms.
    */
    double dot(const VSHTerm & s) const
    {
        if (m != s.m || n != s.n)
            return 0.0;
        else
            return rhohatVr * s.rhohatVr + rhohatVi * s.rhohatVi + rhohatWr * s.rhohatWr + rhohatWi * s.rhohatWi + rhohatXr * s.rhohatXr + rhohatXi * s.rhohatXi;
    }

    void axpy(const VSHTerm& term, const double& alpha)
    {

        rhohatVr +=  term.rhohatVr * alpha;
        rhohatVi += term.rhohatVi * alpha;
        rhohatWr += term.rhohatWr * alpha;
        rhohatWi += term.rhohatWi * alpha;
        rhohatXr += term.rhohatXr * alpha;
        rhohatXi += term.rhohatXi * alpha;

    }

    VSHTerm& operator +=(const VSHTerm& term)
    {
        rhohatVr += term.rhohatVr;
        rhohatVi += term.rhohatVi;
        rhohatWr += term.rhohatWr;
        rhohatWi += term.rhohatWi;
        rhohatXr += term.rhohatXr;
        rhohatXi += term.rhohatXi;

        return *this;
    }


    VSHTerm& operator *=(const double & a)
    {
        rhohatVr *= a;
        rhohatVi *= a;
        rhohatWr *= a;
        rhohatWi *= a;
        rhohatXr *= a;
        rhohatXi *= a;

        return *this;
    }

    /**
     * @brief Evaluates the term at the point x.
     * @return In terms of complex numbers, @f$   \hat{\rho}_{n,m}^{V*}\boldsymbol{V}_n^m(\theta, \phi) +  \hat{\rho}_{n,m}^{W*}\boldsymbol{W}_n^m(\theta, \phi)
    +\hat{\rho}_{n,m}^{X*}\boldsymbol{X}_n^m(\theta, \phi)  @f$. The evaluation uses real functions.
    */
    RectCoord operator()(const SphereCoord & x) const
    {
        if (n == 0)
            return 2 * rhohatVr * vr(x);

        RectCoord temp = rhohatVr * vr(x) + rhohatVi * vi(x);
        temp = temp + (rhohatWr * wr(x) + rhohatWi * wi(x));
        return temp + (rhohatXr * xr(x) + rhohatXi * xi(x));
    }


};

/**
 * @brief The Vector Spherical Harmonic Series Functor. Computes and evaluates the fourier series representation of a SphericalVectorField.
*/
class VSHSeries : public VectorFieldSum
{

private:



    
    QuadConstants* consts;
    

public:
    int N;
    SphereCoordSystem* system;
    Legendre* P;
    std::vector<std::vector<VSHTerm>> terms;



    /**
     * @brief Constructs the series. n is the maximum value of the degree. there will be @f$ n(n+1)/2 @f$ terms total.
     * @param c The center of the coordinate system of the series.
    */
    VSHSeries(int n ,QuadConstants* cts,Legendre*p, SphereCoordSystem* sys)
    {
        N = n;

        terms.resize(n + 1);

        for (int i = 0; i <= n; i++)
            terms[i].resize(i + 1);

        system = sys;
        
        consts = cts;

        P = p;


    }

    VSHSeries(const VSHSeries& vsh)
    {
        N = vsh.N;

        terms = vsh.terms;

        system = vsh.system;

        consts = vsh.consts;

        P = vsh.P;

    }

    VSHSeries(const VSHSeries & vsh , Legendre* p) 
    {
        N = vsh.N;

        terms = vsh.terms;

        system = vsh.system;

        consts = vsh.consts; 

        P = p;

    }

    VSHSeries(const VSHSeries& vsh , QuadConstants * cts , Legendre*p , SphereCoordSystem* sys)
    {
        N = vsh.N;

        terms = vsh.terms;

        system = sys;

        consts = cts;

        P = p;

    }


    VSHSeries & operator =(const VSHSeries& vsh)
    {
        N = vsh.N;

        terms = vsh.terms;

        system = vsh.system;

        return *this;

    }


    /**
     * @brief Approximates the passed SphericalVectorField.
     * @param r the function to approximate.
     * @param n number of terms in the series.
    */
    void approximate(SphericalVectorField* r, int p = 0)
    {

        if (p > 0)
        {
            terms.resize(p + 1);

            for (int i = 0; i <= p; i++)
                terms[i].resize( i + 1);

        }




        std::array<double, 6> rhohats;
        VReal vr(P);
        VImag vi(P);
        WReal wr(P);
        WImag wi(P);
        XReal xr(P);
        XImag xi(P);

        for (int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
            {
                if (m == 0)
                {
                    double n_ = double(n);
                    vr.reset(m, n);
                    rhohats[0] =  L2InnerProduct(r, &vr, system , *consts) / (double)((2 * n_ + 1) * (n_ + 1));

                    vi.reset(m, n);
                    rhohats[1] = 0.0;

                    wr.reset(m, n);
                    rhohats[2] =  L2InnerProduct(r, &wr , system , *consts) / (double)((2 * n_ + 1) * n_);

                    wi.reset(m, n);
                    rhohats[3] = 0.0;

                    xr.reset(m, n);
                    rhohats[4] =  L2InnerProduct(r, &xr, system , *consts) / (double)((n_ + 1) * n_);

                    xi.reset(m, n);
                    rhohats[5] = 0.0;

                    terms[n][m] = VSHTerm(P,m, n, rhohats);
                }
                else if (m > 0)
                {
                    double n_ = n;
                    vr.reset(m, n);
                    rhohats[0] = 2.0 * L2InnerProduct(r, &vr, system , *consts) / (double)((2 * n_ + 1) * (n_ + 1));

                    vi.reset(m, n);
                    rhohats[1] = 2.0 * L2InnerProduct(r, &vi, system , *consts) / (double)((2 * n_ + 1) * (n_ + 1));

                    wr.reset(m, n);
                    rhohats[2] = 2.0 * L2InnerProduct(r, &wr, system , *consts) / (double)((2 * n_ + 1) * n_);

                    wi.reset(m, n);
                    rhohats[3] = 2.0 * L2InnerProduct(r, &wi, system, *consts) / (double)((2 * n_ + 1) * n_);

                    xr.reset(m, n);
                    rhohats[4] = 2.0 * L2InnerProduct(r, &xr, system , *consts) / (double)((n_ + 1) * n_);

                    xi.reset(m, n);
                    rhohats[5] = 2.0 * L2InnerProduct(r, &xi, system , *consts) / (double)((n_ + 1) * n_);

                    terms[n][m] = VSHTerm(P,m, n, rhohats);
                }
                
                append(&terms[n][m]);
            }
    }

    /**
     * @brief Approximates A function defined on the surface using data in quadrature nodes.
     * @param data boundary data.
     * @param n maximum degree of the series.
    */
    void approximate(SphereData & data, int p = 0)
    {

        if (p > 0)
        {
            N = p;
            terms.resize(p + 1);

            for (int i = 0; i <= p; i++)
                terms[i].resize( i + 1);

        }




        std::array<double, 6> rhohats;
        VReal vr(P);
        VImag vi(P);
        WReal wr(P);
        WImag wi(P);
        XReal xr(P);
        XImag xi(P);

        for (int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
            {
                if (m == 0)
                {
                    double n_ = double(n);
                    vr.reset(m, n);
                    rhohats[0] = L2InnerProductDiscrete(data, &vr,system,*consts) / ((2 * n_ + 1) * (n_ + 1));

                    vi.reset(m, n);
                    rhohats[1] = 0.0;

                    wr.reset(m, n);
                    rhohats[2] = L2InnerProductDiscrete(data, &wr, system , *consts) / ((2 * n_ + 1) * n_);

                    wi.reset(m, n);
                    rhohats[3] = 0.0;

                    xr.reset(m, n);
                    rhohats[4] = L2InnerProductDiscrete(data, &xr, system , *consts) / (double)((n_ + 1) * n_);

                    xi.reset(m, n);
                    rhohats[5] = 0.0;

                    terms[n][m] = VSHTerm(P,m, n, rhohats);
                }
                else if (m > 0)
                {
                    double n_ = double(n);
                    vr.reset(m, n);
                    rhohats[0] = 2.0 * L2InnerProductDiscrete(data, &vr, system , *consts) / (double)((2 * n_ + 1) * (n_ + 1));
                    vi.reset(m, n);
                    rhohats[1] = 2.0 * L2InnerProductDiscrete(data, &vi, system , *consts) / (double)((2 * n_ + 1) * (n_ + 1));

                    wr.reset(m, n);
                    rhohats[2] = 2.0 * L2InnerProductDiscrete(data, &wr, system , *consts) / (double)((2 * n_ + 1) * n_);

                    wi.reset(m, n);
                    rhohats[3] = 2.0 * L2InnerProductDiscrete(data, &wi, system , *consts) / (double)((2 * n_ + 1) * n_);
                    xr.reset(m, n);
                    rhohats[4] = 2.0 * L2InnerProductDiscrete(data, &xr, system, *consts) / (double)((n_ + 1) * n_);

                    xi.reset(m, n);
                    rhohats[5] = 2.0 * L2InnerProductDiscrete(data, &xi, system , *consts) / (double)((n_ + 1) * n_);

                    terms[n][m] = VSHTerm(P,m, n, rhohats);
                }
                if (m < 0)
                    std::cout << "Error, m < 0 in VSH approximation function!\n";
                append(&terms[n][m]);
            }
    }

    /**
     * @brief Adds two VSHSeries.
     * @param s the second series to add.
    */
    VSHSeries plus(const VSHSeries & s , Legendre*P) const
    {
        VSHSeries temp(N, consts,P, s.system);

        temp.terms.resize(N + 1);

        for (int i = 0; i <= N; i++)
            temp.terms[i].resize(i + 1);



        for(int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
            {
                temp.terms[n][m] = terms[n][m].plus( s.terms[n][m] , P);
            }
        return temp;
            
    }

    /**
     * @brief subtracts two VSHSeries.
     * @param s the second series to add.
    */
    VSHSeries minus(const VSHSeries& s , Legendre*P) const
    {
        VSHSeries temp(N, consts,P, s.system);

        temp.terms.resize(N + 1);

        for (int i = 0; i <= N; i++)
            temp.terms[i].resize(i + 1);



        for (int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
            {
                temp.terms[n][m] = terms[n][m].minus( s.terms[n][m] , P);
            }
        return temp;

    }

    void axpy(const VSHSeries& s , const double & alpha)
    {



        for (int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
                terms[n][m].axpy(s.terms[n][m] , alpha);
            

    }

    VSHSeries uminus(Legendre*P) const
    {
        VSHSeries temp(N, consts, P, system);

        temp.terms.resize(N + 1);

        for (int i = 0; i <= N; i++)
            temp.terms[i].resize(i + 1);



        for (int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
            {
                temp.terms[n][m] = terms[n][m].uminus(P);
            }
        return temp;

    }

    VSHSeries times(const double & a , Legendre*P) const
    {
        VSHSeries temp(N, consts,P, system);

        temp.terms.resize(N + 1);

        for (int i = 0; i <= N; i++)
            temp.terms[i].resize(i + 1);



        for (int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
            {
                temp.terms[n][m] = terms[n][m].times( a , P);
            }
        return temp;

    }

    VSHSeries& operator +=(const VSHSeries& s)
    {
        for (int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
                terms[n][m] += s.terms[n][m];

        return *this;
    }

    VSHSeries& operator *=(const double & a)
    {
        for (int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
                terms[n][m] *= a;

        return *this;
    }

    /**
     * @brief Computes the @f$ L^2 @f$ inner product of two series.
     * @param series the second series in the product.
    */
    double dot(const VSHSeries & series) const
    {
        double total = 0.0;

        for (int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
                total += terms[n][m].dot(series.terms[n][m]);

        return total;
    }

    /**
     * @brief evaluates the series. Recasts x in the spherical coordinate system of the series.
     * @return 
    */
    RectCoord operator()(const SphereCoord & x) const
    {

        SphereCoord temp = recast(x, system);

        P->populate(cos(temp.s.theta), N, N);

        return VectorFieldSum::operator()(temp);
    }


};
#pragma once
#include "SphericalHarmonics.h"

/**
 * @brief Single term in the SSHSeries class. Contains each Spherical Harmonic Basis function for as given index.
*/
class SSHTerm : public SphericalScalarFunction
{
private:

    int m;
    int n;
    YReal yr;
    YImag yi;
    Legendre* poly;


public:
    double rhohatYr;
    double rhohatYi;



    SSHTerm(Legendre* P = nullptr) :yr(P), yi(P) { poly = P; m = 0; n = 0; }

    /**
     * @brief Constructs a term in the series given the normalized coefficients.
     * @param rhohats in order, the VReal, VImag, WReal, WImag, XReal, XImag coefficients to be assigned.
    */
    SSHTerm(Legendre*P,int m0, int n0, double rhyr, double rhyi) :
        yr(P,m0, n0),
        yi(P,m0, n0)
    {
        m = m0;
        n = n0;

        rhohatYr = rhyr;
        rhohatYi = rhyi;

        poly = P;

    }

    SSHTerm(const SSHTerm& term , Legendre*P) :
        yr(P,term.m, term.n),
        yi(P,term.m, term.n)
    {
        m = term.m;
        n = term.n;

        poly = P;

        rhohatYr = term.rhohatYr;
        rhohatYi = term.rhohatYi;
    }

    SSHTerm operator =(const SSHTerm& term)
    {
        m = term.m;
        n = term.n;

        poly = term.poly;
        yr = YReal(term.poly, m, n);
        yi = YImag(term.poly, m, n);

        rhohatYr = term.rhohatYr;
        rhohatYi = term.rhohatYi;


        return *this;
    }

    /**
     * @brief Addition operator between terms.
     * @param t The term to add.
     * @return A term whose coefficients are sums of the two terms.
    */
    SSHTerm plus(const SSHTerm& t , Legendre*P) const
    {
        double temprhyr;
        double temprhyi;

        temprhyr = rhohatYr + t.rhohatYr;
        temprhyi = rhohatYi + t.rhohatYi;

        return SSHTerm(P,m, n, temprhyr , temprhyi);
    }

    /**
     * @brief Subtraction operator between terms.
     * @param t The term to add.
     * @return A term whose coefficients are the difference of the two terms.
    */
    SSHTerm minus(const SSHTerm& t , Legendre*P) const
    {
        double temprhyr;
        double temprhyi;

        temprhyr = rhohatYr - t.rhohatYr;
        temprhyi = rhohatYi - t.rhohatYi;

        return SSHTerm(P,m, n, temprhyr, temprhyi);
    }

    /**
     * @brief Unry minus operator.
     * @return A term whose coefficients are the negative of the calling term.
    */
    SSHTerm uminus(Legendre*P) const
    {

        return SSHTerm(P,m, n, -rhohatYr, -rhohatYi);
    }


    SSHTerm times(const double& a, Legendre*P) const
    {

        return SSHTerm(P,m, n, rhohatYr * a, rhohatYi * a);
    }


    /**
     * @brief Computes the inner product between terms.
    */
    double dot(const SSHTerm& s) const
    {
        if (m != s.m || n != s.n)
            return 0.0;
        else
            return rhohatYr * s.rhohatYr + rhohatYi * s.rhohatYi;
    }

    void axpy(const SSHTerm& term, const double& alpha)
    {

        rhohatYr += term.rhohatYr * alpha;
        rhohatYi += term.rhohatYi * alpha;

    }

    SSHTerm& operator +=(const SSHTerm& term)
    {
        rhohatYr += term.rhohatYr;
        rhohatYi += term.rhohatYi;

        return *this;
    }

    SSHTerm& operator -=(const SSHTerm& term)
    {
        rhohatYr -= term.rhohatYr;
        rhohatYi -= term.rhohatYi;

        return *this;
    }

    SSHTerm& operator *=(const double& a)
    {
        rhohatYr *= a;
        rhohatYi *= a;

        return *this;
    }

    /**
     * @brief Evaluates the term at the point x.
     * @return In terms of complex numbers, @f$   \hat{\rho}_{n,m}^{V*}\boldsymbol{V}_n^m(\theta, \phi) +  \hat{\rho}_{n,m}^{W*}\boldsymbol{W}_n^m(\theta, \phi)
    +\hat{\rho}_{n,m}^{X*}\boldsymbol{X}_n^m(\theta, \phi)  @f$. The evaluation uses real functions.
    */
    double operator()(const SphereCoord& x) const
    {
        if (m == 0)
            return rhohatYr * yr(x);
        else
            return rhohatYr * yr(x) + rhohatYi * yi(x);
    }

    double dTheta(const SphereCoord& x) const
    {
        if (m == 0)
            return rhohatYr * yr.dTheta(x);
        else
            return rhohatYr * yr.dTheta(x) + rhohatYi * yi.dTheta(x);
    }

    double dPhi(const SphereCoord& x) const
    {
        if (m == 0)
            return rhohatYr * yr.dPhi(x);
        else
            return rhohatYr * yr.dPhi(x) + rhohatYi * yi.dPhi(x);
    }

};

/**
 * @brief The Vector Spherical Harmonic Series Functor. Computes and evaluates the fourier series representation of a SphericalVectorField.
*/
class SSHSeries : public SphereScalFunctionSum
{



    
    


public:
    Legendre* P;
    QuadConstants* consts;
    int N;
    SphereCoordSystem* system;
    std::vector<std::vector<SSHTerm>> terms;



    /**
     * @brief Constructs the series. n is the maximum value of the degree. there will be @f$ n(n+1)/2 @f$ terms total.
     * @param c The center of the coordinate system of the series.
    */
    SSHSeries(int n, QuadConstants* cts,Legendre*p,SphereCoordSystem* sys)
    {
        N = n;

        terms.resize(n + 1);

        for (int i = 0; i <= n; i++)
            terms[i].resize(i + 1);

        for (int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
                append(&terms[n][m]);

        system = sys;

        consts = cts;

        P = p;


    }

    SSHSeries(const SSHSeries& ssh)
    {
        N = ssh.N;

        terms = ssh.terms;

        system = ssh.system;

        consts = ssh.consts;

        
        for (int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
                append(&terms[n][m]);

        P = ssh.P;

    }

    SSHSeries(const SSHSeries& ssh , Legendre*p)
    {
        N = ssh.N;

        terms.resize(N + 1);

        for (int i = 0; i <= N; i++)
            terms[i].resize(i + 1);

        for (int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
            {
                terms[n][m] = SSHTerm(ssh.terms[n][m], p);
                append(&terms[n][m]);
            }


        system = ssh.system;

        consts = ssh.consts;

        P = p;

    }

    SSHSeries(const SSHSeries& ssh, QuadConstants* cts , Legendre*p , SphereCoordSystem* sys)
    {
        N = ssh.N;

        terms.resize(N + 1);

        for (int i = 0; i <= N; i++)
            terms[i].resize(i + 1);

        for (int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
            {
                terms[n][m] = SSHTerm(ssh.terms[n][m], p);
                append(&terms[n][m]);
            }


        system = sys;

        consts = cts;


        P = p;

    }


    SSHSeries& operator =(const SSHSeries& ssh)
    {
        N = ssh.N;

        terms = ssh.terms;

        system = ssh.system;

        consts = ssh.consts;
        P = ssh.P;

        clear();
        for (int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
                append(&terms[n][m]);

        return *this;

    }


    /**
     * @brief Approximates the passed SphericalVectorField.
     * @param r the function to approximate.
     * @param n number of terms in the series.
    */
    void approximate(SphericalScalarFunction* r, int p = 0)
    {
        clear();
        if (p >= 0)
        {
            N = p;
            terms.resize(p + 1);

            for (int i = 0; i <= p; i++)
                terms[i].resize(i + 1);

        }

        double rhohatyr;
        double rhohatyi;
        YReal yr(P);
        YImag yi(P);

        for (int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
            {
            

                    if (m == 0)
                    {
                        yr.reset(m, n);
                        rhohatyr = 2.0 * L2InnerProduct(r, &yr, system, *consts);
                        rhohatyi = 0.0;
                    }
                    else
                    {
                        yr.reset(m, n);
                        rhohatyr = 2.0 * L2InnerProduct(r, &yr, system, *consts);
                        yi.reset(m, n);
                        rhohatyi =  2.0 * L2InnerProduct(r, &yi, system, *consts);
                    }
          
                    terms[n][m] = SSHTerm(P,m, n, rhohatyr, rhohatyi);
                

                append(&terms[n][m]);
            }
    }

    void generateRandomSmooth(int p = 0)
    {
        clear();
        if (p >= 0)
        {
            N = p;
            terms.resize(p + 1);

            for (int i = 0; i <= p; i++)
                terms[i].resize(i + 1);

        }


        std::default_random_engine gen;
        std::uniform_real_distribution<double> r(0.0, 1.0);

        double rhohatyr;
        double rhohatyi;
        YReal yr(P);
        YImag yi(P);

        for (int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
            {
                
                if (m == 0)
                {
                    yr.reset(m, n);
                    rhohatyr =  r(gen) * std::exp(-n);
                    rhohatyi = 0.0;
                }
                else
                {
                    rhohatyr =  r(gen) * std::exp(-n);
                    rhohatyi =  r(gen) * std::exp(-n);
                }

                terms[n][m] = SSHTerm(P,m, n, rhohatyr, rhohatyi);


                append(&terms[n][m]);
            }
    }


    void approximate(ScalSphereData& data, int p = 0)
    {
        clear();
        if (p >= 0)
        {
            N = p;
            terms.resize(p + 1);

            for (int i = 0; i <= p; i++)
                terms[i].resize(i + 1);

        }




        double rhohatyr;
        double rhohatyi;
        YReal yr(P);
        YImag yi(P);

        for (int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
            {

                if (m == 0)
                {
                    yr.reset(m, n);
                    rhohatyr = L2InnerProductDiscrete(data, &yr, system, *consts);
                    rhohatyi = 0.0;
                }
                else
                {
                    yr.reset(m, n);
                    rhohatyr =   2.0 * L2InnerProductDiscrete(data, &yr, system, *consts);
                    yi.reset(m, n);
                    rhohatyi =  2.0 * L2InnerProductDiscrete(data, &yi, system, *consts);
                }

                terms[n][m] = SSHTerm(P, m, n, rhohatyr, rhohatyi);


                append(&terms[n][m]);
            }
    }


    /**
     * @brief Adds two VSHSeries.
     * @param s the second series to add.
    */
    SSHSeries plus(const SSHSeries& s , Legendre*P) const
    {
        SSHSeries temp(N, consts,P, system);

        temp.terms.resize(N + 1);

        for (int i = 0; i <= N; i++)
            temp.terms[i].resize(i + 1);



        for (int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
            {
                temp.terms[n][m] = terms[n][m].plus(s.terms[n][m] , P);
                temp.append(&temp.terms[n][m]);
            }
        return temp;

    }

    /**
     * @brief subtracts two VSHSeries.
     * @param s the second series to add.
    */
    SSHSeries minus(const SSHSeries& s , Legendre*P) const
    {
        SSHSeries temp(N, consts,P, system);

        temp.terms.resize(N + 1);

        for (int i = 0; i <= N; i++)
            temp.terms[i].resize(i + 1);



        for (int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
            {
                temp.terms[n][m] = terms[n][m].minus( s.terms[n][m] , P);
                temp.append(&temp.terms[n][m]);
            }
        return temp;

    }

    void axpy(const SSHSeries& s, const double& alpha)
    {



        for (int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
                terms[n][m].axpy(s.terms[n][m], alpha);


    }

    SSHSeries uminus(Legendre*P) const
    {
        SSHSeries temp(N, consts,P, system);

        temp.terms.resize(N + 1);

        for (int i = 0; i <= N; i++)
            temp.terms[i].resize(i + 1);



        for (int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
            {
                temp.terms[n][m] = terms[n][m].uminus(P);
                temp.append(&temp.terms[n][m]);
            }
        return temp;

    }

    SSHSeries times(const double& a , Legendre*P) const
    {
        SSHSeries temp(N, consts,P, system);

        temp.terms.resize(N + 1);

        for (int i = 0; i <= N; i++)
            temp.terms[i].resize(i + 1);



        for (int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
            {
                temp.terms[n][m] = terms[n][m].times( a , P);
                temp.append(&temp.terms[n][m]);
            }
        return temp;

    }

    SSHSeries& operator +=(const SSHSeries& s)
    {
        for (int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
                terms[n][m] += s.terms[n][m];

        return *this;
    }

    SSHSeries& operator -=(const SSHSeries& s)
    {
        for (int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
                terms[n][m] -= s.terms[n][m];

        return *this;
    }


    SSHSeries& operator *=(const double& a)
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
    double dot(const SSHSeries& series) const
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
    double operator()(const SphereCoord& x) const
    {

        SphereCoord temp = recast(x, system);

        P->populate(cos(temp.s.theta), N, N);

        return SphereScalFunctionSum::operator()(temp);
    }

    double dTheta(const SphereCoord& x) const
    {
        SphereCoord xtemp = recast(x, system);

        P->populate(cos(xtemp.s.theta), N+1, N+1);
        double temp = 0.0;
        for (int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
                temp += terms[n][m].dTheta(xtemp);

        return temp;
    }

    double dPhi(const SphereCoord& x) const
    {
        SphereCoord xtemp = recast(x, system);

        P->populate(cos(xtemp.s.theta), N + 1, N + 1);
        double temp = 0.0;
        for (int n = 0; n <= N; n++)
            for (int m = 0; m <= n; m++)
                temp += terms[n][m].dPhi(xtemp);

        return temp;
    }


};

class SurfaceGrad : public SphericalVectorField
{
private:
    SSHSeries* f;
public:
    SurfaceGrad(SSHSeries* rho): f(rho){}

    RectCoord operator()(const SphereCoord& x) const
    {
			SphereCoord temp = recast(x, f->system);
            return f->dTheta(temp) * e_theta(temp) + f->dPhi(temp) * e_phi(temp) / sin(x.s.theta);
    }

};
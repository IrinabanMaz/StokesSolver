/**
* @file Legendre.h
* 
* Contains the class definition for the Associated Legendre Functions.
* We compute the values of the funtions and their derivatives by recurrance relations.
* Once computed at a point, the values of the functions are stored for all indices until the point is changed.
* 
*/



#pragma once
#include <vector>
#include<cmath>
#include<iostream>
#include<map>



/**
*
* @brief computes the factorial of the argument.
*/
double factorial(int n)
{
    static int maxn = 0;
    static std::vector<double> vals;
    static bool firstcall = true;

    if (firstcall)
    {
        vals.resize(1);
        vals[0] = 1;
        firstcall = false;
    }
    if (n < 0)
    {
        std::cout << "Error: negative integer passed to factorial";
        return -1;
    }
    else if (n <= maxn)
        return vals[n];
    else
    {
        vals.resize(n + 1);
        for (int i = maxn + 1; i <= n; i++)
            vals[i] = vals[i - 1] * i;

        maxn = n;
        return vals[n];
    }
}

double ffactorial(int n)
{
    static int maxn = 1;
    static std::vector<double> vals;
    static bool firstcall = true;

    if (firstcall)
    {
        vals.resize(2);
        vals[0] = 1;
        vals[1] = 1;
        firstcall = false;
    }
    if (n < 0)
    {
        std::cout << "Error: negative integer passed to factorial";
        return -1;
    }
    else if (n <= maxn)
        return vals[n];
    else
    {
        vals.resize(n + 1);
        for (int i = maxn + 1; i <= n; i++)
            vals[i] = vals[i - 2] * i;

        maxn = n;
        return vals[n];
    }
}
/**
* @brief comparison operator for maps with double. Returns true if a < b within relative accuracy of template parameter.
*/
struct LessWithTol
{
    bool operator()(const double& a, const double& b) const
    {
        
            return a  < b - 1e-14 * abs(a);
    }
};

inline bool EqualsWithTol(const double& a, const double& b)
{
    LessWithTol lwt;

    return !lwt(a, b) && !lwt(b, a);
}

typedef std::map<double, std::vector<std::vector<double>>, LessWithTol> PrecomputedVals;

/**
*  @brief Class representing the Associated Legendre Functions.
*/
class Legendre
{


private:
    //the value of the input variable to the function.
    double input; /**< Value of the last argument passed to the function. */
    //the maximum value of the parameter m computed for this value of the input.
    int maxorder; /**< Value of the highest order computed for the value above. */
    //the maximum value of the degree n computed for this value of the input.
    int maxdegree; /**< Value of the highest degree computed for the value above. */
    std::vector<std::vector<double>> values; /**< Values of the Legendre Functions for the value input. */
    //PrecomputedVals prevals;

public:

    Legendre()
    {
        input = 100;
        maxorder = -1;
        maxdegree = -1;
    }

    /** 
       \fn void populate(double val , int M, int L)
       \brief Computes the valeus of the Legendre Functions at the value val.

       The function computes the Legendre functions using the following recurrence relations.
       \f{eqnarray}{
        P_0^0(x) &=& 1 \\
        P_1^0(x) &=& x \\
        P_{n+1}^0(x) &=&  \frac{2 n + 1}{n+1} x P_n^0(x) -  \frac{n}{n+1} P_{n-1}^0(x) \\
        P_n^{m+1}(x) &=& (n-m) \frac{x}{1-x^2} P_n^m(x) - (n+m)\frac{1}{1-x^2}P_{n-1}^m(x) 
       \f}
    */
    void populate(double val, int M, int L)
    {
        if (abs(val - input) < 1e-14)
            if (maxdegree > L)
                if (maxorder > L)
                    return;
/*
        if(prevals.count(val) > 0)
            if (maxdegree > L)
                if (maxorder > L)
                {
                    input = val;
                    values = prevals[input];
                    return;
                }
  */      

        input = val;
        maxorder = std::max(L, 2);
        maxdegree = std::max(L , 2);
        double s = sqrt(1 - input * input);

        


        if (values.size() < maxdegree+1) 
        {
            values.resize(maxdegree + 1);
            for (int l = 0; l <= maxdegree; l++)
                    values[l].resize(l + 1);
                
            
        }

        
        //set values in case of boundary condition.

        //set initial values.
        values[0][0] = 1.0;
        values[1][0] = input;

        //set values for m = 0 using recurrence in terms of two previous i values.
        for (int m = 1; m <= maxdegree; m++)
        {
            values[m][m] = (-2.0 * (m%2) + 1.0) * ffactorial(2 * m - 1) * pow(s,m);
        }

        for (int m = 0; m < maxdegree; m++)
        {
                values[m+1][m] = input * (2.0 * (double)m + 1.0) * values[m][m];
        }

        //computes values for larger values of m in terms of previous value of m. for i = 0, uses the values from i = 1, i = 2, instead.
        for (int n = 2; n <=maxdegree; n++)
           for (int m = n-1; m > 0; m--)
                values[n][m - 1] = -(values[n][m+1] + 2.0 * (double)m * input / s * values[n][m])/(((double)n - (double)m + 1.0) * ((double)n + (double)m));

            
      //  prevals[val] = values;
        
        
    }

private:
    /** 
    \fn double getValue(int m, int i)
    \brief Retrieves the values computed by populate at the given indices,
    */
    double getValue(double x, int m, int i)
    {/*
        std::vector<std::vector<double>>* vals;
        
        if (EqualsWithTol(x , input))
            vals = &values;
        else
            vals = &prevals[x];
            */
        if (m >= 0)
            return values[i][m];
        else 
            return ((m % 2 == 0) ? 1.0 : -1.0) * (double)factorial(i + m) / (double)factorial(i - m) * values[i][-m];

        
    }
public:

    /**
    \fn double operator()(double x, int m, int l)

     Performs populate(double val, int M, int L) and getValue(int m, int i) in order for the arguments.
    */
    double operator()(double x, int m, int l)
    {
        
        populate(x, m, l);

        return getValue(x , m, l);
    }

    /**
    \fn double dTheta(int m, int n, double theta, int order = 1)
    \brief Computes the first or second derivative with respect to \f$ \theta \f$ of the Legendre Function at the value \f$ \cos(\theta) \f$.
    
    For the first derivative, we use the recurrence relations
    \f{eqnarray}{
    \frac{dP_n^m(x)}{dx} &=& -(n + m)(n - m + 1))P_n^{m-1}(x)
            -m \frac{x}{\sqrt{1-x^2}}P_n^m(x)\\
    \frac{dP_n^m(x)}{dx} &=&   P_n^{m+1}(x) + m  \frac{x}{1-x^2}  P_n^m(x)
    \f}
    */
    double dTheta(int m, int n, double theta, int order = 1)
    {
        double x = cos(theta);
        populate(x, m +1, n+1 );

        double temp = 0.0;
        if (order == 0)
            return getValue(x , m, n);
        else if(order == 1)
        {
            if (n == 0)
                temp = 0;
            else if (m >= 1)
                temp = -(double)((n + m) * (n - m + 1)) * getValue(x ,m - 1, n)
                - (double)m * x / sin(theta) * getValue(x , m, n);
            else
                temp = getValue(x , m + 1, n) + m * x / sin(theta) * getValue(x , m, n);


            return temp;
        }
        else if (order == 2)
        {
            if (n == 0)
                return 0;
            else if (m >= 0)
            {
                double y = sin(theta);
                return -(double)((n + m) * (n - m + 1)) * dTheta(m - 1, n, theta, 1)+ (double)m/ (y * y) * getValue(x ,m, n)
                    - (double)m * x / y * dTheta(m, n, theta, 1);
            }
            else
            {
                double y = sin(theta);
                return -(double)(n + m) * (n - m + 1) * dTheta(m + 1, n, theta, 1) - (double)m * x / (y * y) * getValue(x ,m, n)
                    + (double)m * x / y * dTheta(m, n, theta, 1);
            }
        }
     return 0.0;
    }
   

   

};
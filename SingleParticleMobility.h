#pragma once
#include "StokesOperators.h"
#include "GMRES.h"

const int MAXITER = 100;
int numcount = 0;
class StokesParticle
{
public:
	int numSeriesTerms;

	SphereCoordSystem system;
	Legendre P;
	
	SphereData data;
	VSHSeries fourierdata;


	QuadConstants consts;
	SingleLayerDirectDiscrete quadflow;
	StokesTractionDirectDiscrete quadtraction;

	SingleLayerHarmonic flow;
	StokesPressure pressure;
	StokesTractionHarmonic traction;
	

	StokesParticle(int n = 6 , SphereCoordSystem sys = defaultsystem  ) :
		  data(2 * n + 1,n + 1),
		  consts(2 * n+1 , n + 1),
		  fourierdata(n, &consts,&P, &system),
		  quadflow(&data,&consts, &system),
		  quadtraction(&data,&consts , &system),
		  flow(fourierdata,&P) ,
		  pressure(fourierdata,&P) , 
		  traction(&flow , &pressure),
		  system(sys),
		  numSeriesTerms(n)
	{
		fourierdata.approximate(data, numSeriesTerms);
		flow.solve(fourierdata);
		pressure.solve(fourierdata);
	}

	StokesParticle(const SphereData& d, SphereCoordSystem sys, int nst) :
		data(d),
		consts(2 * nst + 1, nst + 1),
		numSeriesTerms(nst),
		fourierdata(nst, &consts, &P, &system),
		quadflow(&data, &consts, &system),
		quadtraction(&data, &consts, &system),
		flow(fourierdata,&P),
		pressure(fourierdata ,&P),
		traction(&flow, &pressure),
		system(sys)
	{
		fourierdata.approximate(data, numSeriesTerms);
		flow.solve(fourierdata);
		pressure.solve(fourierdata);
	}


	

	StokesParticle(const StokesParticle& particle) :
		numSeriesTerms(particle.numSeriesTerms) ,
		data(particle.data),
		consts(2*particle.numSeriesTerms + 1 , particle.numSeriesTerms + 1),
		quadflow(&data,&consts, &system),
		quadtraction(&data,&consts, &system),
		fourierdata(particle.fourierdata , &consts , &P , &system),
		flow(fourierdata , &P), 
		pressure(fourierdata , &P),
		traction(&flow, &pressure),
		system(particle.system)
	{

	}

	StokesParticle operator +(const StokesParticle & particle) const
	{
		StokesParticle temp(numSeriesTerms , system);

		
			for (int t = 0; t < (int)consts.NUMGLNODES; t++)
				for (int p = 0; p < (int)consts.NUMTRAPNODES; p++)
					temp.data[t][p] = data[t][p] + particle.data[t][p];
		

		temp.fourierdata = fourierdata.plus(particle.fourierdata , &temp.P);


		temp.numSeriesTerms = numSeriesTerms;

		temp.flow.solve(temp.fourierdata);
		temp.pressure.solve(temp.fourierdata);

		temp.system = system;

		return temp;

	}

    StokesParticle operator -(const StokesParticle & particle) const
	{
		StokesParticle temp(numSeriesTerms , system);

		for (int t = 0; t < consts.NUMGLNODES; t++)
			for (int p = 0; p < consts.NUMTRAPNODES; p++)
				temp.data[t][p] = data[t][p] - particle.data[t][p];

		temp.fourierdata = fourierdata.minus( particle.fourierdata , &temp.P);

		temp.numSeriesTerms = numSeriesTerms;

		temp.flow.solve(temp.fourierdata);
		temp.pressure.solve(temp.fourierdata);

		temp.system = system;
		return temp;

	}

	StokesParticle operator -() const
	{
		StokesParticle temp(numSeriesTerms , system);

		for (int t = 0; t < consts.NUMGLNODES; t++)
			for (int p = 0; p < consts.NUMTRAPNODES; p++)
				temp.data[t][p] = - data[t][p];

		temp.fourierdata = fourierdata.uminus(&temp.P);

		temp.numSeriesTerms = numSeriesTerms;

		temp.flow.solve(temp.fourierdata);
		temp.pressure.solve(temp.fourierdata);

		temp.system = system;
		return temp;
	}

    StokesParticle & operator +=(const StokesParticle & particle)
	{
		for (int t = 0; t < (int)consts.NUMGLNODES; t++)
			for (int p = 0; p < (int)consts.NUMTRAPNODES; p++)
				data[t][p] += particle.data[t][p];


		fourierdata += particle.fourierdata;



		flow.solve(fourierdata);
		pressure.solve(fourierdata);
		return *this;
	}

	StokesParticle & operator *=(const double & a)
	{
		for (int t = 0; t < (int)consts.NUMGLNODES; t++)
			for (int p = 0; p < (int)consts.NUMTRAPNODES; p++)
				data[t][p] *= a;


		fourierdata *= a;



		flow.solve(fourierdata);
		pressure.solve(fourierdata);
		return *this;
	}

	StokesParticle operator *(const double& a) const
	{
		StokesParticle temp(numSeriesTerms , system);
		for (int t = 0; t < (int)consts.NUMGLNODES; t++)
			for (int p = 0; p < (int)consts.NUMTRAPNODES; p++)
				temp.data[t][p] = data[t][p] * a;

		temp.fourierdata = fourierdata.times( a , &temp.P);

		temp.numSeriesTerms = numSeriesTerms;

		temp.flow.solve(temp.fourierdata);
		temp.pressure.solve(temp.fourierdata);

		temp.system = system;

		return temp;
	}

	StokesParticle& operator =(const StokesParticle & particle)
	{
		data = particle.data;
		fourierdata = VSHSeries(particle.fourierdata , &consts , &P , &system);
		numSeriesTerms = particle.numSeriesTerms;
		consts = particle.consts;
		flow.solve(fourierdata);
		pressure.solve(fourierdata);
		system = particle.system;

		return *this;

	}



	double norm() const
	{
		return sqrt(dot(*this));
	}

	double dot(const StokesParticle& particle) const
	{
		return L2InnerProductDiscrete(data , particle.data , system.radius , consts );
	}

	static double dot(const StokesParticle &particle, const StokesParticle& qarticle)
	{
		return particle.dot(qarticle);
	}

	int getN() const
	{
		return (int)consts.NUMGLNODES * (int)consts.NUMTRAPNODES * 3;
	}

	void axpy(StokesParticle* particle, double scale)
	{


		for (int t = 0; t < (int)consts.NUMGLNODES; t++)
			for (int p = 0; p < (int)consts.NUMTRAPNODES; p++)
				data[t][p] +=particle->data[t][p] * scale;


		fourierdata.axpy( particle->fourierdata , scale);



		flow.solve(fourierdata);
		pressure.solve(fourierdata);

	}

	

	void refreshData()
	{
		fourierdata.approximate(data, numSeriesTerms);
		

		flow.solve(fourierdata);
		pressure.solve(fourierdata);
	}

	
};

RectCoord rotationCoefficient(const SphereData& data ,SphereCoordSystem system, const QuadConstants & consts)
{
	

		RectCoord total = 0.0;

		QuadConstants temp;
		if (consts.NUMGLNODES != data.size() || consts.NUMTRAPNODES != data[0].size())
			temp = QuadConstants(data[0].size(), data.size());
		else
			temp = consts;

		for (int p = 0; p < temp.NUMTRAPNODES; p++)
			for (int i = 0; i < temp.NUMGLNODES; i++)
			{
				SurfaceCoord s(temp.GLnodes[i], 2.0 * temp.PI * (double)p / (double)temp.NUMTRAPNODES);
				SphereCoord x(s, &system);
				total = total + system.radius * system.radius * temp.GLweights[i] * cross(SphereToRect(x) - system.center,data[i][p]) * (2.0 * temp.PI / temp.NUMTRAPNODES);


			}

		return total;
	
}

class LHSOperator
{
public:
	StokesParticle operator *(const StokesParticle & particle) const
	{
		StokesParticle temp(particle.numSeriesTerms , particle.system);

		RectCoord average = 0.25 / particle.consts.PI * IntegrateDiscrete(particle.data ,&particle.system, particle.consts );
		RectCoord rcoef = 0.375/particle.consts.PI * rotationCoefficient(particle.data, particle.system, particle.consts);

		for(int t = 0; t < particle.consts.NUMGLNODES; t++)
			for (int p = 0; p < particle.consts.NUMTRAPNODES; p++)
			{
				SurfaceCoord s(particle.consts.GLnodes[t], 2.0 * particle.consts.PI * (double)p / (double)particle.consts.NUMTRAPNODES);
				SphereCoord x(s , &particle.system);

				temp.data[t][p] =  particle.data[t][p] + particle.traction(x) + average + cross(rcoef, x - x.system->center);

			}

		temp.refreshData();

		return temp;
	}


};

#include <iomanip>
class RHSOperator
{
public:
	StokesParticle operator *(const StokesParticle & particle)
	{
	    StokesParticle temp(particle.numSeriesTerms , particle.system);


		for (int t = 0; t < particle.consts.NUMGLNODES; t++)
		{
			for (int p = 0; p < particle.consts.NUMTRAPNODES; p++)
			{
				SurfaceCoord s(particle.consts.GLnodes[t], 2.0 * particle.consts.PI * (double)p / (double)particle.consts.NUMTRAPNODES);
				SphereCoord x(s , &particle.system);

				temp.data[t][p] =  -particle.data[t][p] - particle.traction(x);

			}
		}

		temp.refreshData();
/*
		std::cout << "Real part of Rhohat_0^1V " << std::setprecision(16) << temp.fourierdata.terms[1][0].rhohatVr << "\n";
		std::cout << "Imag part of Rhohat_0^1V " << std::setprecision(16) << temp.fourierdata.terms[1][0].rhohatVi << "\n";
		std::cout << "Real part of Rhohat_0^1W " << std::setprecision(16) << temp.fourierdata.terms[1][0].rhohatWr << "\n";
		std::cout << "Imag part of Rhohat_0^1W " << std::setprecision(16) << temp.fourierdata.terms[1][0].rhohatWi << "\n";
		std::cout << "Real part of Rhohat_0^1X " << std::setprecision(16) << temp.fourierdata.terms[1][0].rhohatXr << "\n";
		std::cout << "Imag part of Rhohat_0^1X " << std::setprecision(16) << temp.fourierdata.terms[1][0].rhohatWi << "\n";
*/	

	return temp;
	}
};

class StokesIdentityPreconditioner
{
public:
	StokesParticle solve(const StokesParticle& particle)
	{
		return particle;
	}
};


class ForceBalance : public SphericalVectorField
{
private:
	RectCoord force;
	RectCoord torque;
	double PI = 4.0 * atan(1.0);
public:

	ForceBalance(RectCoord f , RectCoord t) : force(f) , torque(t){}

	RectCoord operator()(const SphereCoord & x) const
	{
		return 0.25 / PI *  force + 0.375 / PI * cross(torque , x - x.system->center);
	}
};



StokesParticle SolveMobilitySingleSphere(RectCoord F , RectCoord T , int seriesSize)
{

		
	
	LHSOperator L;
	ForceBalance rho(F,T);
	
	const double PI = 4.0 * atan(1.0);

	QuadConstants consts(2 * seriesSize + 1, seriesSize + 1);

	StokesParticle rh(discretize(&rho, &defaultsystem , consts ), defaultsystem, seriesSize);
	RHSOperator R;
	StokesParticle rhs = R * rh;
	
	StokesIdentityPreconditioner I;
    StokesParticle solution = -rh;
	
	GMRES(&L, &solution, &rhs, &I, 100, 5, 1e-10);


	
	StokesParticle soln = solution + rh;

	std::cout <<"Average of BIE solution: " << IntegrateDiscrete(solution.data) * 0.25 / PI << std::endl;
	std::cout <<"Average of rhs: " << IntegrateDiscrete(rh.data) * 0.25 / PI << std::endl;
	std::cout <<"Integral of flow: " <<Integrate(&soln.flow) * 0.25 / PI << std::endl;
	std::cout << "Integral of rho: " << Integrate(&rh.flow) * 0.25 / PI << std::endl;
	std::cout << "Angular Velocity of rho: " << rotationCoefficient(discretize(&rh.flow),rh.system , rh.consts) * 0.375 / PI << "\n";
	return soln;

}

StokesParticle SolveBIEWithKnown(int seriesSize)
{



	LHSOperator L;

	std::default_random_engine gen;
	std::uniform_real_distribution<double> rand;

	SphereData data(2 * seriesSize + 1,  seriesSize + 1);

	

	for (size_t t = 0; t < data.size(); t++)
	{
		for (size_t p = 0; p < data[t].size(); p++)
			data[t][p] = RectCoord(rand(gen), rand(gen), rand(gen));
	}
	StokesParticle rh(data, RectCoord(), seriesSize);
	
	StokesParticle rhs = L * rh;

	StokesIdentityPreconditioner I;
	StokesParticle solution;

	GMRES(&L, &solution, &rhs, &I, 100, 5, 1e-10);

	StokesParticle err = solution - rh;

	std::cout << "Error in solving: " << err.norm()<< std::endl;

	
	std::cout << "[ ";
	for (int t = 0; t < err.consts.NUMGLNODES - 1; t++)
	{
		std::cout << "[";
		for (int p = 0; p < err.consts.NUMTRAPNODES - 1; p++)
		  std::cout << norm(err.data[t][p]) << ", ";
		std::cout << norm(err.data[t][err.consts.NUMTRAPNODES - 1]) << " ], ";
	}
	std::cout << "[";
	for (int p = 0; p < err.consts.NUMTRAPNODES - 1; p++)
		std::cout << norm(err.data[(int)(err.consts.NUMTRAPNODES - 1)][p]) << ", ";
	std::cout << norm(err.data[(int)(err.consts.NUMGLNODES - 1)][err.consts.NUMTRAPNODES - 1]) <<" ]]" << std::endl;


	return solution;

}
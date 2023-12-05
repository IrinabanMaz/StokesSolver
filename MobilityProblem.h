#pragma once
#include "SingleParticleMobility.h"
#include <ctime>

class StokesParticleSystem : public SphericalVectorField
{
public:

	std::vector<StokesParticle> particles;
	VectorFieldSum traction;

	StokesParticleSystem(std::vector<SphereCoordSystem> systems = std::vector<SphereCoordSystem>(), int n = 6) 
	{
		particles.resize(systems.size());
		for (size_t i = 0; i < particles.size(); i++)
		{
			particles[i] = StokesParticle( n , systems[i]);

			traction.append(&particles[i].traction);
		}
	}

	StokesParticleSystem(std::vector<SphereData> ds, std::vector<SphereCoordSystem> systems, int nst)
	{

		if (ds.size() != systems.size())
			std::cout << "Error: Data size and centers size in paticle system are not equal.\n";

		particles.resize(systems.size());
		for (size_t i = 0; i < particles.size(); i++)
		{
			particles[i] = StokesParticle(ds[i], systems[i], nst);

			traction.append(&particles[i].traction);
		}

	
	}


	StokesParticleSystem(const StokesParticleSystem& ps)
	{

		traction.clear();

		particles.resize(ps.particles.size());
		for (int i = 0; i < particles.size(); i++)
		{
			particles[i] = ps.particles[i];
			traction.append(&particles[i].traction);
		}
	}


	StokesParticleSystem(const int& n)
	{
		particles.resize(n);
		for (int i = 0; i < particles.size(); i++)
		{
			traction.append(&particles[i].traction);
		}
	}

	void push_back(const StokesParticle & p)
	{
		particles.push_back(p);
		traction.clear();
		for (int i = 0; i < particles.size(); i++)
		{
			traction.append(&particles[i].traction);
		}
	}

	StokesParticleSystem operator +(const StokesParticleSystem& ps) const
	{
		StokesParticleSystem temp(*this);

		if (particles.size() != ps.particles.size())
			std::cout << "Error:Summing different numbers of particles.\n";

		for (size_t i = 0; i < particles.size(); i++)
			temp.particles[i] = particles[i] + ps.particles[i];

		return temp;

	}

	StokesParticleSystem operator -(const StokesParticleSystem& ps) const
	{
		StokesParticleSystem temp(*this);

		if (particles.size() != ps.particles.size())
			std::cout << "Error:Summing different numbers of particles.\n";

		for (size_t i = 0; i < particles.size(); i++)
			temp.particles[i] = particles[i] - ps.particles[i];

		return temp;

	}

	StokesParticleSystem operator -() const
	{
		StokesParticleSystem temp(*this);


		for (size_t i = 0; i < particles.size(); i++)
			temp.particles[i] = -particles[i];

		return temp;

	}

	StokesParticleSystem & operator +=(const StokesParticleSystem& ps)
	{
		for (int i = 0; i < particles.size(); i++)
			particles[i] += ps.particles[i];

		return *this;
	}

	StokesParticleSystem & operator *=(const double& a)
	{
		for (int i = 0; i < particles.size(); i++)
			particles[i] *= a;
		
		return *this;
	}

	StokesParticleSystem operator *(const double& a) const
	{
		StokesParticleSystem temp(*this);


		for (size_t i = 0; i < particles.size(); i++)
			temp.particles[i] *=a ;

		return temp;

	}

	StokesParticleSystem& operator =(const StokesParticleSystem& ps)
	{
		particles = ps.particles;

		traction.clear();

		for (int i = 0; i < particles.size(); i++)
		{
				traction.append(&particles[i].traction);
		}

		return *this;

	}

	


	double norm() const
	{
		double temp = 0.0;
		for (int i = 0; i < particles.size();i++)
		{
			double normj = particles[i].norm();
			temp += normj * normj;
		}
		return sqrt(temp);

	}

	double dot(const StokesParticleSystem& ps) const
	{
		double temp = 0.0;
		for (size_t i =0; i < particles.size(); i++)
		{
			temp += particles[i].dot(ps.particles[i]);
		}
		return temp;

	}

	static double dot(const StokesParticleSystem& particles, const StokesParticleSystem& qarticles)
	{
		return particles.dot(qarticles);
	}

	int getN() const
	{
		return particles.size();
	}

	void axpy(StokesParticleSystem* ps, double scale)
	{
		for (int i = 0; i < particles.size(); i++)
			particles[i].axpy(&(ps->particles[i]), scale);
		
	}


	RectCoord operator()(const SphereCoord& x) const
	{
		RectCoord temp;

		for (int i = 0; i < particles.size(); i++)
		{
			RectCoord r = x.system->center - particles[i].system.center;
			if ( sqrt(r.x * r.x + r.y * r.y + r.z * r.z) < 4.0)
				temp += particles[i].flow(x);
			else
				temp += particles[i].quadflow(x);
		}

		return temp;
	}

};


class PSLHSOperator
{
public:
	StokesParticleSystem operator *(const StokesParticleSystem& ps)
	{
		LHSOperator L;

		StokesParticleSystem temp(ps.getN());

		for (size_t i = 0; i < ps.particles.size(); i++)
		{
			std::time_t ti = std::time(0);
			temp.particles[i] = L * ps.particles[i];
			SphereCoordSystem system = ps.particles[i].system;
			StokesParticle tr( ps.particles[i].numSeriesTerms , system);
			
			SphereData t(2 * ps.particles[i].numSeriesTerms + 1 , ps.particles[i].numSeriesTerms + 1);

			

			for (size_t j = 0; j < ps.particles.size(); j++)
			{
				if (i != j)
				{
					if(norm(ps.particles[i].system.center - ps.particles[j].system.center) < 4.0)
						t = t + discretize(&ps.particles[j].traction, &system , ps.particles[j].consts);
					else
						t = t + discretize(&ps.particles[j].quadtraction, &system , ps.particles[j].consts);
				}
			}
			temp.particles[i] += StokesParticle(t , system , ps.particles[i].numSeriesTerms);
			std::cout << "Operator on particle " << i + 1 << " computed. Time elapsed: " << std::time(0) - ti <<"\n";
		}
		return temp;
	}
};

	class PSRHSOperator
	{
	public:
		StokesParticleSystem operator *(StokesParticleSystem& ps)
		{
			RHSOperator R;

			StokesParticleSystem temp(ps.getN());

			for (size_t i = 0; i < ps.particles.size(); i++)
			{
				temp.particles[i] = R * ps.particles[i];
				SphereCoordSystem system = ps.particles[i].system;
				StokesParticle tr(ps.particles[i].numSeriesTerms , system);


				SphereData t(2 * ps.particles[i].numSeriesTerms + 1, ps.particles[i].numSeriesTerms + 1);

				for (size_t j = 0; j < ps.particles.size(); j++)
				{
					if (i != j)
					{
						if (norm(ps.particles[i].system.center - ps.particles[j].system.center) < 4.0)
							t = t - discretize(&ps.particles[j].traction, &system , ps.particles[j].consts);
						else
							t = t - discretize(&ps.particles[j].quadtraction, &system, ps.particles[j].consts);
					}
				}
				temp.particles[i] += StokesParticle(t, system, ps.particles[i].numSeriesTerms);
				std::cout << "RHS Operator on particle " << i + 1 << " computed.\n";
			
			}
			return temp;
		}

		
	};
class SPSIdentityPreconditioner
{
public:
	StokesParticleSystem solve(const StokesParticleSystem& ps)
	{
		return ps;
	}
};


//Solve single particle problem with linear combinations of F,T, compare to exact solution off sphere.
// evaulate solution for single sphere problem off sphere, compare to exact solution off sphere.
//Two particles widely seperated, different forces and torques, compare solutions to isolated sphere problem.
// Show difference in data using sup norm.

StokesParticleSystem SolveMobilityManySphere(std::vector<SphereCoordSystem> systems,std::vector<RectCoord> Fs, std::vector<RectCoord> Ts, int seriesSize ,const StokesParticleSystem & initial= StokesParticleSystem(0), int numIters = 100, double r = 1e-9)
{

	if ((systems.size() != Fs.size()) || (Fs.size() != Ts.size()))
	{
		std::cout << "Error: size mismatch between mobility problem input vectors.\n";
		return StokesParticleSystem(0);
	}

	int size = systems.size();
	StokesParticleSystem BIEsolution(systems , seriesSize);
	if (initial.particles.size() == size)
		BIEsolution = initial;
	
		
	PSLHSOperator L;



    StokesParticleSystem rh(systems , seriesSize);
	
	for (int i = 0; i < size; i++)
	{
		ForceBalance rho = ForceBalance(Fs[i], Ts[i]);
		rh.particles[i] = StokesParticle(discretize(&rho , &systems[i] , rh.particles[i].consts), systems[i], seriesSize);
	}


	BIEsolution = -rh;
	PSRHSOperator R;
	StokesParticleSystem rhs = R * rh;

	SPSIdentityPreconditioner I;

	GMRES(&L, &BIEsolution, &rhs, &I, numIters, 1, r);



	StokesParticleSystem soln = BIEsolution + rh;

	//std::cout << Integrate(&BIEsolution, 1.0, BIEsolution.particles[0].center) * 0.25 / PI << std::endl;
	//std::cout << Integrate(&rh, 1.0, rh.particles[0].center) * 0.25 / PI << std::endl;

	return soln;

}
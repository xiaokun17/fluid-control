#ifndef __TimeStepPBF_h__
#define __TimeStepPBF_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/TimeStep.h"
#include "SimulationDataPBF.h"
#include "SPlisHSPlasH/SPHKernels.h"

namespace SPH
{
	class SimulationDataPBF;


	/** \brief This class implements the position-based fluids approach introduced 
	* by Macklin and Mueller \cite Macklin:2013:PBF, \cite BMOTM2014, \cite BMM2015.
	*/
	class TimeStepPBF : public TimeStep
	{
	protected:
		SimulationDataPBF m_simulationData;
		unsigned int m_counter;
		int m_velocityUpdateMethod;
		Real alpha;
		Real beta;
		Real gamma;
		bool enableTarget;
		/** Perform a position-based correction step for the following density constraint:\n
		*  \f$C(\mathbf{x}) = \left (\frac{\rho_i}{\rho_0} - 1 \right )= 0\f$\n
		*/
		void pressureSolve();
		void pressureSolveIteration(const unsigned int fluidModelIndex, Real &avg_density_err);

		//control
		void controlPressureSolve();
		void computeLambda(const unsigned int controlModelIndex);
		void controlConstraints();
		void updateSpring(const unsigned int fluidModelIndex);
		void updateDensityControl(const unsigned int fluidModelIndex);
		void updateVelocityControl(const unsigned int fluidModelIndex);

		void goToTarget();

		/** Perform the neighborhood search for all fluid particles. 
		*/
		void performNeighborhoodSearch();

		virtual void emittedParticles(FluidModel *model, const unsigned int startIndex);

		virtual void initParameters();

	public:
		static int VELOCITY_UPDATE_METHOD;
		static int ENUM_PBF_FIRST_ORDER;
		static int ENUM_PBF_SECOND_ORDER;
		static int DENSITY_ALPHA;
		static int SPRING_BETA;
		static int VELOCITY_GAMMA;
		static int ENABLE_TARGET;


		/** Initialize the simulation data required for this method. */
		TimeStepPBF();
		virtual ~TimeStepPBF(void);

		/** Perform a simulation step. */
		virtual void step();

		/** Reset the simulation method. */
		virtual void reset();
		virtual void resize();
	};
}

#endif

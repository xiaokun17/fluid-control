#ifndef __SimulationDataPBF_h__
#define __SimulationDataPBF_h__

#include "SPlisHSPlasH/Common.h"
#include <vector>
#include "SPlisHSPlasH/FluidModel.h"

namespace SPH 
{	
	/** \brief Simulation data which is required by the method Position-Based Fluids introduced
	* by Macklin and Mueller \cite Macklin:2013:PBF, \cite BMOTM2014, \cite BMM2015.
	*/
	class SimulationDataPBF
	{
		public:
			SimulationDataPBF();
			virtual ~SimulationDataPBF();

		protected:	
			std::vector<std::vector<Real>> m_lambda;		
			std::vector<std::vector<Vector3r>> m_deltaX;
			std::vector<std::vector<Vector3r>> m_oldX;
			std::vector<std::vector<Vector3r>> m_lastX;
			std::vector<std::vector<Real>> m_control_lambda;
			
			std::vector<std::vector<Vector3r>> c_deltaSpring;
			std::vector<std::vector<Vector3r>> c_deltaVelocity;
			std::vector<std::vector<Vector3r>> c_deltaDensity;
		public:
			/** Initialize the arrays containing the particle data.
			*/
			virtual void init();
			virtual void initControl();

			/** Release the arrays containing the particle data.
			*/
			virtual void cleanup();
			virtual void cleanupControl();

			/** Reset the particle data.
			*/
			virtual void reset();
			virtual void resetControl();

			/** Important: First call m_model->performNeighborhoodSearchSort()
			* to call the z_sort of the neighborhood search.
			*/
			void performNeighborhoodSearchSort();

			void controlNeighborhoodSearchSort();


			void emittedParticles(FluidModel *model, const unsigned int startIndex);

			FORCE_INLINE const Real& getLambda(const unsigned int fluidIndex, const unsigned int i) const
			{
				return m_lambda[fluidIndex][i];
			}

			FORCE_INLINE Real& getLambda(const unsigned int fluidIndex, const unsigned int i)
			{
				return m_lambda[fluidIndex][i];
			}

			FORCE_INLINE void setLambda(const unsigned int fluidIndex, const unsigned int i, const Real &val)
			{
				m_lambda[fluidIndex][i] = val;
			}

			FORCE_INLINE const Real& getControlLambda(const unsigned int controlIndex, const unsigned int i) const
			{
				return m_control_lambda[controlIndex][i];
			}

			FORCE_INLINE Real& getControlLambda(const unsigned int controlIndex, const unsigned int i)
			{
				return m_control_lambda[controlIndex][i];
			}

			FORCE_INLINE void setControlLambda(const unsigned int controlIndex, const unsigned int i, const Real& val)
			{
				m_control_lambda[controlIndex][i] = val;
			}

			FORCE_INLINE Vector3r& getDeltaX(const unsigned int fluidIndex, const unsigned int i)
			{
				return m_deltaX[fluidIndex][i];
			}

			FORCE_INLINE const Vector3r& getDeltaX(const unsigned int fluidIndex, const unsigned int i) const
			{
				return m_deltaX[fluidIndex][i];
			}

			FORCE_INLINE void setDeltaX(const unsigned int fluidIndex, const unsigned int i, const Vector3r& val)
			{
				m_deltaX[fluidIndex][i] = val;
			}
			//Spring
			FORCE_INLINE Vector3r& getDeltaSp(const unsigned int fluidIndex, const unsigned int i)
			{
				return c_deltaSpring[fluidIndex][i];
			}

			FORCE_INLINE const Vector3r& getDeltaSp(const unsigned int fluidIndex, const unsigned int i) const
			{
				return c_deltaSpring[fluidIndex][i];
			}

			FORCE_INLINE void setDeltaSp(const unsigned int fluidIndex, const unsigned int i, const Vector3r& val)
			{
				c_deltaSpring[fluidIndex][i] = val;
			}
			//Density
			FORCE_INLINE Vector3r& getDeltaDe(const unsigned int fluidIndex, const unsigned int i)
			{
				return c_deltaDensity[fluidIndex][i];
			}

			FORCE_INLINE const Vector3r& getDeltaDe(const unsigned int fluidIndex, const unsigned int i) const
			{
				return c_deltaDensity[fluidIndex][i];
			}

			FORCE_INLINE void setDeltaDe(const unsigned int fluidIndex, const unsigned int i, const Vector3r& val)
			{
				c_deltaDensity[fluidIndex][i] = val;
			}
			//Velocity
			FORCE_INLINE Vector3r& getDeltaVe(const unsigned int fluidIndex, const unsigned int i)
			{
				return c_deltaVelocity[fluidIndex][i];
			}

			FORCE_INLINE const Vector3r& getDeltaVe(const unsigned int fluidIndex, const unsigned int i) const
			{
				return c_deltaVelocity[fluidIndex][i];
			}

			FORCE_INLINE void setDeltaVe(const unsigned int fluidIndex, const unsigned int i, const Vector3r& val)
			{
				c_deltaVelocity[fluidIndex][i] = val;
			}
			//Lp
			FORCE_INLINE Vector3r &getLastPosition(const unsigned int fluidIndex, const unsigned int i)
			{
				return m_lastX[fluidIndex][i];
			}

			FORCE_INLINE const Vector3r &getLastPosition(const unsigned int fluidIndex, const unsigned int i) const
			{
				return m_lastX[fluidIndex][i];
			}

			FORCE_INLINE void setLastPosition(const unsigned int fluidIndex, const unsigned int i, const Vector3r &pos)
			{
				m_lastX[fluidIndex][i] = pos;
			}

			FORCE_INLINE Vector3r &getOldPosition(const unsigned int fluidIndex, const unsigned int i)
			{
				return m_oldX[fluidIndex][i];
			}

			FORCE_INLINE const Vector3r &getOldPosition(const unsigned int fluidIndex, const unsigned int i) const
			{
				return m_oldX[fluidIndex][i];
			}

			FORCE_INLINE void setOldPosition(const unsigned int fluidIndex, const unsigned int i, const Vector3r &pos)
			{
				m_oldX[fluidIndex][i] = pos;
			}

	};
}

#endif
#ifndef __VorticityRefinement_Liu2020_h__
#define __VorticityRefinement_Liu2020_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "VorticityBase.h"

namespace SPH
{
	class VorticityRefinement_Liu2020 : public VorticityBase
	{
	protected:
		std::vector<Vector3r> omega;
		std::vector<Vector3r> dOmega;
		Real m_alpha_adjustment;


		virtual void initParameters();

	public:
		static int ALPHA_ADJUSTMENT;

		VorticityRefinement_Liu2020(FluidModel* model);
		virtual ~VorticityRefinement_Liu2020(void);
		virtual void step();
		virtual void reset();

		virtual void performNeighborhoodSearchSort();
		//omega
		FORCE_INLINE const Vector3r& getOmega(const unsigned int i) const
		{
			return omega[i];
		}

		FORCE_INLINE Vector3r& getOmega(const unsigned int i)
		{
			return omega[i];
		}

		FORCE_INLINE void setOmega(const unsigned int i, const Vector3r& val)
		{
			omega[i] = val;
		}
		//dOmega
		FORCE_INLINE const Vector3r& getdOmega(const unsigned int i) const
		{
			return dOmega[i];
		}

		FORCE_INLINE Vector3r& getdOmega(const unsigned int i)
		{
			return dOmega[i];
		}

		FORCE_INLINE void setdOmega(const unsigned int i, const Vector3r& val)
		{
			dOmega[i] = val;
		}
	};
}

#endif
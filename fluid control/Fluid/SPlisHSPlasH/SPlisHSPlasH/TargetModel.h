#ifndef __TargetModel_
#define __TargetModel_
#include "Common.h"
#include <vector>

#include "SPHKernels.h"
#include "ParameterObject.h"
namespace SPH
{
	class TimeStep;
	class TargetModel :public GenParam::ParameterObject
	{
	public:

		TargetModel();
		virtual ~TargetModel();

		std::string getId() const { return m_id; }
	protected:
		std::string m_id;
		unsigned int t_pointSetIndex;

		unsigned int m_numActiveParticles;
		unsigned int m_numActiveParticles0;

		virtual void resizeTargetParticles(const unsigned int newSize);
		virtual void releaseTargetParticles();

		std::vector<Vector3r> m_x0;//getPosition0
		std::vector<Vector3r> m_x;//getPosition

	public:
		unsigned int numberOfParticles() const { return static_cast<unsigned int>(m_x.size()); }

		virtual void reset();

		void performNeighborhoodSearchSort();

		void initModel(const std::string& id, const unsigned int numControlBoundaryParticles, Vector3r* controlBoundaryParticles);

		unsigned int getPiontSetIndex() const { return t_pointSetIndex; }

		void setNumActiveParticles(const unsigned int num);
		unsigned int numActiveParticles() const;
		//x0
		FORCE_INLINE Vector3r& getPosition0(const unsigned int i)
		{
			return m_x0[i];
		}

		FORCE_INLINE const Vector3r& getPosition0(const unsigned int i) const
		{
			return m_x0[i];
		}

		FORCE_INLINE void setPosition0(const unsigned int i, const Vector3r& pos)
		{
			m_x0[i] = pos;
		}
		//m_x
		FORCE_INLINE Vector3r& getPosition(const unsigned int i)
		{
			return m_x[i];
		}

		FORCE_INLINE const Vector3r& getPosition(const unsigned int i) const
		{
			return m_x[i];
		}

		FORCE_INLINE void setPosition(const unsigned int i, const Vector3r& pos)
		{
			m_x[i] = pos;
		}
	};
}
#endif // __TargetModel_

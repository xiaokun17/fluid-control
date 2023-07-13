#ifndef __ControlBoundaryModel_h__
#define __ControlBoundaryModel_h__
#include "Common.h"
#include <vector>

#include "RigidBodyObject.h"
#include "SPHKernels.h"
#include "ParameterObject.h"
namespace SPH
{
	class TimeStep;
	class ControlBoundaryModel:public GenParam::ParameterObject
	{
		public:
			static int NUM_PARTICLES;
			static int NUM_REUSED_PARTICLES;
			static int DENSITY0;

			ControlBoundaryModel();
			virtual ~ControlBoundaryModel();

			std::string getId() const { return m_id; }
		protected:
			std::string m_id;

			Real m_density0;
			unsigned int c_pointSetIndex;
			unsigned int t_pointSetIndex;

			unsigned int m_numActiveParticles;
			unsigned int m_numActiveParticles0;

			virtual void resizeControlParticles(const unsigned int newSize);
			virtual void releaseControlParticles();

			RigidBodyObject* m_rigidBody;
			std::vector<Vector3r> m_x0;//getPosition0
			std::vector<Vector3r> m_x;//getPosition
			std::vector<Vector3r> __x;//_x
			std::vector<Real> m_density;//getDensity
			std::vector<Real> m_V;//getVolume
			std::vector<Real> m_C;//colorValue
			std::vector<int> t_x;//target index
			std::vector<int> edge;//是否为浅层
			std::vector<int> face;//第几层索引
	public:
		unsigned int numberOfParticles() const { return static_cast<unsigned int>(m_V.size()); }
		
		void computeControlBoundaryVolume();

		virtual void reset();

		void performNeighborhoodSearchSort();

		void initModel(const std::string& id, const unsigned int numControlBoundaryParticles, Vector3r* controlBoundaryParticles);
		void initModel(RigidBodyObject* rbo,const unsigned int id, const unsigned int numControlBoundaryParticles, Vector3r* controlBoundaryParticles);
		RigidBodyObject* getRigidBodyObject() { return m_rigidBody; }

		FORCE_INLINE Real getDensity0() const { return m_density0; }
		void setDensity0(const Real v) { m_density0 = v; }

		unsigned int getPiontSetIndex() const { return c_pointSetIndex; }
		unsigned int getTargetPointSetIndex() const { return t_pointSetIndex; }

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
		//__X
				//m_x
		FORCE_INLINE Vector3r& getPosition_(const unsigned int i)
		{
			return __x[i];
		}

		FORCE_INLINE const Vector3r& getPosition_(const unsigned int i) const
		{
			return __x[i];
		}

		FORCE_INLINE void setPosition_(const unsigned int i, const Vector3r& pos)
		{
			__x[i] = pos;
		}
		//density
		FORCE_INLINE const Real& getDensity(const unsigned int i) const
		{
			return m_density[i];
		}

		FORCE_INLINE Real& getDensity(const unsigned int i)
		{
			return m_density[i];
		}

		FORCE_INLINE void setDensity(const unsigned int i, const Real& val)
		{
			m_density[i] = val;
		}
		//m_V
		FORCE_INLINE const Real& getVolume(const unsigned int i) const
		{
			return m_V[i];
		}

		FORCE_INLINE Real& getVolume(const unsigned int i)
		{
			return m_V[i];
		}

		FORCE_INLINE void setVolume(const unsigned int i, const Real& val)
		{
			m_V[i] = val;
		}
		//m_C
		FORCE_INLINE const Real& getC(const unsigned int i) const
		{
			return m_C[i];
		}

		FORCE_INLINE Real& getC(const unsigned int i)
		{
			return m_C[i];
		}

		FORCE_INLINE void setC(const unsigned int i, const Real& val)
		{
			m_C[i] = val;
		}
		//edge
		FORCE_INLINE const int& getEdge(const unsigned int i) const
		{
			return edge[i];
		}

		FORCE_INLINE int& getEdge(const unsigned int i)
		{
			return edge[i];
		}

		FORCE_INLINE void setEdge(const unsigned int i, const int& val)
		{
			edge[i] = val;
		}
		//facce
		FORCE_INLINE const int& getFace(const unsigned int i) const
		{
			return face[i];
		}

		FORCE_INLINE int& getFace(const unsigned int i)
		{
			return face[i];
		}

		FORCE_INLINE void setFace(const unsigned int i, const int& val)
		{
			face[i] = val;
		}
		//target position
		FORCE_INLINE int& getTargetPos(const unsigned int i)
		{
			return t_x[i];
		}

		FORCE_INLINE const int& getTargetPos(const unsigned int i) const
		{
			return t_x[i];
		}
	};
}
#endif // !__ControlModel_Boundary_h__
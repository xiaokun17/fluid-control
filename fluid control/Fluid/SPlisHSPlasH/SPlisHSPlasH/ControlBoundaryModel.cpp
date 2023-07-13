#include "ControlBoundaryModel.h"
#include "SPHKernels.h"
#include <iostream>
#include "TimeManager.h"
#include "TimeStep.h"
#include "Utilities/Logger.h"
#include "NeighborhoodSearch.h"
#include "Simulation.h"

using namespace SPH;

ControlBoundaryModel::ControlBoundaryModel() :
	m_x0(),
	m_x(),
	__x(),
	m_V(),
	m_C(),
	t_x(),
	edge(),
	face(),
	m_density()
{
	m_density0 = 200;
	c_pointSetIndex = 0;
	t_pointSetIndex = 0;
}


ControlBoundaryModel::~ControlBoundaryModel(void)
{
	releaseControlParticles();
	delete m_rigidBody;
}


void ControlBoundaryModel::reset()
{
	setNumActiveParticles(m_numActiveParticles0);
	const unsigned int nPoints = numActiveParticles();
	
	for (unsigned int i = 0; i < nPoints; i++)
	{
		const Vector3r& x0 = getPosition0(i);
		getPosition(i) = x0;
		getPosition_(i).setZero();
		getTargetPos(i) = -1;
		m_density[i] = 0.0;
		edge[i] = 0;
		face[i] = -1;
	}
	NeighborhoodSearch* cneighborhoodSearch = Simulation::getCurrent()->getCNeighborhoodSearch();
	if (cneighborhoodSearch->point_set(c_pointSetIndex).n_points() != nPoints)
		cneighborhoodSearch->resize_point_set(c_pointSetIndex, &getPosition(0)[0], nPoints);
	//if (Simulation::getCurrent()->ifTarget())
	//{
	//	NeighborhoodSearch* tneighborhoodSearch = Simulation::getCurrent()->getTNeighborhoodSearch();
	//	if (tneighborhoodSearch->point_set(t_pointSetIndex).n_points() != nPoints)
	//		tneighborhoodSearch->resize_point_set(t_pointSetIndex, &getPosition(0)[0], nPoints);
	//}

}

void ControlBoundaryModel::computeControlBoundaryVolume()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundarys = sim->numberOfBoundaryModels();
	NeighborhoodSearch* neighborhoodSearch = Simulation::getCurrent()->getNeighborhoodSearch();

	const Real particleRadius = sim->getParticleRadius();
	const Real diam = static_cast<Real>(2.0)* particleRadius;
	
	const unsigned int numControlBoundaryParticles = numberOfParticles();
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numControlBoundaryParticles; i++)
		{
			m_V[i] = static_cast<Real>(0.8)* diam* diam* diam;

			//Real delta = sim->W_zero();
			//for (unsigned int pid = nFluids; pid < sim->numberOfPointSets(); pid++)
			//{
			//	
			//	ControlBoundaryModel* cbm_neighbor = sim->getControlBoundaryModelFromPointSet(pid);
			//	for (unsigned int j = 0; j < neighborhoodSearch->point_set(c_pointSetIndex).n_neighbors(pid, i); j++)
			//	{
			//		const unsigned int neighborIndex = neighborhoodSearch->point_set(c_pointSetIndex).neighbor(pid, i, j);
			//		delta += sim->W(getPosition(i) - cbm_neighbor->getPosition(neighborIndex));
			//	}
			//}
			//const Real volume = static_cast<Real>(1.0) / delta;
			//m_V[i] = volume;
		}
	}
}

void ControlBoundaryModel::initModel(const std::string& id, const unsigned int numControlBoundaryParticles, Vector3r* controlBoundaryParticles)
{
	m_id = id;
	releaseControlParticles();
	resizeControlParticles(numControlBoundaryParticles);

#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)
		for (int i = 0; i < (int)numControlBoundaryParticles; i++)
		{
			getPosition0(i) = controlBoundaryParticles[i];
			getPosition(i) = controlBoundaryParticles[i];
			getTargetPos(i) = -1;
			m_density[i] = 0;
			edge[i] = 0;
			face[i] = -1;
			//m_C[i] = 0;
		}
	}
	NeighborhoodSearch* cneighborhoodSearch = Simulation::getCurrent()->getCNeighborhoodSearch();
	c_pointSetIndex = cneighborhoodSearch->add_point_set(&getPosition(0)[0], numControlBoundaryParticles, false, true, true, this);
	m_numActiveParticles0 = numControlBoundaryParticles;
	m_numActiveParticles = m_numActiveParticles0;
	//if (Simulation::getCurrent()->ifTarget())
	//{
	//	NeighborhoodSearch* tneighborhoodSearch = Simulation::getCurrent()->getTNeighborhoodSearch();
	//	t_pointSetIndex = tneighborhoodSearch->add_point_set(&getPosition(0)[0], numControlBoundaryParticles, false, true, true, this);
	//}
}

void ControlBoundaryModel::initModel(RigidBodyObject* rbo,const unsigned int id, const unsigned int numControlBoundaryParticles, Vector3r* controlBoundaryParticles)
{
	releaseControlParticles();
	resizeControlParticles(numControlBoundaryParticles);

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)
		for (int i = 0; i < (int)numControlBoundaryParticles; i++)
		{
			getPosition0(i) = controlBoundaryParticles[i];
			getPosition(i) = controlBoundaryParticles[i];
			getPosition_(i).setZero();
			getTargetPos(i) = -1;
			m_density[i] = 0;
		}
	}
	m_rigidBody = rbo;
	NeighborhoodSearch* cneighborhoodSearch = Simulation::getCurrent()->getCNeighborhoodSearch();
	c_pointSetIndex = cneighborhoodSearch->add_point_set(&getPosition(0)[0], numControlBoundaryParticles, false , true, true, this);
	m_numActiveParticles0 = numControlBoundaryParticles;
	m_numActiveParticles = m_numActiveParticles0;
	//if (Simulation::getCurrent()->ifTarget())
	//{
	//	NeighborhoodSearch* tneighborhoodSearch = Simulation::getCurrent()->getTNeighborhoodSearch();
	//	t_pointSetIndex = tneighborhoodSearch->add_point_set(&getPosition(0)[0], numControlBoundaryParticles, false, true, true, this);
	//}
}

void ControlBoundaryModel::performNeighborhoodSearchSort()
{
	const unsigned int numPart = numActiveParticles();
	if (numPart == 0)
		return;

	NeighborhoodSearch* cneighborhoodSearch = Simulation::getCurrent()->getCNeighborhoodSearch();

	auto const& d = cneighborhoodSearch->point_set(c_pointSetIndex);
	d.sort_field(&m_x[0]);
	d.sort_field(&__x[0]);
	d.sort_field(&m_V[0]);
	d.sort_field(&m_C[0]);
	d.sort_field(&edge[0]);
	d.sort_field(&face[0]);
	d.sort_field(&m_density[0]);
	d.sort_field(&t_x[0]);
}

void ControlBoundaryModel::setNumActiveParticles(const unsigned int num)
{
	m_numActiveParticles = num;
}

unsigned int ControlBoundaryModel::numActiveParticles() const
{
	return m_numActiveParticles;
}

void ControlBoundaryModel::releaseControlParticles()
{
	m_x0.clear();
	m_x.clear();
	__x.clear();
	m_V.clear();
	m_C.clear();
	edge.clear();
	face.clear();
	m_density.clear();
	t_x.clear();
}

void ControlBoundaryModel::resizeControlParticles(const unsigned int newSize)
{
	m_x0.resize(newSize);
	m_x.resize(newSize);
	__x.resize(newSize);
	m_V.resize(newSize);
	m_C.resize(newSize);
	edge.resize(newSize);
	face.resize(newSize);
	m_density.resize(newSize);
	t_x.resize(newSize);
}

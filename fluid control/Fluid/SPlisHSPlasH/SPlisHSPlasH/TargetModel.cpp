#include "TargetModel.h"
#include "SPHKernels.h"
#include <iostream>
#include "TimeManager.h"
#include "TimeStep.h"
#include "Utilities/Logger.h"
#include "NeighborhoodSearch.h"
#include "Simulation.h"

using namespace SPH;

TargetModel::TargetModel() :
	m_x0(),
	m_x()
{
	t_pointSetIndex = 0;
}

TargetModel::~TargetModel(void)
{
	releaseTargetParticles();
}

void TargetModel::reset()
{
	setNumActiveParticles(m_numActiveParticles0);
	const unsigned int nPoints = numActiveParticles();

	for (unsigned int i = 0; i < nPoints; i++)
	{
		const Vector3r& x0 = getPosition0(i);
		getPosition(i) = x0;
	}
	NeighborhoodSearch* tneighorhoodSearch = Simulation::getCurrent()->getTNeighborhoodSearch();
	if (tneighorhoodSearch->point_set(t_pointSetIndex).n_points() != nPoints)
		tneighorhoodSearch->resize_point_set(t_pointSetIndex, &getPosition(0)[0], nPoints);
}

void TargetModel::initModel(const std::string& id, const unsigned int numTargetParticles, Vector3r* TargetParticles)
{
	m_id = id;
	releaseTargetParticles();
	resizeTargetParticles(numTargetParticles);

#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)
		for (int i = 0; i < (int)numTargetParticles; i++)
		{
			getPosition0(i) = TargetParticles[i];
			getPosition(i) = TargetParticles[i];
		}
	}
	NeighborhoodSearch* tneighborhoodSearch = Simulation::getCurrent()->getTNeighborhoodSearch();
	t_pointSetIndex = tneighborhoodSearch->add_point_set(&getPosition(0)[0], numTargetParticles, false, true, false, this);
	m_numActiveParticles0 = numTargetParticles;
	m_numActiveParticles = 0;
}

void TargetModel::performNeighborhoodSearchSort()
{
	const unsigned int numPart = numActiveParticles();
	if (numPart == 0)
		return;

	NeighborhoodSearch* tneighborhoodSearch = Simulation::getCurrent()->getTNeighborhoodSearch();

	auto const& d = tneighborhoodSearch->point_set(t_pointSetIndex);
	d.sort_field(&m_x[0]);
}
void TargetModel::setNumActiveParticles(const unsigned int num)
{
	m_numActiveParticles = num;
}

unsigned int TargetModel::numActiveParticles() const
{
	return m_numActiveParticles;
}

void TargetModel::releaseTargetParticles()
{
	m_x0.clear();
	m_x.clear();
}

void TargetModel::resizeTargetParticles(const unsigned int newSize)
{
	m_x0.resize(newSize);
	m_x.resize(newSize);
}

#include "SimulationDataPBF.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SPlisHSPlasH/Simulation.h"

using namespace SPH;

SimulationDataPBF::SimulationDataPBF() :
	m_deltaX(),
	m_control_lambda(),
	m_lambda(),
	m_lastX(),
	m_oldX(),
	c_deltaDensity(),
	c_deltaSpring(),
	c_deltaVelocity()
{
}

SimulationDataPBF::~SimulationDataPBF(void)
{
	cleanup();
}

void SimulationDataPBF::initControl()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nControlModels = sim->numberOfControlBoundaryModels();
	m_control_lambda.resize(nControlModels);

	for (unsigned int i = 0; i < nControlModels; i++)
	{
		ControlBoundaryModel *cbm = sim->getControlBoundaryModel(i);
		m_control_lambda[i].resize(cbm->numberOfParticles(), 0.0);
	}
}
void SimulationDataPBF::init()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	const unsigned int nControls = sim->numberOfControlBoundaryModels();

	m_lambda.resize(nModels);
	m_deltaX.resize(nModels);
	m_oldX.resize(nModels);
	m_lastX.resize(nModels);

	c_deltaDensity.resize(nModels);
	c_deltaSpring.resize(nModels);
	c_deltaVelocity.resize(nModels);
	for (unsigned int i = 0; i < nModels; i++)
	{
		FluidModel *fm = sim->getFluidModel(i);
		m_lambda[i].resize(fm->numParticles(), 0.0);
		m_deltaX[i].resize(fm->numParticles(), Vector3r::Zero());
		m_oldX[i].resize(fm->numParticles(), Vector3r::Zero());
		m_lastX[i].resize(fm->numParticles(), Vector3r::Zero());
		c_deltaDensity[i].resize(fm->numberOfParticles(), Vector3r::Zero());
		c_deltaSpring[i].resize(fm->numberOfParticles(), Vector3r::Zero());
		c_deltaVelocity[i].resize(fm->numberOfParticles(), Vector3r::Zero());
	}
	initControl();
	reset();
}
void SimulationDataPBF::cleanupControl()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nControlModels = sim->numberOfBoundaryModels();
	for (unsigned int i = 0; i < nControlModels; i++)
	{
		m_control_lambda[i].clear();

	}
	m_control_lambda.clear();

}
void SimulationDataPBF::cleanup()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int i = 0; i < nModels; i++)
	{
		m_lambda[i].clear();
		m_deltaX[i].clear();
		m_oldX[i].clear();
		m_lastX[i].clear();
		c_deltaDensity[i].clear();
		c_deltaSpring[i].clear();
		c_deltaVelocity[i].clear();
	}
	m_lambda.clear();
	m_deltaX.clear();
	m_oldX.clear();
	m_lastX.clear();
	c_deltaDensity.clear();
	c_deltaSpring.clear();
	c_deltaVelocity.clear();
	cleanupControl();
}

void SimulationDataPBF::resetControl()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nControlModels = sim->numberOfControlBoundaryModels();
	for (unsigned int i = 0; i < nControlModels; i++)
	{
		ControlBoundaryModel* cbm = sim->getControlBoundaryModel(i);
		for (unsigned int j = 0; j < cbm->numActiveParticles(); j++)
		{
			m_control_lambda[i][j] = 0.0;
		}
	}
}

void SimulationDataPBF::reset()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int i = 0; i < nModels; i++)
	{
		FluidModel *fm = sim->getFluidModel(i);
		for (unsigned int j = 0; j < fm->numActiveParticles(); j++)
		{
			m_deltaX[i][j].setZero();
			c_deltaDensity[i][j].setZero();
			c_deltaSpring[i][j].setZero();
			c_deltaVelocity[i][j].setZero();
			m_lambda[i][j] = 0.0;
			getLastPosition(i, j) = fm->getPosition(j);
			getOldPosition(i, j) = fm->getPosition(j);
		}
	}
	resetControl();
}

void SimulationDataPBF::performNeighborhoodSearchSort()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int i = 0; i < nModels; i++)
	{
		FluidModel *fm = sim->getFluidModel(i);
		const unsigned int numPart = fm->numActiveParticles();
		if (numPart != 0)
		{
			auto const& d = sim->getNeighborhoodSearch()->point_set(fm->getPointSetIndex());
			d.sort_field(&m_lambda[i][0]);
			d.sort_field(&m_deltaX[i][0]);
			d.sort_field(&m_oldX[i][0]);
			d.sort_field(&m_lastX[i][0]);
			d.sort_field(&c_deltaDensity[i][0]);
			d.sort_field(&c_deltaSpring[i][0]);
			d.sort_field(&c_deltaVelocity[i][0]);

			auto const& cd = sim->getCNeighborhoodSearch()->point_set(fm->getCponitSetIndex());
			cd.sort_field(&m_lambda[i][0]);
			cd.sort_field(&m_deltaX[i][0]);
			cd.sort_field(&m_oldX[i][0]);
			cd.sort_field(&m_lastX[i][0]);
			cd.sort_field(&c_deltaDensity[i][0]);
			cd.sort_field(&c_deltaSpring[i][0]);
			cd.sort_field(&c_deltaVelocity[i][0]);
		}
	}
	controlNeighborhoodSearchSort();
}

void SimulationDataPBF::controlNeighborhoodSearchSort()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nControls = sim->numberOfControlBoundaryModels();
	for (unsigned int i = 0; i < nControls; i++)
	{
		ControlBoundaryModel* cbm = sim->getControlBoundaryModel(i);
		const unsigned int numPart = cbm->numActiveParticles();
		if (numPart != 0)
		{
			auto const& d = sim->getCNeighborhoodSearch()->point_set(cbm->getPiontSetIndex());
			d.sort_field(&m_control_lambda[i][0]);
			//to do
		}
	}
}
void SimulationDataPBF::emittedParticles(FluidModel *model, const unsigned int startIndex)
{
	// initialize lastX values for new particles
	Simulation *sim = Simulation::getCurrent();
	const unsigned int fluidModelIndex = model->getPointSetIndex();
	for (unsigned int j = startIndex; j < model->numActiveParticles(); j++)
	{
		m_lastX[fluidModelIndex][j] = model->getPosition(j);
	}
}


#include "VorticityRefinement_Liu2020.h"
#include <iostream>
#include "../TimeManager.h"
#include "../Simulation.h"

using namespace SPH;
using namespace GenParam;

int VorticityRefinement_Liu2020::ALPHA_ADJUSTMENT = -1;

VorticityRefinement_Liu2020::VorticityRefinement_Liu2020(FluidModel *model) :
	VorticityBase(model)
{
	omega.resize(model->numParticles(), Vector3r::Zero());
	dOmega.resize(model->numParticles(), Vector3r::Zero());

	m_alpha_adjustment = 0.6;

	model->addField({ "Omega Liu", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &omega[i][0]; } });
}

VorticityRefinement_Liu2020::~VorticityRefinement_Liu2020(void)
{
	m_model->removeFieldByName("Omega Liu");

	omega.clear();
	dOmega.clear();
}

void VorticityRefinement_Liu2020::initParameters()
{
	VorticityBase::initParameters();

	ALPHA_ADJUSTMENT = createNumericParameter("vorticityAlpha", "adjustment parameter alpha", &m_alpha_adjustment);
	setGroup(ALPHA_ADJUSTMENT, "Vorticity");
	setDescription(ALPHA_ADJUSTMENT, "VR adjustment parameter");
	RealParameter* rparam = static_cast<RealParameter*>(getParameter(ALPHA_ADJUSTMENT));
	rparam->setMinValue(0.0);
}

void VorticityRefinement_Liu2020::step()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int numParticles = m_model->numActiveParticles();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	const unsigned int nFluids = sim->numberOfFluidModels();
	FluidModel* model = m_model;
	const Real density0 = model->getDensity0();

	const Real dt = TimeManager::getCurrent()->getTimeStepSize();
	const Real invDt = static_cast<Real>(1.0) / dt;

	const Real alphaVR = m_alpha_adjustment;
	const Real nu = m_vorticityCoeff;

	const Real h = sim->getSupportRadius();
	const Real h2 = h * h;
	//d ?? 
	Real d = 10.0;
	if (sim->is2DSimulation())
		d = 6.0;

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r& xi = m_model->getPosition(i);
			const Vector3r& vi = m_model->getVelocity(i);
			Vector3r& omegai = omega[i];
			Vector3r& dOmegai = dOmega[i];
			dOmegai.setZero();
			const Real density_i = m_model->getDensity(i);

			Vector3r tmp = Vector3r::Zero();
			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors_in_same_phase(
				const Vector3r &vj = m_model->getVelocity(neighborIndex);
				const Vector3r xij = xi - xj;
				const Vector3r gradW = sim->gradW(xij);
				const Real density_j = m_model->getDensity(neighborIndex);
				const Real mass_j = m_model->getMass(neighborIndex);
				tmp += ((vi - vj)*(1 / density_j*mass_j)).cross(gradW);
			);
			dOmegai = omegai - tmp;
			omegai = tmp;
		}
	}

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r& xi = m_model->getPosition(i);
			Vector3r vi = m_model->getVelocity(i);
			forall_fluid_neighbors_in_same_phase(
				const Vector3r &dOmegaj = dOmega[neighborIndex];
				const Vector3r xij = xi - xj;
				const Real rc = sqrt(4 * 1.25643 * dt * nu) * alphaVR;
				Vector3r tmp = (dOmegaj / 2).cross(xij) * pow(rc / xij.norm(), 2);
				vi += tmp;
			);
		}
	}
}

void VorticityRefinement_Liu2020::reset()
{
	for (unsigned int i = 0; i < m_model->numActiveParticles(); i++)
	{
		omega[i].setZero();
		dOmega[i].setZero();// to do ...
	}
}

void SPH::VorticityRefinement_Liu2020::performNeighborhoodSearchSort()
{
	const unsigned int numPart = m_model->numActiveParticles();
	if (numPart == 0)
		return;

	Simulation* sim = Simulation::getCurrent();
	auto const& d = sim->getNeighborhoodSearch()->point_set(m_model->getPointSetIndex());
	d.sort_field(&omega[0]);
}
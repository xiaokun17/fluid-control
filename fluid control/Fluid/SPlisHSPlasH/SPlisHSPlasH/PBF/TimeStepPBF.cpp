#include "TimeStepPBF.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "TimeIntegration.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SimulationDataPBF.h"
#include <iostream>
#include "Utilities/Timing.h"
#include "SPlisHSPlasH/Simulation.h"
#include "Utilities/Counting.h"

using namespace SPH;
using namespace std;
using namespace GenParam;

int TimeStepPBF::VELOCITY_UPDATE_METHOD = -1;
int TimeStepPBF::ENUM_PBF_FIRST_ORDER = -1;
int TimeStepPBF::ENUM_PBF_SECOND_ORDER = -1;
int TimeStepPBF::DENSITY_ALPHA = -1;
int TimeStepPBF::SPRING_BETA = -1;
int TimeStepPBF::VELOCITY_GAMMA = -1;
int TimeStepPBF::ENABLE_TARGET = -1;



TimeStepPBF::TimeStepPBF() :
	TimeStep()
{
	m_simulationData.init();
	m_counter = 0;
	m_velocityUpdateMethod = 0;

	alpha = 0.0;//Density 2021/3 0.3 0.5可行
	beta = 0.0;//Spring
	gamma = 0.0;//Velocity
	enableTarget = false;

	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		model->addField({ "lambda", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getLambda(fluidModelIndex, i); } });
		model->addField({ "deltaX", FieldType::Vector3, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getDeltaX(fluidModelIndex, i)[0]; } });
	}


}

TimeStepPBF::~TimeStepPBF(void)
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		model->removeFieldByName("lambda");
		model->removeFieldByName("deltaX");
	}
}

void TimeStepPBF::initParameters()
{
	TimeStep::initParameters();

	VELOCITY_UPDATE_METHOD = createEnumParameter("velocityUpdateMethod", "Velocity update method", &m_velocityUpdateMethod);
	setGroup(VELOCITY_UPDATE_METHOD, "PBF");
	setDescription(VELOCITY_UPDATE_METHOD, "Method for the velocity integration.");
	EnumParameter *enumParam = static_cast<EnumParameter*>(getParameter(VELOCITY_UPDATE_METHOD));
	enumParam->addEnumValue("First Order Update", ENUM_PBF_FIRST_ORDER);
	enumParam->addEnumValue("Second Order Update", ENUM_PBF_SECOND_ORDER);


	if (Simulation::getCurrent()->ifControl())
	{
		DENSITY_ALPHA = createNumericParameter("densityConstraint", "Density Constraint", &alpha);
		setGroup(DENSITY_ALPHA, "PBF");
		setDescription(DENSITY_ALPHA, "Density constraint coefficient alpha.");
		static_cast<RealParameter*>(getParameter(DENSITY_ALPHA))->setMinValue(1e-6);

		SPRING_BETA = createNumericParameter("springConstraint", "Spring Constraint", &beta);
		setGroup(SPRING_BETA, "PBF");
		setDescription(SPRING_BETA, "Spring constraint coefficient beta.");
		static_cast<RealParameter*>(getParameter(SPRING_BETA))->setMinValue(1e-6);

		VELOCITY_GAMMA = createNumericParameter("velocityConstraint", "VelocityConstraint", &gamma);
		setGroup(VELOCITY_GAMMA, "PBF");
		setDescription(VELOCITY_GAMMA, "Velocity constraint coefficient gamma.");
		static_cast<RealParameter*>(getParameter(VELOCITY_GAMMA))->setMinValue(1e-6);
		//if (Simulation::getCurrent()->ifTarget())
		//{
		//	ENABLE_TARGET = createBoolParameter("EnableTarget", "Enabel Target", &enableTarget);
		//	setGroup(ENABLE_TARGET, "PBF");
		//	setDescription(ENABLE_TARGET, "Enable/disable Target.");
		//}
	}


}

void TimeStepPBF::step()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	TimeManager *tm = TimeManager::getCurrent ();
	const Real h = tm->getTimeStepSize();

	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		clearAccelerations(fluidModelIndex);

		// Time integration
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < (int) model->numActiveParticles(); i++)
			{
				m_simulationData.getLastPosition(fluidModelIndex, i) = m_simulationData.getOldPosition(fluidModelIndex, i);
				m_simulationData.getOldPosition(fluidModelIndex, i) = model->getPosition(i);
				TimeIntegration::semiImplicitEuler(h, model->getMass(i), model->getPosition(i), model->getVelocity(i), model->getAcceleration(i));
			}
		}
	}
	
	// Perform neighborhood search
	
	performNeighborhoodSearch();
	
	if (sim->ifControl())
	{
		//if control particles move then
		//controlNeighborhoodSearchSort();
		controlPressureSolve();
		controlConstraints();
	}

	// Solve density constraint
	START_TIMING("pressureSolve");
	pressureSolve();
	STOP_TIMING_AVG;
	
	// Update velocities	
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		if (m_velocityUpdateMethod == ENUM_PBF_FIRST_ORDER)
		{
			#pragma omp parallel default(shared)
			{
				#pragma omp for schedule(static)  
				for (int i = 0; i < (int)model->numActiveParticles(); i++)
				{
					TimeIntegration::velocityUpdateFirstOrder(h, model->getMass(i), model->getPosition(i), m_simulationData.getOldPosition(fluidModelIndex, i), model->getVelocity(i));
					model->getAcceleration(i).setZero();
				}
			}
		}
		else
		{
			#pragma omp parallel default(shared)
			{
				#pragma omp for schedule(static)  
				for (int i = 0; i < (int)model->numActiveParticles(); i++)
				{
					TimeIntegration::velocityUpdateSecondOrder(h, model->getMass(i), model->getPosition(i), m_simulationData.getOldPosition(fluidModelIndex, i), m_simulationData.getLastPosition(fluidModelIndex, i), model->getVelocity(i));
					model->getAcceleration(i).setZero();
				}
			}
		}
	}

	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
		computeDensities(fluidModelIndex);
	sim->computeNonPressureForces();

	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < (int)model->numActiveParticles(); i++)
			{
				model->getVelocity(i) += h * model->getAcceleration(i);
			}
		}
	}

	sim->updateTimeStepSize();

	const Real stopEmit =  sim->getValue<Real>(sim->STOP_EMIT_AT);
	if ((stopEmit < 0.0) || (stopEmit > TimeManager::getCurrent()->getTime()))
		sim->emitParticles();
	else
		sim->stopEmit();
		

	// Compute new time	
	tm->setTime (tm->getTime () + h);
}


void TimeStepPBF::reset()
{
	TimeStep::reset();
	m_simulationData.reset();
	m_counter = 0;
	alpha = 0.0;
	beta = 0.0;
	gamma = 0.0;
}
/*
void TimeStepPBF::goToTarget()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nControls = sim->numberOfControlBoundaryModels();
	for (unsigned int controlIndex = 0; controlIndex < nControls; controlIndex++)
	{
		ControlBoundaryModel* model = sim->getControlBoundaryModel(controlIndex);
		const unsigned int numParticles = model->numActiveParticles();
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)
			for (int i = 0; i < numParticles; i++)
			{
				const Vector3r& xi = model->getPosition(i);
				int& ti = model->getTargetPos(i);
				for (int pid = nControls; pid < sim->numberOfTargetPointSets(); pid++)
				{
					ControlBoundaryModel* cbm_neighbor = sim->getControlBoundaryModelFromPointSet(pid);
					TargetModel* tm_neighbor = sim->getTargetModelFromPointSet(pid);
					//Real ri_dis = supportradius;
					for (unsigned int j = 0; j<sim->numberOfTargetNeighbors(controlIndex, pid, i); j++)
					{
						const unsigned int neighborIndex = sim->getTargetNeighbor(controlIndex, pid, i, j);
						const Vector3r& xj = tm_neighbor->getPosition(neighborIndex);
						const Real xj_dis = (xi - xj).norm();
					}

				}
			}
		}
	}

}
*/
void TimeStepPBF::controlPressureSolve()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nControls = sim->numberOfControlBoundaryModels();
	for (unsigned int i = 0; i < nControls; i++)
	{
		computeLambda(i);
	}
}
void TimeStepPBF::pressureSolve()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();

	Real avg_density_err = 0;
	m_iterations = 0;
	bool chk = false;
	while ((!chk || (m_iterations < 2)) && (m_iterations < m_maxIterations))
	{
		chk = true;
		for (unsigned int i = 0; i < nFluids; i++)
		{
			FluidModel *model = sim->getFluidModel(i);
			const Real density0 = model->getDensity0();

			avg_density_err = 0.0;
			
			pressureSolveIteration(i, avg_density_err);
			
			// Maximal allowed density fluctuation
			const Real eta = m_maxError * static_cast<Real>(0.01) * density0;  // maxError is given in percent
			chk = chk && (avg_density_err <= eta);
		}

		m_iterations++;
	}
	INCREASE_COUNTER("PBF - iterations", m_iterations);
}

void TimeStepPBF::computeLambda(const unsigned int controlModelIndex)
{
	Simulation* sim = Simulation::getCurrent();
	ControlBoundaryModel* model = sim->getControlBoundaryModel(controlModelIndex);
	const unsigned int numParticles = model->numActiveParticles();
	const unsigned int nControls = sim->numberOfControlBoundaryModels();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundary = sim->numberOfBoundaryModels();
	const Real supportradius = sim->getSupportRadius();
	const Real eps = 1.0e-6;
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)
		for (int i = 0; i < (int)numParticles; i++)
		{
			Real& density = model->getDensity(i);
			//density = 0.0;
			density = sim->CW_zero();
			const Real inv = 1.0 / density;
			const Vector3r& xi = model->getPosition(i);
			const Real& Ci = model->getC(i);
			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int pid = 0; pid < nFluids; pid++)
			{
				FluidModel* fm_neighbor = sim->getFluidModelFromPointSetControl(pid);
				for (unsigned int j = 0; j < sim->cnumberOfNeighbors(controlModelIndex+nFluids, pid, i); j++)
				{
					const unsigned int neighborIndex = sim->cgetNeighbor(controlModelIndex + nFluids, pid, i, j);
					const Vector3r &xj = fm_neighbor->getPosition(neighborIndex);
					//if((xi-xj).norm()<=supportradius)
						//density += fm_neighbor->getVolume(neighborIndex)* sim->CW(xi - xj);
						density += sim->CW(xi - xj);
						//density += fm_neighbor->getMass(i) * sim->CW(xi - xj);
				}
				
			}
			//for (unsigned int pid = nFluids; pid < sim->numberOfPointSets(); pid++)
			//{
			//	BoundaryModel* bm_neighbor = sim->getBoundaryModelFromPointSet(pid);
			//	for (unsigned int j = 0; j < sim->numberOfNeighbors(fluidModelIndex, pid, i); j++)
			//	{
			//		const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, pid, i, j);
			//		const Vector3r& xj = bm_neighbor->getPosition(neighborIndex);
			//		density += bm_neighbor->getVolume(neighborIndex) * sim->W(xi - xj);
			//	}
			//}
			//const Real C = std::max(density/model->getDensity0() - static_cast<Real>(1.0), static_cast<Real>(0.0));
			//const Real C = std::max(density*inv - static_cast<Real>(1.0), static_cast<Real>(0.0))*Ci;
			const Real C = std::max(density * inv - static_cast<Real>(1.0), static_cast<Real>(0.0));
			if (C != 0)
			{
				// Compute gradients dC/dx_j 
				Real sum_grad_C2 = 0.0;
				Vector3r gradC_i(0.0, 0.0, 0.0);
				//////////////////////////////////////////////////////////////////////////
				// Fluid
				//////////////////////////////////////////////////////////////////////////
				for (unsigned int pid = 0; pid < nFluids; pid++)
				{
					FluidModel* fm_neighbor = sim->getFluidModelFromPointSetControl(pid);
					for (unsigned int j = 0; j < sim->cnumberOfNeighbors(controlModelIndex + nFluids, pid, i); j++)
					{
						const unsigned int neighborIndex = sim->cgetNeighbor(controlModelIndex + nFluids, pid, i, j);
						const Vector3r& xj = fm_neighbor->getPosition(neighborIndex);
						//if ((xi - xj).norm() <= supportradius)
						//{
							const Vector3r gradC_j = sim->gradCW(xi - xj) * inv;
							//const Vector3r gradC_j = fm_neighbor->getVolume(neighborIndex)*sim->gradCW(xi - xj);
							sum_grad_C2 -= gradC_j.squaredNorm();
							gradC_i += gradC_j;
						//}

					}
				}

				//sum_grad_C2 += gradC_i.squaredNorm();
				// Compute lambda
				Real& lambda = m_simulationData.getControlLambda(controlModelIndex, i);
				lambda = -C / (sum_grad_C2 + eps);
			}
			else
				m_simulationData.getControlLambda(controlModelIndex, i) = 0.0;
		}
	}
}

void TimeStepPBF::pressureSolveIteration(const unsigned int fluidModelIndex, Real &avg_density_err)
{
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = sim->getFluidModel(fluidModelIndex);
	const unsigned int numParticles = model->numActiveParticles();
	const Real invH = static_cast<Real>(1.0) / TimeManager::getCurrent()->getTimeStepSize();
	const Real invH2 = invH*invH;
	const unsigned int nFluids = sim->numberOfFluidModels();
	//const unsigned int nControls = sim->numberOfControlBoundaryModels();
	const Real eps = 1.0e-6;

	const Real density0 = model->getDensity0();
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int) numParticles; i++)
		{
			Real &density = model->getDensity(i);				

			// Compute current density for particle i
			density = model->getVolume(i) * sim->W_zero();
			const Vector3r &xi = model->getPosition(i);

			////////////////////////////////////////////////////////////////////////
			// Fluid
			////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors(
				density += fm_neighbor->getVolume(neighborIndex) * sim->W(xi - xj);
			)
			//for (unsigned int pid = 0; pid < nFluids; pid++) 
			//{ 
			//	FluidModel* fm_neighbor = sim->getFluidModelFromPointSet(pid); 
			//	for (unsigned int j = 0; j < sim->numberOfNeighbors(fluidModelIndex, pid, i); j++) 
			//	{ 
			//		const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, pid, i, j); 
			//		const Vector3r& xj = fm_neighbor->getPosition(neighborIndex); 
			//		density += fm_neighbor->getVolume(neighborIndex) * sim->W(xi - xj);
			//	} 
			//}

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			forall_boundary_neighbors(
				// Boundary: Akinci2012
				density += bm_neighbor->getVolume(neighborIndex) * sim->W(xi - xj);
			)

			//for (unsigned int pid = nFluids; pid < sim->numberOfPointSets(); pid++)
			//{
			//	BoundaryModel* bm_neighbor = sim->getBoundaryModelFromPointSet(pid);
			//	for (unsigned int j = 0; j < sim->numberOfNeighbors(fluidModelIndex, pid, i); j++) 
			//	{ 
			//		const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, pid, i, j); 
			//		const Vector3r& xj = bm_neighbor->getPosition(neighborIndex); 
			//		density += bm_neighbor->getVolume(neighborIndex) * sim->W(xi - xj);
			//	} 
			//}

			const Real density_err = density0 * (max(density, static_cast<Real>(1.0)) - static_cast<Real>(1.0));
			#pragma omp atomic
			avg_density_err += density_err;

			// Evaluate constraint function
			const Real C = std::max(density - static_cast<Real>(1.0), static_cast<Real>(0.0));			// clamp to prevent particle clumping at surface

			if (C != 0.0)
			{
				// Compute gradients dC/dx_j 
				Real sum_grad_C2 = 0.0;
				Vector3r gradC_i(0.0, 0.0, 0.0);

				//////////////////////////////////////////////////////////////////////////
				// Fluid
				//////////////////////////////////////////////////////////////////////////
				forall_fluid_neighbors(
					const Vector3r gradC_j = -fm_neighbor->getVolume(neighborIndex) * sim->gradW(xi - xj);
					sum_grad_C2 += gradC_j.squaredNorm();
					gradC_i -= gradC_j;
				)

				//////////////////////////////////////////////////////////////////////////
				// Boundary
				//////////////////////////////////////////////////////////////////////////
				forall_boundary_neighbors(
					// Boundary: Akinci2012
					const Vector3r gradC_j = -bm_neighbor->getVolume(neighborIndex) * sim->gradW(xi - xj);
					sum_grad_C2 += gradC_j.squaredNorm();
					gradC_i -= gradC_j;
				)
				/*for (unsigned int pid = nFluids; pid < sim->numberOfPointSets(); pid++)
				{
					BoundaryModel* bm_neighbor = sim->getBoundaryModelFromPointSet(pid);
					for (unsigned int j = 0; j < sim->numberOfNeighbors(fluidModelIndex, pid, i); j++)
					{
						const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, pid, i, j);
						const Vector3r& xj = bm_neighbor->getPosition(neighborIndex);
						const Vector3r gradC_j = -bm_neighbor->getVolume(neighborIndex) * sim->gradW(xi - xj);
						sum_grad_C2 += gradC_j.squaredNorm();
						gradC_i -= gradC_j;
					}
				}*/


				sum_grad_C2 += gradC_i.squaredNorm();

				// Compute lambda
				Real &lambda = m_simulationData.getLambda(fluidModelIndex, i);
				lambda = -C / (sum_grad_C2 + eps);
			}
			else
				m_simulationData.getLambda(fluidModelIndex, i) = 0.0;
		}
	}
	avg_density_err /= numParticles;

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			Vector3r &corr = m_simulationData.getDeltaX(fluidModelIndex, i);

			// Compute position correction
			corr.setZero();
			const Vector3r &xi = model->getPosition(i);

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors(
				const Vector3r gradC_j = -fm_neighbor->getVolume(neighborIndex) * sim->gradW(xi - xj);
				corr -= (m_simulationData.getLambda(fluidModelIndex, i) + (fm_neighbor->getDensity0() / density0) * m_simulationData.getLambda(pid, neighborIndex)) * gradC_j;
				)
				//for (unsigned int pid = 0; pid < nFluids; pid++)
				//{
				//	FluidModel* fm_neighbor = sim->getFluidModelFromPointSet(pid);
				//	for (unsigned int j = 0; j < sim->numberOfNeighbors(fluidModelIndex, pid, i); j++)
				//	{
				//		const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, pid, i, j);
				//		const Vector3r& xj = fm_neighbor->getPosition(neighborIndex);
				//		const Vector3r gradC_j = -fm_neighbor->getVolume(neighborIndex) * sim->gradW(xi - xj);
				//		corr -= (m_simulationData.getLambda(fluidModelIndex, i) + (fm_neighbor->getDensity0() / density0) * m_simulationData.getLambda(pid, neighborIndex)) * gradC_j;
				//	}
				//}
				//////////////////////////////////////////////////////////////////////////
				// Boundary
				//////////////////////////////////////////////////////////////////////////
				//	forall_boundary_neighbors(
				//		// Boundary: Akinci2012
				//		const Vector3r gradC_j = -bm_neighbor->getVolume(neighborIndex) * sim->gradW(xi - xj);
				//const Vector3r dx = 2.0 * m_simulationData.getLambda(fluidModelIndex, i) * gradC_j;
				//corr -= dx;

				//bm_neighbor->addForce(xj, model->getMass(i)* dx* invH2);
				//);
				for (unsigned int pid = nFluids; pid < sim->numberOfPointSets(); pid++)
				{
					BoundaryModel* bm_neighbor = sim->getBoundaryModelFromPointSet(pid);
					for (unsigned int j = 0; j < sim->numberOfNeighbors(fluidModelIndex, pid, i); j++)
					{
						const unsigned int neighborIndex = sim->getNeighbor(fluidModelIndex, pid, i, j);
						const Vector3r& xj = bm_neighbor->getPosition(neighborIndex);
						const Vector3r gradC_j = -bm_neighbor->getVolume(neighborIndex) * sim->gradW(xi - xj);
						const Vector3r dx = 2.0 * m_simulationData.getLambda(fluidModelIndex, i) * gradC_j;
						corr -= dx;
						bm_neighbor->addForce(xj, model->getMass(i)* dx* invH2);
					}
				}
		}
	}
	
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			model->getPosition(i) += m_simulationData.getDeltaX(fluidModelIndex, i);
		}
	}
}

void TimeStepPBF::controlConstraints()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++)
	{
		Simulation* sim = Simulation::getCurrent();
		FluidModel* model = sim->getFluidModel(fluidModelIndex);
		const unsigned int numParticles = model->numActiveParticles();
		updateSpring(fluidModelIndex);
		updateDensityControl(fluidModelIndex);
		updateVelocityControl(fluidModelIndex);

#pragma omp parallel default(shared)
		{
#pragma omp for schedule(static)  
			for (int i = 0; i < (int)numParticles; i++)
			{
				model->getPosition(i) += m_simulationData.getDeltaSp(fluidModelIndex, i);
				//m_simulationData.getDeltaSp(fluidModelIndex, i).setZero();
				model->getPosition(i) += m_simulationData.getDeltaDe(fluidModelIndex, i);
				//m_simulationData.getDeltaDe(fluidModelIndex, i).setZero();
				model->getPosition(i) += m_simulationData.getDeltaVe(fluidModelIndex, i);
				//model->getPosition(i) += m_simulationData.getDeltaX(fluidModelIndex, i);
			}
		}
	}

}

void TimeStepPBF::updateVelocityControl(const unsigned int fluidModelIndex)
{
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	const unsigned int numParticles = model->numActiveParticles();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundarys = sim->numberOfBoundaryModels();
	const Real h = TimeManager::getCurrent()->getTimeStepSize();
	const Real supportradius = sim->getSupportRadius();
	const Real gammaValue = TimeStepPBF::getValue<Real>(TimeStepPBF::VELOCITY_GAMMA);
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r& vi = model->getVelocity(i);
			const Vector3r& xi = model->getPosition(i);
			Vector3r& corr = m_simulationData.getDeltaVe(fluidModelIndex, i);
			corr.setZero();
			for (unsigned int pid = nFluids; pid < sim->cnumberOfPointSets(); pid++)
			{
				ControlBoundaryModel* cbm_neighbor = sim->getControlBoundaryModelFromPointSet(pid);
				const int num_nei = sim->cnumberOfNeighbors(fluidModelIndex, pid, i);
				const int numOfNeighbors = sim->cnumberOfNeighbors(fluidModelIndex, pid, i);
				for (unsigned int j = 0; j < numOfNeighbors; j++)
				{
					const unsigned int neighborIndex = sim->cgetNeighbor(fluidModelIndex, pid, i, j);
					const Vector3r& xj = cbm_neighbor->getPosition(neighborIndex);
					//if ((xi - xj).norm() <= supportradius)
					//{
						//if (numOfNeighbors > 0 && numOfNeighbors <= 10)
						//	corr = -0.1 * h * vi;
						//else if (numOfNeighbors > 10 && numOfNeighbors <= 20)
						//	corr = -1.0 * h * vi;
						//else if (numOfNeighbors > 20)
						//	corr = -1.0 * h * vi;
						//else
							corr = -gammaValue * h * vi;
					//}
				}
			}
			
		}
	}
}

void TimeStepPBF::updateDensityControl(const unsigned int fluidModelIndex)
{
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	const unsigned int numParticles = model->numActiveParticles();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundarys = sim->numberOfBoundaryModels();
	const Real fluidDensity0 = model->getDensity0();
	const Real supportradius = sim->getSupportRadius();
	const Real alphaValue = TimeStepPBF::getValue<Real>(TimeStepPBF::DENSITY_ALPHA);
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r& xi = model->getPosition(i);

			// Compute position correction
			Vector3r& corr = m_simulationData.getDeltaDe(fluidModelIndex, i);
			corr.setZero();
			
			for (unsigned int pid = nFluids; pid < sim->cnumberOfPointSets(); pid++) 
			{ 
				
				ControlBoundaryModel* cbm_neighbor = sim->getControlBoundaryModelFromPointSet(pid); 
				for (unsigned int j = 0; j < sim->cnumberOfNeighbors(fluidModelIndex, pid, i); j++) 
				{ 
					const unsigned int neighborIndex = sim->cgetNeighbor(fluidModelIndex, pid, i, j); 
					const Vector3r& xj = cbm_neighbor->getPosition(neighborIndex);
					//if ((xi - xj).norm() <= supportradius)
					//{
						const Vector3r gradC_j = sim->gradCW(xj - xi);
						//const Vector3r gradC_j = cbm_neighbor->getVolume(neighborIndex) * sim->gradW(xj - xi);
						//corr += alphaValue/ cbm_neighbor->getDensity0() * m_simulationData.getControlLambda(pid - nFluids, neighborIndex) * gradC_j;//2.0代表alpha的值
						corr -= (alphaValue / 1000 )* m_simulationData.getControlLambda(pid - nFluids, neighborIndex) * gradC_j;//暂时÷流体1000的密度，按道理是控制粒子密度
					//}
					
				}

			}
		}
	}

}

void TimeStepPBF::updateSpring(const unsigned int fluidModelIndex)
{
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	const unsigned int numParticles = model->numActiveParticles();
	const Real supportradius = sim->getCSupportRadius();
	//const Real invH = static_cast<Real>(1.0) / TimeManager::getCurrent()->getTimeStepSize();
	//const Real invH2 = invH * invH;
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundarys = sim->numberOfBoundaryModels();
	const Real betaValue = TimeStepPBF::getValue<Real>(TimeStepPBF::SPRING_BETA);
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r& xi = model->getPosition(i);

			// Compute position correction:spring
			Vector3r &corr = m_simulationData.getDeltaSp(fluidModelIndex, i);
			corr.setZero();

			int& ri = model->getri(i);
			int ri_tmp = 0;
			int& timelife = model->getTimelife(i);
			bool chk = false;
			for (unsigned int pid = nFluids; pid < sim->cnumberOfPointSets(); pid++)
			{ 
				ControlBoundaryModel* cbm_neighbor = sim->getControlBoundaryModelFromPointSet(pid);
				Real ri_dis = supportradius;
				for (unsigned int j = 0; j < sim->cnumberOfNeighbors(fluidModelIndex, pid, i); j++) 
				{ 

					const unsigned int neighborIndex = sim->cgetNeighbor(fluidModelIndex, pid, i, j);
					const Vector3r& xj = cbm_neighbor->getPosition(neighborIndex);
					const Real xj_dis = (xi - xj).norm();
					if (xj_dis <= ri_dis)// && timelife==0)
					{
						ri_dis = xj_dis;
						ri_tmp = neighborIndex;
						chk = true;
					}
				}
			}

			for (unsigned int pid = nFluids; pid < sim->cnumberOfPointSets(); pid++)
			{
				ControlBoundaryModel* cbm_neughbor = sim->getControlBoundaryModelFromPointSet(pid);
				
				const Vector3r& xj = cbm_neughbor->getPosition(ri_tmp);
				const Real dis = (xi - xj).norm();
				if (chk && dis < supportradius)// d/2
				{
					timelife = (rand() % 3 + 1);
					ri = ri_tmp;
				}
				
				
				if (timelife >0)
				{
					const Real& Ci = cbm_neughbor->getC(ri);
					const Vector3r xj_ri = cbm_neughbor->getPosition(ri);
					//const Real& Ci = cbm_neughbor->getC(ri);
					const Real dist = (xj_ri - xi).norm();
					timelife = timelife - 1;
					corr = betaValue * (xj_ri - xi) / dist / (1 + exp(80 * (supportradius - dist)));//2.0 代表调整弹簧力系数Beta
				}
				else
				{
					corr.setZero();
					ri = -1;
				}
			}
		}
	}
}
void TimeStepPBF::performNeighborhoodSearch()
{
	if (m_counter % 500 == 0)
	{
		Simulation::getCurrent()->performNeighborhoodSearchSort();
		m_simulationData.performNeighborhoodSearchSort();
	}
	m_counter++;

	Simulation::getCurrent()->performNeighborhoodSearch();
}

void TimeStepPBF::emittedParticles(FluidModel *model, const unsigned int startIndex)
{
	m_simulationData.emittedParticles(model, startIndex);
}

void TimeStepPBF::resize()
{
	m_simulationData.init();
}

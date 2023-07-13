#include "Simulation.h"
#include "TimeManager.h"
#include "Utilities/Timing.h"
#include "TimeStep.h"
#include "EmitterSystem.h"
#include "SPlisHSPlasH/WCSPH/TimeStepWCSPH.h"
#include "SPlisHSPlasH/PCISPH/TimeStepPCISPH.h"
#include "SPlisHSPlasH/PBF/TimeStepPBF.h"
#include "SPlisHSPlasH/IISPH/TimeStepIISPH.h"
#include "SPlisHSPlasH/DFSPH/TimeStepDFSPH.h"
#include "SPlisHSPlasH/PF/TimeStepPF.h"
#include <set>



using namespace SPH;
using namespace std;
using namespace GenParam;

Simulation* Simulation::current = nullptr;
int Simulation::SIM_2D = -1;
int Simulation::PARTICLE_RADIUS = -1;
int Simulation::GRAVITATION = -1;
int Simulation::START_EMIT_AT = -1;
int Simulation::STOP_EMIT_AT = -1;
int Simulation::CFL_METHOD = -1;
int Simulation::CFL_FACTOR = -1;
int Simulation::CFL_MAX_TIMESTEPSIZE = -1;
int Simulation::KERNEL_METHOD = -1;
int Simulation::GRAD_KERNEL_METHOD = -1;
int Simulation::ENUM_KERNEL_CUBIC = -1;
int Simulation::ENUM_KERNEL_WENDLANDQUINTICC2 = -1;
int Simulation::ENUM_KERNEL_POLY6 = -1;
int Simulation::ENUM_KERNEL_SPIKY = -1;
int Simulation::ENUM_KERNEL_PRECOMPUTED_CUBIC = -1;
int Simulation::ENUM_KERNEL_CUBIC_2D = -1;
int Simulation::ENUM_KERNEL_WENDLANDQUINTICC2_2D = -1;
int Simulation::ENUM_GRADKERNEL_CUBIC = -1;
int Simulation::ENUM_GRADKERNEL_WENDLANDQUINTICC2 = -1;
int Simulation::ENUM_GRADKERNEL_POLY6 = -1;
int Simulation::ENUM_GRADKERNEL_SPIKY = -1;
int Simulation::ENUM_GRADKERNEL_PRECOMPUTED_CUBIC = -1;
int Simulation::ENUM_GRADKERNEL_CUBIC_2D = -1;
int Simulation::ENUM_GRADKERNEL_WENDLANDQUINTICC2_2D = -1;
int Simulation::SIMULATION_METHOD = -1;
int Simulation::ENUM_CFL_NONE = -1;
int Simulation::ENUM_CFL_STANDARD = -1;
int Simulation::ENUM_CFL_ITER = -1;
int Simulation::ENUM_SIMULATION_WCSPH = -1;
int Simulation::ENUM_SIMULATION_PCISPH = -1;
int Simulation::ENUM_SIMULATION_PBF = -1;
int Simulation::ENUM_SIMULATION_IISPH = -1;
int Simulation::ENUM_SIMULATION_DFSPH = -1;
int Simulation::ENUM_SIMULATION_PF = -1;


Simulation::Simulation () 
{
	m_cflMethod = 1;
	stopEm = -1;
	startEm = -1;
	m_cflFactor = 0.5;
	m_cflMaxTimeStepSize = 0.005;
	m_gravitation = Vector3r(0.0, -9.81, 0.0);

	m_kernelMethod = -1;
	m_gradKernelMethod = -1;

	m_neighborhoodSearch = nullptr;
	if(_CONTROL)
		c_neighborhoodSearch = nullptr;
	//if(_TARGET&_CONTROL)
	//	t_neighborhoodSearch = nullptr;
	m_timeStep = nullptr;
	m_simulationMethod = SimulationMethods::NumSimulationMethods;
	m_simulationMethodChanged = NULL;

	m_sim2D = false;
	if_control = _CONTROL;
	//if_target = _TARGET && _CONTROL;
}

Simulation::~Simulation () 
{
	delete m_timeStep;
	delete m_neighborhoodSearch;
	delete c_neighborhoodSearch;
	delete t_neighborhoodSearch;
	delete TimeManager::getCurrent();

	for (unsigned int i = 0; i < m_fluidModels.size(); i++)
		delete m_fluidModels[i];
	m_fluidModels.clear();

	for (unsigned int i = 0; i < m_boundaryModels.size(); i++)
		delete m_boundaryModels[i];
	m_boundaryModels.clear();

	for (unsigned int i = 0; i < m_controlBoundaryModels.size(); i++)
		delete m_controlBoundaryModels[i];
	m_controlBoundaryModels.clear();

	for (unsigned int i = 0; i < m_targetModels.size(); i++)
		delete m_targetModels[i];
	m_targetModels.clear();

	current = nullptr;
}

Simulation* Simulation::getCurrent ()
{
	if (current == nullptr)
	{
		current = new Simulation ();
	}
	return current;
}

void Simulation::setCurrent (Simulation* tm)
{
	current = tm;
}

bool Simulation::hasCurrent()
{
	return (current != nullptr);
}

void Simulation::init(const Real particleRadius, const bool sim2D)
{
	m_sim2D = sim2D;
	initParameters();



	// init kernel
	setParticleRadius(particleRadius);
	
	setValue(Simulation::KERNEL_METHOD, Simulation::ENUM_KERNEL_PRECOMPUTED_CUBIC);
	setValue(Simulation::GRAD_KERNEL_METHOD, Simulation::ENUM_GRADKERNEL_PRECOMPUTED_CUBIC);

	// Initialize neighborhood search
	if (m_neighborhoodSearch == NULL)
#ifdef GPU_NEIGHBORHOOD_SEARCH
		m_neighborhoodSearch = new NeighborhoodSearch(m_supportRadius);
#else
		m_neighborhoodSearch = new NeighborhoodSearch(m_supportRadius, false);
#endif
	m_neighborhoodSearch->set_radius(m_supportRadius);
	if (c_neighborhoodSearch == NULL && ifControl())
	{
#ifdef GPU_NEIGHBORHOOD_SEARCH
		c_neighborhoodSearch = new NeighborhoodSearch(c_supportRadius);
#else
		c_neighborhoodSearch = new NeighborhoodSearch(c_supportRadius, false);
#endif
		c_neighborhoodSearch->set_radius(c_supportRadius);
	}
	
	//if (t_neighborhoodSearch == NULL && ifTarget())
	//{
	//	t_neighborhoodSearch = new NeighborhoodSearch(m_supportRadius, false);
	//	t_neighborhoodSearch->set_radius(m_supportRadius);
	//}
		
}

void Simulation::initParameters()
{
	ParameterObject::initParameters();

	SIM_2D = createBoolParameter("sim2D", "2D Simulation", &m_sim2D);
	setGroup(SIM_2D, "Simulation");
	setDescription(SIM_2D, "2D/3D simulation.");
	getParameter(SIM_2D)->setReadOnly(true);

	ParameterBase::GetFunc<Real> getRadiusFct = std::bind(&Simulation::getParticleRadius, this);
	ParameterBase::SetFunc<Real> setRadiusFct = std::bind(&Simulation::setParticleRadius, this, std::placeholders::_1);
	PARTICLE_RADIUS = createNumericParameter("particleRadius", "Particle radius", getRadiusFct, setRadiusFct);
	setGroup(PARTICLE_RADIUS, "Simulation");
	setDescription(PARTICLE_RADIUS, "Radius of the fluid particles.");
	getParameter(PARTICLE_RADIUS)->setReadOnly(true);

 	GRAVITATION = createVectorParameter("gravitation", "Gravitation", 3u, m_gravitation.data());
 	setGroup(GRAVITATION, "Simulation");
 	setDescription(GRAVITATION, "Vector to define the gravitational acceleration.");

	START_EMIT_AT = createNumericParameter("startEmit", "start Emit at", &startEm);
	setGroup(START_EMIT_AT, "Simulation");
	setDescription(START_EMIT_AT, "start Emit Particles at time");

	STOP_EMIT_AT = createNumericParameter("stopEmit", "Stop Emit at", &stopEm);
	setGroup(STOP_EMIT_AT, "Simulation");
	setDescription(STOP_EMIT_AT, "Stop Emit Particles at time");

	CFL_METHOD = createEnumParameter("cflMethod", "CFL - method", &m_cflMethod);
	setGroup(CFL_METHOD, "CFL");
	setDescription(CFL_METHOD, "CFL method used for adaptive time stepping.");
	EnumParameter *enumParam = static_cast<EnumParameter*>(getParameter(CFL_METHOD));
	enumParam->addEnumValue("None", ENUM_CFL_NONE);
	enumParam->addEnumValue("CFL", ENUM_CFL_STANDARD);
	enumParam->addEnumValue("CFL - iterations", ENUM_CFL_ITER);

	CFL_FACTOR = createNumericParameter("cflFactor", "CFL - factor", &m_cflFactor);
	setGroup(CFL_FACTOR, "CFL");
	setDescription(CFL_FACTOR, "Factor to scale the CFL time step size.");
	static_cast<RealParameter*>(getParameter(CFL_FACTOR))->setMinValue(1e-6);

	CFL_MAX_TIMESTEPSIZE = createNumericParameter("cflMaxTimeStepSize", "CFL - max. time step size", &m_cflMaxTimeStepSize);
	setGroup(CFL_MAX_TIMESTEPSIZE, "CFL");
	setDescription(CFL_MAX_TIMESTEPSIZE, "Max. time step size.");
	static_cast<RealParameter*>(getParameter(CFL_MAX_TIMESTEPSIZE))->setMinValue(1e-6);

	ParameterBase::GetFunc<int> getKernelFct = std::bind(&Simulation::getKernel, this);
	ParameterBase::SetFunc<int> setKernelFct = std::bind(&Simulation::setKernel, this, std::placeholders::_1);
	KERNEL_METHOD = createEnumParameter("kernel", "Kernel", getKernelFct, setKernelFct);
	setGroup(KERNEL_METHOD, "Kernel");
	setDescription(KERNEL_METHOD, "Kernel function used in the SPH model.");
	enumParam = static_cast<EnumParameter*>(getParameter(KERNEL_METHOD));
	if (!m_sim2D)
	{
		enumParam->addEnumValue("Cubic spline", ENUM_KERNEL_CUBIC);
		enumParam->addEnumValue("Wendland quintic C2", ENUM_KERNEL_WENDLANDQUINTICC2);
		enumParam->addEnumValue("Poly6", ENUM_KERNEL_POLY6);
		enumParam->addEnumValue("Spiky", ENUM_KERNEL_SPIKY);
		enumParam->addEnumValue("Precomputed cubic spline", ENUM_KERNEL_PRECOMPUTED_CUBIC);
	}
	else
	{
		enumParam->addEnumValue("Cubic spline 2D", ENUM_KERNEL_CUBIC_2D);
		enumParam->addEnumValue("Wendland quintic C2 2D", ENUM_KERNEL_WENDLANDQUINTICC2_2D);
	}

	ParameterBase::GetFunc<int> getGradKernelFct = std::bind(&Simulation::getGradKernel, this);
	ParameterBase::SetFunc<int> setGradKernelFct = std::bind(&Simulation::setGradKernel, this, std::placeholders::_1);
	GRAD_KERNEL_METHOD = createEnumParameter("gradKernel", "Gradient of kernel", getGradKernelFct, setGradKernelFct);
	setGroup(GRAD_KERNEL_METHOD, "Kernel");
	setDescription(GRAD_KERNEL_METHOD, "Gradient of the kernel function used in the SPH model.");
	enumParam = static_cast<EnumParameter*>(getParameter(GRAD_KERNEL_METHOD));
	if (!m_sim2D)
	{
		enumParam->addEnumValue("Cubic spline", ENUM_GRADKERNEL_CUBIC);
		enumParam->addEnumValue("Wendland quintic C2", ENUM_GRADKERNEL_WENDLANDQUINTICC2);
		enumParam->addEnumValue("Poly6", ENUM_GRADKERNEL_POLY6);
		enumParam->addEnumValue("Spiky", ENUM_GRADKERNEL_SPIKY);
		enumParam->addEnumValue("Precomputed cubic spline", ENUM_GRADKERNEL_PRECOMPUTED_CUBIC);
	}
	else
	{
		enumParam->addEnumValue("Cubic spline 2D", ENUM_GRADKERNEL_CUBIC_2D);
		enumParam->addEnumValue("Wendland quintic C2 2D", ENUM_GRADKERNEL_WENDLANDQUINTICC2_2D);
	}

	ParameterBase::GetFunc<int> getSimulationFct = std::bind(&Simulation::getSimulationMethod, this);
	ParameterBase::SetFunc<int> setSimulationFct = std::bind(&Simulation::setSimulationMethod, this, std::placeholders::_1);
	SIMULATION_METHOD = createEnumParameter("simulationMethod", "Simulation method", getSimulationFct, setSimulationFct);
	setGroup(SIMULATION_METHOD, "Simulation");
	setDescription(SIMULATION_METHOD, "Simulation method.");
	enumParam = static_cast<EnumParameter*>(getParameter(SIMULATION_METHOD));
	enumParam->addEnumValue("WCSPH", ENUM_SIMULATION_WCSPH);
	enumParam->addEnumValue("PCISPH", ENUM_SIMULATION_PCISPH);
	enumParam->addEnumValue("PBF", ENUM_SIMULATION_PBF);
	enumParam->addEnumValue("IISPH", ENUM_SIMULATION_IISPH);
	enumParam->addEnumValue("DFSPH", ENUM_SIMULATION_DFSPH);
	enumParam->addEnumValue("Projective Fluids", ENUM_SIMULATION_PF);
}


void Simulation::setParticleRadius(Real val)
{
	m_particleRadius = val;
	m_supportRadius = static_cast<Real>(4.0)*m_particleRadius;//0.1
	c_supportRadius = 2.2 * m_supportRadius;//0.22
	// init kernel
	Poly6Kernel::setRadius(m_supportRadius);
	SpikyKernel::setRadius(m_supportRadius);
	CubicKernel::setRadius(m_supportRadius);
	WendlandQuinticC2Kernel::setRadius(m_supportRadius);
	PrecomputedCubicKernel::setRadius(m_supportRadius);
	CohesionKernel::setRadius(m_supportRadius);
	AdhesionKernel::setRadius(m_supportRadius);
	CubicKernel2D::setRadius(m_supportRadius);
	WendlandQuinticC2Kernel2D::setRadius(m_supportRadius);

	//init control kernel
	Poly6Kernel::setCRadius(c_supportRadius);
	SpikyKernel::setCRadius(c_supportRadius);
	CubicKernel::setCRadius(c_supportRadius);
	WendlandQuinticC2Kernel::setCRadius(c_supportRadius);
	PrecomputedCubicKernel::setCRadius(c_supportRadius);
	CohesionKernel::setCRadius(c_supportRadius);
	AdhesionKernel::setCRadius(c_supportRadius);
	CubicKernel2D::setCRadius(c_supportRadius);
	WendlandQuinticC2Kernel2D::setCRadius(c_supportRadius);
}

void Simulation::setGradKernel(int val)
{
	m_gradKernelMethod = val;

	if (!m_sim2D)
	{
		if ((m_gradKernelMethod < 0) || (m_gradKernelMethod > 4))
			m_gradKernelMethod = 0;

		if (m_gradKernelMethod == 0)
		{
			m_gradKernelFct = CubicKernel::gradW;
			c_gradKernelFct = CubicKernel::gradCW;
		}
		else if (m_gradKernelMethod == 1)
		{
			m_gradKernelFct = WendlandQuinticC2Kernel::gradW;
			c_gradKernelFct = WendlandQuinticC2Kernel::gradCW;
		}
		else if (m_gradKernelMethod == 2)
		{
			m_gradKernelFct = Poly6Kernel::gradW;
			c_gradKernelFct = Poly6Kernel::gradCW;
		}
		else if (m_gradKernelMethod == 3)
		{
			m_gradKernelFct = SpikyKernel::gradW;
			c_gradKernelFct = SpikyKernel::gradCW;
		}
		else if (m_gradKernelMethod == 4)
		{
			m_gradKernelFct = Simulation::PrecomputedCubicKernel::gradW;
			c_gradKernelFct = Simulation::PrecomputedCubicKernel::gradCW;
		}
			
	}
	else
	{
		if ((m_gradKernelMethod < 0) || (m_gradKernelMethod > 1))
			m_gradKernelMethod = 0;

		if (m_gradKernelMethod == 0)
		{
			m_gradKernelFct = CubicKernel2D::gradW;
			c_gradKernelFct = CubicKernel2D::gradCW;
		}
			
		else if (m_gradKernelMethod == 1)
		{
			m_gradKernelFct = WendlandQuinticC2Kernel2D::gradW;
			c_gradKernelFct = WendlandQuinticC2Kernel2D::gradCW;
		}
			
	}
}

void Simulation::setKernel(int val)
{
	if (val == m_kernelMethod)
		return;

	m_kernelMethod = val;
	if (!m_sim2D)
	{
		if ((m_kernelMethod < 0) || (m_kernelMethod > 4))
			m_kernelMethod = 0;

		if (m_kernelMethod == 0)
		{
			m_W_zero = CubicKernel::W_zero();
			m_kernelFct = CubicKernel::W;
			c_W_zero = CubicKernel::CW_zero();
			c_kernelFct = CubicKernel::CW;
		}
		else if (m_kernelMethod == 1)
		{
			m_W_zero = WendlandQuinticC2Kernel::W_zero();
			m_kernelFct = WendlandQuinticC2Kernel::W;
			c_W_zero = WendlandQuinticC2Kernel::CW_zero();
			c_kernelFct = WendlandQuinticC2Kernel::CW;
		}
		else if (m_kernelMethod == 2)
		{
			m_W_zero = Poly6Kernel::W_zero();
			m_kernelFct = Poly6Kernel::W;
			c_W_zero = Poly6Kernel::CW_zero();
			c_kernelFct = Poly6Kernel::CW;
		}
		else if (m_kernelMethod == 3)
		{
			m_W_zero = SpikyKernel::W_zero();
			m_kernelFct = SpikyKernel::W;
			c_W_zero = SpikyKernel::CW_zero();
			c_kernelFct = SpikyKernel::CW;
		}
		else if (m_kernelMethod == 4)
		{
			m_W_zero = Simulation::PrecomputedCubicKernel::W_zero();
			m_kernelFct = Simulation::PrecomputedCubicKernel::W;
			c_W_zero = Simulation::PrecomputedCubicKernel::CW_zero();
			c_kernelFct = Simulation::PrecomputedCubicKernel::CW;
		}
	}
	else
	{
		if ((m_kernelMethod < 0) || (m_kernelMethod > 1))
			m_kernelMethod = 0;

		if (m_kernelMethod == 0)
		{
			m_W_zero = CubicKernel2D::W_zero();
			m_kernelFct = CubicKernel2D::W;
			c_W_zero = CubicKernel2D::CW_zero();
			c_kernelFct = CubicKernel2D::CW;
		}
		else if (m_kernelMethod == 1)
		{
			m_W_zero = WendlandQuinticC2Kernel2D::W_zero();
			m_kernelFct = WendlandQuinticC2Kernel2D::W;
			c_W_zero = WendlandQuinticC2Kernel2D::CW_zero();
			c_kernelFct = WendlandQuinticC2Kernel2D::CW;
		}
	}
	updateBoundaryVolume();
	if(_CONTROL)
		updateControlVloume();

	
	//ActivateOnlyFluids();
}

void Simulation::updateTimeStepSize()
{
	if (m_cflMethod == 1)
		updateTimeStepSizeCFL(0.0001);
	else if (m_cflMethod == 2)
	{
		Real h = TimeManager::getCurrent()->getTimeStepSize();
		updateTimeStepSizeCFL(0.0001);
		const unsigned int iterations = m_timeStep->getValue<unsigned int>(TimeStep::SOLVER_ITERATIONS);
		if (iterations > 10)
			h *= 0.9;
		else if (iterations < 5)
			h *= 1.1;
		h = min(h, TimeManager::getCurrent()->getTimeStepSize());
		TimeManager::getCurrent()->setTimeStepSize(h);
	}
}

void Simulation::updateTimeStepSizeCFL(const Real minTimeStepSize)
{
	const Real radius = m_particleRadius;
	Real h = TimeManager::getCurrent()->getTimeStepSize();

	// Approximate max. position change due to current velocities
	Real maxVel = 0.1;
	const Real diameter = static_cast<Real>(2.0)*radius;

	// fluid particles
	for (unsigned int i = 0; i < numberOfFluidModels(); i++)
	{
		FluidModel *fm = getFluidModel(i);
		const unsigned int numParticles = fm->numActiveParticles();
		for (unsigned int i = 0; i < numParticles; i++)
		{
			const Vector3r &vel = fm->getVelocity(i);
			const Vector3r &accel = fm->getAcceleration(i);
			const Real velMag = (vel + accel*h).squaredNorm();
			if (velMag > maxVel)
				maxVel = velMag;
		}
	}

	// boundary particles
	for (unsigned int i = 0; i < numberOfBoundaryModels(); i++)
	{
		BoundaryModel *bm = getBoundaryModel(i);
		if (bm->getRigidBodyObject()->isDynamic())
		{
			for (unsigned int j = 0; j < bm->numberOfParticles(); j++)
			{
				const Vector3r &vel = bm->getVelocity(j);
				const Real velMag = vel.squaredNorm();
				if (velMag > maxVel)
					maxVel = velMag;
			}
		}
	}


	// Approximate max. time step size 		
	h = m_cflFactor * static_cast<Real>(0.4) * (diameter / (sqrt(maxVel)));

	h = min(h, m_cflMaxTimeStepSize);
	h = max(h, minTimeStepSize);

	TimeManager::getCurrent()->setTimeStepSize(h);
}

void Simulation::computeNonPressureForces()
{
	START_TIMING("computeNonPressureForces")
	for (unsigned int i = 0; i < numberOfFluidModels(); i++)
	{
		FluidModel *fm = getFluidModel(i);
		fm->computeSurfaceTension();
		fm->computeViscosity();
		fm->computeVorticity();
		fm->computeDragForce();
		fm->computeElasticity();
	}
	STOP_TIMING_AVG
}

void Simulation::reset()
{
	// reset fluid models
	for (unsigned int i = 0; i < numberOfFluidModels(); i++)
		getFluidModel(i)->reset();

	// reset boundary models
	for (unsigned int i = 0; i < numberOfBoundaryModels(); i++)
		getBoundaryModel(i)->reset();
	for (unsigned int i = 0; i < numberOfControlBoundaryModels(); i++)
		getControlBoundaryModel(i)->reset();
	updateBoundaryVolume();
	if(_CONTROL)
		updateControlVloume();
	ActivateNeighborFun();


	if (m_timeStep)
		m_timeStep->reset();

	performNeighborhoodSearchSort();

	TimeManager::getCurrent()->setTime(0.0);
	//TimeManager::getCurrent()->setTimeStepSize(0.001);
}

void Simulation::setSimulationMethod(const int val)
{
	SimulationMethods method = static_cast<SimulationMethods>(val);
	if ((method < SimulationMethods::WCSPH) || (method >= SimulationMethods::NumSimulationMethods))
		method = SimulationMethods::PBF;


	if (method == m_simulationMethod)
		return;

	delete m_timeStep;
	m_timeStep = nullptr;

	m_simulationMethod = method;

	if (method == SimulationMethods::WCSPH)
	{
		m_timeStep = new TimeStepWCSPH();
		m_timeStep->init();
		setValue(Simulation::CFL_METHOD, Simulation::ENUM_CFL_NONE);
		setValue(Simulation::KERNEL_METHOD, Simulation::ENUM_KERNEL_CUBIC);
		setValue(Simulation::GRAD_KERNEL_METHOD, Simulation::ENUM_GRADKERNEL_CUBIC);
		//TimeManager::getCurrent()->setTimeStepSize(0.001);
	}
	else if (method == SimulationMethods::PCISPH)
	{
		m_timeStep = new TimeStepPCISPH();
		m_timeStep->init();
		setValue(Simulation::CFL_METHOD, Simulation::ENUM_CFL_STANDARD);
		setValue(Simulation::KERNEL_METHOD, Simulation::ENUM_KERNEL_CUBIC);
		setValue(Simulation::GRAD_KERNEL_METHOD, Simulation::ENUM_GRADKERNEL_CUBIC);
	}
	else if (method == SimulationMethods::PBF)
	{
		m_timeStep = new TimeStepPBF();
		m_timeStep->init();
		setValue(Simulation::CFL_METHOD, Simulation::ENUM_CFL_STANDARD);
		//setValue(Simulation::KERNEL_METHOD, Simulation::ENUM_KERNEL_POLY6);
		//setValue(Simulation::GRAD_KERNEL_METHOD, Simulation::ENUM_GRADKERNEL_SPIKY);
		setValue(Simulation::KERNEL_METHOD, Simulation::ENUM_KERNEL_PRECOMPUTED_CUBIC);
		setValue(Simulation::GRAD_KERNEL_METHOD, Simulation::ENUM_GRADKERNEL_PRECOMPUTED_CUBIC);
	}
	else if (method == SimulationMethods::IISPH)
	{
		m_timeStep = new TimeStepIISPH();
		m_timeStep->init();
		setValue(Simulation::CFL_METHOD, Simulation::ENUM_CFL_STANDARD);
		setValue(Simulation::KERNEL_METHOD, Simulation::ENUM_KERNEL_CUBIC);
		setValue(Simulation::GRAD_KERNEL_METHOD, Simulation::ENUM_GRADKERNEL_CUBIC);
	}
	else if (method == SimulationMethods::DFSPH)
	{
		m_timeStep = new TimeStepDFSPH();
		m_timeStep->init();
		setValue(Simulation::CFL_METHOD, Simulation::ENUM_CFL_STANDARD);
		setValue(Simulation::KERNEL_METHOD, Simulation::ENUM_KERNEL_PRECOMPUTED_CUBIC);
		setValue(Simulation::GRAD_KERNEL_METHOD, Simulation::ENUM_GRADKERNEL_PRECOMPUTED_CUBIC);
	}
	else if (method == SimulationMethods::PF)
	{
		m_timeStep = new TimeStepPF();
		m_timeStep->init();
		setValue(Simulation::CFL_METHOD, Simulation::ENUM_CFL_STANDARD);
		setValue(Simulation::KERNEL_METHOD, Simulation::ENUM_KERNEL_PRECOMPUTED_CUBIC);
		setValue(Simulation::GRAD_KERNEL_METHOD, Simulation::ENUM_GRADKERNEL_PRECOMPUTED_CUBIC);
	}

	if (m_simulationMethodChanged != nullptr)
		m_simulationMethodChanged();
}



void Simulation::performNeighborhoodSearch()
{
	START_TIMING("neighborhood_search");
	m_neighborhoodSearch->find_neighbors();
	if(_CONTROL)
		c_neighborhoodSearch->find_neighbors();
	//if(_TARGET)
	//	t_neighborhoodSearch->find_neighbors();
	STOP_TIMING_AVG;
}

void Simulation::performNeighborhoodSearchSort()
{
	m_neighborhoodSearch->z_sort();
	if(_CONTROL)
		c_neighborhoodSearch->z_sort();

	for (unsigned int i = 0; i < numberOfFluidModels(); i++)
	{
		FluidModel *fm = getFluidModel(i);
		fm->performNeighborhoodSearchSort();
	}
}

void Simulation::setSimulationMethodChangedCallback(std::function<void()> const& callBackFct)
{
	m_simulationMethodChanged = callBackFct;
}

void Simulation::emittedParticles(FluidModel *model, const unsigned int startIndex)
{
	model->emittedParticles(startIndex);
	m_timeStep->emittedParticles(model, startIndex);
}

void Simulation::emitParticles()
{
	for (unsigned int i = 0; i < numberOfFluidModels(); i++)
	{
		FluidModel *fm = getFluidModel(i);
		fm->getEmitterSystem()->step();
	}
}

void Simulation::stopEmit()
{
	for (unsigned int i = 0; i < numberOfFluidModels(); i++)
	{
		FluidModel* fm = getFluidModel(i);
		fm->getEmitterSystem()->disableReuseParticles();
	}
}

void Simulation::addTargetModel(const std::string& id, const unsigned int nTargetParticles, Vector3r* targetParticles)
{
	TargetModel* tm = new TargetModel();
	tm->initModel(id, nTargetParticles, targetParticles);
	m_targetModels.push_back(tm);
}

void Simulation::addControlModelf(const std::string& id, const unsigned int nControlParticles, Vector3r* controlParticles)
{
	ControlBoundaryModel* cbm = new ControlBoundaryModel();
	cbm->initModel(id, nControlParticles, controlParticles);
	m_controlBoundaryModels.push_back(cbm);
}

void Simulation::addControlModel(RigidBodyObject* rbo,unsigned int i, const unsigned int numControlParticles, Vector3r* controlParticles)
{
	ControlBoundaryModel* cbm = new ControlBoundaryModel();
	cbm->initModel(rbo,i, numControlParticles, controlParticles);
	m_controlBoundaryModels.push_back(cbm);
}

void Simulation::addBoundaryModel(RigidBodyObject *rbo, const unsigned int numBoundaryParticles, Vector3r *boundaryParticles)
{
	BoundaryModel *bm = new BoundaryModel();
	bm->initModel(rbo, numBoundaryParticles, boundaryParticles);
	m_boundaryModels.push_back(bm);
}

void Simulation::addFluidModel(const std::string &id, const unsigned int nFluidParticles, Vector3r* fluidParticles, Vector3r* fluidVelocities, const unsigned int nMaxEmitterParticles)
{
	FluidModel *fm = new FluidModel();
	fm->initModel(id, nFluidParticles, fluidParticles, fluidVelocities, nMaxEmitterParticles);
	m_fluidModels.push_back(fm);
}

void Simulation::updateTargetSearch()
{
	if (t_neighborhoodSearch == nullptr)
		return;
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nControls = sim->numberOfControlBoundaryModels();
	LOG_INFO << "Initialize target search";
	t_neighborhoodSearch->set_active(false);
	for (unsigned int i = 0; i < numberOfTargetModels(); i++)
	{
		t_neighborhoodSearch->set_active(i + nControls, true, true);
	}
	t_neighborhoodSearch->find_neighbors();
}

void Simulation::updateControlVloume()
{
	if (c_neighborhoodSearch == nullptr)
		return;
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBonudarys = sim->numberOfBoundaryModels();
	// Activate only static boundaries
	LOG_INFO << "Initialize control volume";
	c_neighborhoodSearch->set_active(false);
	for (unsigned int i = 0; i < numberOfControlBoundaryModels(); i++)
	{
		//if(!getControlBoundaryModel(i)->getRigidBodyObject()->isDynamic())
		c_neighborhoodSearch->set_active(i + nFluids, true, true);
	}
	c_neighborhoodSearch->find_neighbors();

	//if (if_target&&t_neighborhoodSearch != nullptr)
	//{
	//	t_neighborhoodSearch->set_active(false);
	//	for (unsigned int i = 0; i < numberOfControlBoundaryModels(); i++)
	//	{
	//		t_neighborhoodSearch->set_active(i, true, true);
	//	}
	//	LOG_INFO << "before fnid t_n";
	//	t_neighborhoodSearch->find_neighbors();
	//	LOG_INFO << "after find t_n";
	//}
	//Boundary objects
	for (unsigned int body = 0; body < numberOfControlBoundaryModels(); body++)
	{
		//if (!getControlBoundaryModel(body)->getRigidBodyObject()->isDynamic())
		ControlBoundaryModel* model = getControlBoundaryModel(body);
		model->computeControlBoundaryVolume();
		//表面
		const unsigned int numParticles = model->numActiveParticles();
		const unsigned int nControls = sim->numberOfControlBoundaryModels();

		//#pragma omp parallel default(shared)
		{
			//#pragma omp for schedule(static)
			for (int i = 0; i < (int)numParticles; i++)
			{
				Real rho = 0.0;
				Real& Ci = model->getC(i);
				const Vector3r& xi = model->getPosition(i);
				Vector3r& _x = model->getPosition_(i);
				int& edge = model->getEdge(i);
				int& face = model->getFace(i);
				//const Real& Vi = model->getVolume(i);
				int numberOfNei = 0;
				for (unsigned int pid = nFluids; pid < nFluids + nControls; pid++)
				{
					ControlBoundaryModel* cbm_neighbor = sim->getControlBoundaryModel(pid - nFluids);
					for (unsigned int j = 0; j < sim->cnumberOfNeighbors(body + nFluids, pid, i); j++)
					{
						const unsigned int neighborIndex = sim->cgetNeighbor(body + nFluids, pid, i, j);
						const Vector3r& xj = cbm_neighbor->getPosition(neighborIndex);
						if ((xi - xj).norm() <= c_supportRadius / 2.2)
						{
							rho += sim->W(xi - xj);
							numberOfNei++;
						}
						
					}
				}

				//Ci = std::max(0.66553 / rho,1.0);//small bunny 0.025
				//Ci = std::max(0.0765218 / rho, 1.0);//big bunny 0.05
				//Ci = std::max(0.775934 / rho,1.0);//Dragon50k -r 0.025 -s 2 alpha=0.002 beta= 0.04 gamma=0.085
				//Ci = std::max(0.77593 / rho, 1.0);//cat -r 0.025 alpha=0.002 beta= 0.04 gamma=0.085
				//Ci = std::max(0.717324 / rho, 1.0);//dog 0.0025 0.045 0.06
				//Ci = (xi - _x).norm() - 0.025;
				for (unsigned int pid = nFluids; pid < nFluids + nControls; pid++)
				{
					ControlBoundaryModel* cbm_neighbor = sim->getControlBoundaryModel(pid - nFluids);
					for (unsigned int j = 0; j < sim->cnumberOfNeighbors(body + nFluids, pid, i); j++)
					{
						const unsigned int neighborIndex = sim->cgetNeighbor(body + nFluids, pid, i, j);
						const Vector3r& xj = cbm_neighbor->getPosition(neighborIndex);
						if ((xi - xj).norm() <= c_supportRadius / 2.2)
						{
							Ci = sim->W(xi - xj) / (rho + 1.0e-6);
							_x += Ci * xj;
						}
					}
				}

				if (numberOfNei <= 40)
				{
					Ci = -0.025;//-0.85*r=-0.85*0.025=0.02125
					//edge = 1;
					face = 0;
				}
				else
				{
					Ci = (xi - _x).norm() - 0.025;
					//edge = 0;
					face = -1;
				}
					
				
			}
		}

		int inter = 0;
		bool chk = true;
		Real minC = -0.025;
		while (chk)
		{
//#pragma omp parallel default(shared)
			{
//#pragma omp for schedule(static)
				for (int i = 0; i < (int)numParticles; i++)//while's for
				{
					const Real& Ci = model->getC(i);
					const Vector3r& xi = model->getPosition(i);
					//const int& edge = model->getEdge(i);
					const int& face = model->getFace(i);
					if (face==inter)
					{
						for (unsigned int pid = nFluids; pid < nFluids + nControls; pid++)
						{
							ControlBoundaryModel* cbm_neighbor = sim->getControlBoundaryModel(pid - nFluids);

							for (unsigned int j = 0; j < sim->cnumberOfNeighbors(body + nFluids, pid, i); j++)
							{
								const unsigned int neighborIndex = sim->cgetNeighbor(body + nFluids, pid, i, j);
								const Vector3r& xj = cbm_neighbor->getPosition(neighborIndex);
								if ((xi - xj).norm() <= c_supportRadius / 2.2)
								{
									//int& edge_nei = cbm_neighbor->getEdge(neighborIndex);
									int& face_nei = cbm_neighbor->getFace(neighborIndex);
									Real& Cij = cbm_neighbor->getC(neighborIndex);
									if (face_nei == -1)
									{
										Cij = Ci - (xi - xj).norm();// std::max(Cij - (xi - xj).norm(), -0.45);//0.45=-2*9*r
										//edge_nei = inter+1;
										face_nei = inter + 1;
									}
								}
							}
						}
						
						if (Ci < minC)
							minC = Ci;
						//LOG_INFO << Ci;
					}
				}//while's for
				int numFace = 0;
				for (int i = 0; i < (int)numParticles; i++)
				{
					bool tmp = false;
					//const int& edge = model->getEdge(i);
					const int& face = model->getFace(i);
					if (face!=-1)
						numFace++;
				}
				//LOG_INFO << numEdge;
				if (numParticles-numFace==0)
					chk = false;
					
			}
			LOG_INFO << inter;
			inter++;
		}
		int f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10;
		f0 = f1 = f2 = f3 = f4 = f5 = f6 = f7 = f8 = f9 = f10 = 0;
		for (int i = 0; i < (int)numParticles; i++)
		{
			Real& Ci = model->getC(i);
			Ci = 1.0;
			const int face = model->getFace(i);
			switch (face)
			{
			case 0:
				f0++;
				//Ci +=1;
				//Ci = Ci / (minC + 1);
				Ci = 1.8;
				break;
			case 1:
				f1++;
				//Ci += 1;
				//Ci = Ci / (minC + 1);
				Ci = 1.4;
				break;
			case 2:
				f2++;
				//Ci = 1.6;
				break;
			case 3:
				f3++;
				//Ci = 1.4;
				break;
			case 4:
				f4++;
				//Ci = 1.2;
				break;
			case 5:
				f5++;
				break;
			case 6:
				f6++;
				break;
			case 7:
				f7++;
				break;
			case 8:
				f8++;
				break;
			case 9:
				f9++;
				break;
			case 10:
				f10++;
				break;
			default:
				break;
			}
		}
		LOG_INFO << "=============";
		LOG_INFO << f0;
		LOG_INFO << f1;
		LOG_INFO << f2;
		LOG_INFO << f3;
		LOG_INFO << f4;
		LOG_INFO << f5;
		LOG_INFO << f6;
		LOG_INFO << f7;
		LOG_INFO << f8;
		LOG_INFO << f9;
		LOG_INFO << f10;
		LOG_INFO << "----------------";
		LOG_INFO << f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8 + f9 + f10;
	}

	//todo 2020-10; zj
	////////////////////////////////////////////////////////////////////////// 
	// Compute for dynamic bodies
	//////////////////////////////////////////////////////////////////////////
	//ActivateOnlyFluids();
}

void Simulation::updateBoundaryVolume()
{
	if (m_neighborhoodSearch == nullptr)
		return;

	Simulation *sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();

	//////////////////////////////////////////////////////////////////////////
	// Compute value psi for boundary particles (boundary handling)
	// (see Akinci et al. "Versatile rigid - fluid coupling for incompressible SPH", Siggraph 2012
	//////////////////////////////////////////////////////////////////////////

	// Search boundary neighborhood

	// Activate only static boundaries
	LOG_INFO << "Initialize boundary volume";
	m_neighborhoodSearch->set_active(false);
	for (unsigned int i = 0; i < numberOfBoundaryModels(); i++)
	{
		if (!getBoundaryModel(i)->getRigidBodyObject()->isDynamic())
			m_neighborhoodSearch->set_active(i + nFluids, true, true);
	}

	m_neighborhoodSearch->find_neighbors();
	
	// Boundary objects
	for (unsigned int body = 0; body < numberOfBoundaryModels(); body++)
	{
		if (!getBoundaryModel(body)->getRigidBodyObject()->isDynamic())
			getBoundaryModel(body)->computeBoundaryVolume();
	}
	////////////////////////////////////////////////////////////////////////// 
	// Compute boundary psi for all dynamic bodies
	//////////////////////////////////////////////////////////////////////////
	for (unsigned int body = 0; body < numberOfBoundaryModels(); body++)
	{
		// Deactivate all
		m_neighborhoodSearch->set_active(false);

		// Only activate next dynamic body
		if (getBoundaryModel(body)->getRigidBodyObject()->isDynamic())
		{
			m_neighborhoodSearch->set_active(body + nFluids, true, true);
			m_neighborhoodSearch->find_neighbors();
			getBoundaryModel(body)->computeBoundaryVolume();
		}
	}
	//ActivateOnlyFluids();
}
void Simulation::ActivateNeighborFun()
{
	// Activate only fluids 
	m_neighborhoodSearch->set_active(false);
	if(ifControl())
		c_neighborhoodSearch->set_active(false);
	const unsigned int nFluids = numberOfFluidModels();
	for (unsigned int i = 0; i < numberOfFluidModels(); i++)
	{
		for (unsigned int j = 0; j < nFluids; j++)
			m_neighborhoodSearch->set_active(i, j, true);
		for (unsigned int j = nFluids; j < m_neighborhoodSearch->point_sets().size(); j++)
			m_neighborhoodSearch->set_active(i, j, true);
		if(ifControl())
			for (unsigned int j = nFluids; j < c_neighborhoodSearch->point_sets().size(); j++)
				c_neighborhoodSearch->set_active(i, j, true);			
	}
	if(ifControl())
		for (unsigned int i = 0; i < numberOfControlBoundaryModels(); i++)
		{
			for (unsigned int j = 0; j < numberOfFluidModels(); j++)
				c_neighborhoodSearch->set_active(i+nFluids, j, true);
		}
}
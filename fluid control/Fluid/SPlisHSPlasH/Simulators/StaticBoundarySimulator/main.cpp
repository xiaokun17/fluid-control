#include "SPlisHSPlasH/Common.h"
#include "GL/glew.h"
#include "Visualization/MiniGL.h"
#include "GL/glut.h"
#include "SPlisHSPlasH/TimeManager.h"
#include <Eigen/Dense>
#include <iostream>
#include "Utilities/Timing.h"
#include "Utilities/PartioReaderWriter.h"
#include "SPlisHSPlasH/StaticRigidBody.h"
#include "SPlisHSPlasH/ControlBody.h"//2020-10-11    zj
#include "Utilities/OBJLoader.h"
#include "SPlisHSPlasH/Utilities/PoissonDiskSampling.h"
#include "Simulators/Common/SimulatorBase.h"
#include "Utilities/FileSystem.h"
#include "Utilities/Version.h"
#include <fstream>
#include "Utilities/Logger.h"
#include "Utilities/Counting.h"
#include "SPlisHSPlasH/Simulation.h"
#include "Simulators/Common/TweakBarParameters.h"
#include "GL/freeglut_ext.h"

// Enable memory leak detection
#ifdef _DEBUG
#ifndef EIGEN_ALIGN
	#define new DEBUG_NEW 
#endif
#endif

using namespace SPH;
using namespace Eigen;
using namespace std;
using namespace Utilities;
using namespace GenParam;

void timeStep ();
void initBoundaryData();
void initControlData();
void initControlModel();
void initTargetModel();
void render ();
void renderBoundary();
void reset();
void initParameters();
void loadObj(const std::string &filename, TriangleMesh &geo, const Vector3r &scale);
void TW_CALL setCurrentFluidModel(const void *value, void *clientData);
void TW_CALL getCurrentFluidModel(void *value, void *clientData);
void TW_CALL setColorField(const void *value, void *clientData);
void TW_CALL getColorField(void *value, void *clientData);
void TW_CALL setRenderMaxValue(const void *value, void *clientData);
void TW_CALL getRenderMaxValue(void *value, void *clientData);
void TW_CALL setRenderMinValue(const void *value, void *clientData);
void TW_CALL getRenderMinValue(void *value, void *clientData);
void TW_CALL setColorMapType(const void *value, void *clientData);
void TW_CALL getColorMapType(void *value, void *clientData);

SimulatorBase *base;
unsigned int currentFluidModel = 0;
std::vector<std::string> colorFieldNames;

// main 
int main( int argc, char **argv )
{
	REPORT_MEMORY_LEAKS;

	base = new SimulatorBase();
	base->init(argc, argv, "StaticBoundarySimulator");

	Simulation *sim = Simulation::getCurrent();
	sim->init(base->getScene().particleRadius, base->getScene().sim2D);
	
	base->buildModel();

	initBoundaryData();
	
	unsigned int nBoundaryParticles = 0;
	for (unsigned int i = 0; i < sim->numberOfBoundaryModels(); i++)
		nBoundaryParticles += sim->getBoundaryModel(i)->numberOfParticles();

	LOG_INFO << "Number of boundary particles: " << nBoundaryParticles;
	//initControlData();
	if (sim->ifControl())
	{
		initControlModel();
		unsigned int nControlParticles = 0;
		for (unsigned int i = 0; i < sim->numberOfControlBoundaryModels(); i++)
			nControlParticles += sim->getControlBoundaryModel(i)->numActiveParticles();
		LOG_INFO << "Number of control particles : " << nControlParticles;
		//if (sim->ifTarget())
		//{
		//	initTargetModel();
		//}
	}
	sim->ActivateNeighborFun();
	for (unsigned int i = 0; i < sim->numberOfFluidModels(); i++)
	{
		FluidModel *model = sim->getFluidModel(i);
		const std::string &key = model->getId();
		model->setDragMethodChangedCallback([&]() { initParameters(); base->getSceneLoader()->readParameterObject(key, (ParameterObject*)model->getDragBase()); });
		model->setSurfaceMethodChangedCallback([&]() { initParameters(); base->getSceneLoader()->readParameterObject(key, (ParameterObject*)model->getSurfaceTensionBase()); });
		model->setViscosityMethodChangedCallback([&]() { initParameters(); base->getSceneLoader()->readParameterObject(key, (ParameterObject*)model->getViscosityBase()); });
		model->setVorticityMethodChangedCallback([&]() { initParameters(); base->getSceneLoader()->readParameterObject(key, (ParameterObject*)model->getVorticityBase()); });
		model->setElasticityMethodChangedCallback([&]() { reset(); initParameters(); base->getSceneLoader()->readParameterObject(key, (ParameterObject*)model->getElasticityBase()); });
	}

	initParameters();
	
	Simulation::getCurrent()->setSimulationMethodChangedCallback([&]() { reset(); initParameters(); base->getSceneLoader()->readParameterObject("Configuration", Simulation::getCurrent()->getTimeStep()); });
	base->readParameters();
	
	MiniGL::setClientIdleFunc(50, timeStep);
	MiniGL::addKeyFunc('r', reset);
	MiniGL::setClientSceneFunc(render);


	glutMainLoop ();	

	base->cleanup ();

	Utilities::Timing::printAverageTimes();
	Utilities::Timing::printTimeSums();

	Utilities::Counting::printAverageCounts();
	Utilities::Counting::printCounterSums();

	delete Simulation::getCurrent();
	delete base;
	
	return 0;
}

void initParameters()
{
	TwRemoveAllVars(MiniGL::getTweakBar());
	TweakBarParameters::cleanup();

	MiniGL::initTweakBarParameters();

	Simulation *sim = Simulation::getCurrent();
	TweakBarParameters::createParameterGUI();
	TweakBarParameters::createParameterObjectGUI(base);
	TweakBarParameters::createParameterObjectGUI(sim);
	TweakBarParameters::createParameterObjectGUI(sim->getTimeStep());
	// Enum for all fluid models
	if (sim->numberOfFluidModels() > 0)
	{
		TwType enumType = TwDefineEnum("CurrentFluidModelEnum", NULL, 0);
		std::ostringstream oss;
		oss << 0 << " {" << sim->getFluidModel(0)->getId().c_str() << "}";
		for (unsigned int j = 1; j < sim->numberOfFluidModels(); j++)
		{
			oss << ", " << j << " {" << sim->getFluidModel(j)->getId().c_str() << "}";
		}
		std::string enumStr = " label='Current fluid model' enum='" + oss.str() + "' group='Fluid model'";
		TwAddVarCB(MiniGL::getTweakBar(), "CurrentFluidModel", enumType, setCurrentFluidModel, getCurrentFluidModel, &currentFluidModel, enumStr.c_str());
	}
	// show GUI only for currently selected fluid model
	unsigned int i = currentFluidModel;
	FluidModel *model = sim->getFluidModel(currentFluidModel);
	colorFieldNames.clear();
	colorFieldNames.resize(model->numberOfFields());
	TwType enumType = TwDefineEnum("ColorFieldEnum", NULL, 0);
	std::ostringstream oss;
	for (unsigned int j = 0; j < model->numberOfFields(); j++)
	{
		const FieldDescription &field = model->getField(j);
		if ((field.type == FieldType::Scalar) || (field.type == FieldType::Vector3))
		{
			if (j != 0)
				oss << ", ";
			oss << j << " {" << field.name.c_str() << "}";
			colorFieldNames[j] = field.name;
		}
	}
	std::string enumStr = " label='Color field' enum='" + oss.str() + "' group='" + model->getId() + "'";
	enumStr = enumStr + " help='Choose vector or scalar field for particle coloring.'";
	TwAddVarCB(MiniGL::getTweakBar(), "ColorField", enumType, setColorField, getColorField, nullptr, enumStr.c_str());
	
	TwType enumType2 = TwDefineEnum("ColorMapTypeEnum", NULL, 0);
	std::string str = " label='Color map' enum='0 {None}, 1 {Jet}, 2 {Plasma}' group='" + model->getId() + "' help='Choose a color map.'";
	TwAddVarCB(MiniGL::getTweakBar(), "ColorMapType", enumType2, setColorMapType, getColorMapType, nullptr, str.c_str());
	str = " label='Min. value (shader)' step=0.001 precision=3 group='" + model->getId() + "' help='Minimal value used for color-coding the color field in the rendering process.'";
	TwAddVarCB(MiniGL::getTweakBar(), "RenderMinValue", TW_TYPE_REAL, setRenderMinValue, getRenderMinValue, nullptr, str.c_str());
	str = " label='Max. value (shader)' step=0.001 precision=3 group='" + model->getId() + "' help='Maximal value used for color-coding the color field in the rendering process.'";
	TwAddVarCB(MiniGL::getTweakBar(), "RenderMaxValue", TW_TYPE_REAL, setRenderMaxValue, getRenderMaxValue, nullptr, str.c_str());

	TweakBarParameters::createParameterObjectGUI(model);
	TweakBarParameters::createParameterObjectGUI((GenParam::ParameterObject*) model->getDragBase());
	TweakBarParameters::createParameterObjectGUI((GenParam::ParameterObject*) model->getSurfaceTensionBase());
	TweakBarParameters::createParameterObjectGUI((GenParam::ParameterObject*) model->getViscosityBase());
	TweakBarParameters::createParameterObjectGUI((GenParam::ParameterObject*) model->getVorticityBase());
	TweakBarParameters::createParameterObjectGUI((GenParam::ParameterObject*) model->getElasticityBase());
	TwDefine((std::string("TweakBar/FluidModel group='") + model->getId() + "'").c_str());
	TwDefine((std::string("TweakBar/'Drag force' group='") + model->getId() + "'").c_str());
	TwDefine((std::string("TweakBar/'Surface tension' group='") + model->getId() + "'").c_str());
	TwDefine((std::string("TweakBar/Viscosity group='") + model->getId() + "'").c_str());
	TwDefine((std::string("TweakBar/Vorticity group='") + model->getId() + "'").c_str());
	TwDefine((std::string("TweakBar/'Elasticity' group='") + model->getId() + "'").c_str());
}

void reset()
{
	Utilities::Timing::printAverageTimes();
	Utilities::Timing::reset();

	Utilities::Counting::printAverageCounts();
	Utilities::Counting::reset();

	Simulation::getCurrent()->reset();
	base->reset();
	base->getSelectedParticles().clear();
}

void timeStep ()
{
	const Real stopAt = base->getValue<Real>(SimulatorBase::STOP_AT);
	if ((stopAt > 0.0) && (stopAt < TimeManager::getCurrent()->getTime()))
		glutLeaveMainLoop();

	const Real pauseAt = base->getValue<Real>(SimulatorBase::PAUSE_AT);
	if ((pauseAt > 0.0) && (pauseAt < TimeManager::getCurrent()->getTime()))
		base->setValue(SimulatorBase::PAUSE, true);

	if (base->getValue<bool>(SimulatorBase::PAUSE))
		return;

	// Simulation code
	Simulation *sim = Simulation::getCurrent();
	const bool sim2D = sim->is2DSimulation();
	const unsigned int numSteps = base->getValue<unsigned int>(SimulatorBase::NUM_STEPS_PER_RENDER);
	for (unsigned int i = 0; i < numSteps; i++)
	{
		START_TIMING("SimStep");
		Simulation::getCurrent()->getTimeStep()->step();
		STOP_TIMING_AVG;

		base->step();

		INCREASE_COUNTER("Time step size", TimeManager::getCurrent()->getTimeStepSize());

		// Make sure that particles stay in xy-plane in a 2D simulation
		if (sim2D)
		{
			for (unsigned int i = 0; i < sim->numberOfFluidModels(); i++)
			{
				FluidModel *model = sim->getFluidModel(i);
				for (unsigned int i = 0; i < model->numActiveParticles(); i++)
				{
					model->getPosition(i)[2] = 0.0;
					model->getVelocity(i)[2] = 0.0;
				}
			}
		}
	}
}

void render()
{
	float gridColor[4] = { 0.2f, 0.2f, 0.2f, 1.0f };
	const bool sim2D = Simulation::getCurrent()->is2DSimulation();
	if (sim2D)
		MiniGL::drawGrid_xy(gridColor);
	else
		MiniGL::drawGrid_xz(gridColor);

	MiniGL::coordinateSystem();
	MiniGL::drawTime(TimeManager::getCurrent()->getTime());

	Simulation *sim = Simulation::getCurrent();
	for (unsigned int i = 0; i < sim->numberOfFluidModels(); i++)
	{
		float fluidColor[4] = { 0.3f, 0.5f, 0.9f, 1.0f };
		MiniGL::hsvToRgb(0.61f-0.1f*i, 0.66f, 0.9f, fluidColor);
		base->renderFluid(i, fluidColor);
	}
	renderBoundary();
}

void renderBoundary()
{
	Simulation* sim = Simulation::getCurrent();
	Shader& shader = base->getShaderScalar();
	Shader& meshShader = base->getMeshShader();
	SceneLoader::Scene& scene = base->getScene();
	const int renderWalls = base->getValue<int>(SimulatorBase::RENDER_WALLS);
	GLint context_major_version = base->getContextMajorVersion();

	if ((renderWalls == 1) || (renderWalls == 2))
	{
		if (context_major_version > 3)
		{
			shader.begin();
			for (int body = sim->numberOfBoundaryModels() - 1; body >= 0; body--)
			{
				if ((renderWalls == 1) || (!scene.boundaryModels[body]->isWall))
				{
					glUniform3fv(shader.getUniform("color"), 1, scene.boundaryModels[body]->color.data());
					glEnableVertexAttribArray(0);

					BoundaryModel* bm = sim->getBoundaryModel(body);
					glVertexAttribPointer(0, 3, GL_REAL, GL_FALSE, 0, &bm->getPosition(0));
					glDrawArrays(GL_POINTS, 0, bm->numberOfParticles());
					glDisableVertexAttribArray(0);
				}
			}
			shader.end();
		}
		else
		{
			glDisable(GL_LIGHTING);
			glPointSize(4.0f);

			glBegin(GL_POINTS);
			for (int body = sim->numberOfBoundaryModels() - 1; body >= 0; body--)
			{
				if ((renderWalls == 1) || (!scene.boundaryModels[body]->isWall))
				{
					BoundaryModel* bm = sim->getBoundaryModel(body);
					for (unsigned int i = 0; i < bm->numberOfParticles(); i++)
					{
						glColor3fv(scene.boundaryModels[body]->color.data());
						glVertex3v(&bm->getPosition(i)[0]);
					}
				}
			}
			glEnd();
			glEnable(GL_LIGHTING);
		}
	}
	else if ((renderWalls == 3) || (renderWalls == 4))
	{
		for (int body = sim->numberOfBoundaryModels() - 1; body >= 0; body--)
		{
			if ((renderWalls == 3) || (!scene.boundaryModels[body]->isWall))
			{
				meshShader.begin();
				glUniform1f(meshShader.getUniform("shininess"), 5.0f);
				glUniform1f(meshShader.getUniform("specular_factor"), 0.2f);

				GLfloat matrix[16];
				glGetFloatv(GL_MODELVIEW_MATRIX, matrix);
				glUniformMatrix4fv(meshShader.getUniform("modelview_matrix"), 1, GL_FALSE, matrix);
				GLfloat pmatrix[16];
				glGetFloatv(GL_PROJECTION_MATRIX, pmatrix);
				glUniformMatrix4fv(meshShader.getUniform("projection_matrix"), 1, GL_FALSE, pmatrix);

				glUniform3fv(meshShader.getUniform("surface_color"), 1, scene.boundaryModels[body]->color.data());

				BoundaryModel* bm = sim->getBoundaryModel(body);
				MiniGL::drawMesh(((StaticRigidBody*)bm->getRigidBodyObject())->getGeometry(), scene.boundaryModels[body]->color.data());

				meshShader.end();
			}
		}
	}
	else if (renderWalls == 5)
	{
		for (int body = sim->numberOfBoundaryModels() - 1; body >= 0; body--)
		{
			if (!scene.boundaryModels[body]->isRig && !scene.boundaryModels[body]->isWall)
			{
				meshShader.begin();
				glUniform1f(meshShader.getUniform("shininess"), 5.0f);
				glUniform1f(meshShader.getUniform("specular_factor"), 0.2f);

				GLfloat matrix[16];
				glGetFloatv(GL_MODELVIEW_MATRIX, matrix);
				glUniformMatrix4fv(meshShader.getUniform("modelview_matrix"), 1, GL_FALSE, matrix);
				GLfloat pmatrix[16];
				glGetFloatv(GL_PROJECTION_MATRIX, pmatrix);
				glUniformMatrix4fv(meshShader.getUniform("projection_matrix"), 1, GL_FALSE, pmatrix);

				glUniform3fv(meshShader.getUniform("surface_color"), 1, scene.boundaryModels[body]->color.data());

				BoundaryModel* bm = sim->getBoundaryModel(body);
				MiniGL::drawMesh(((StaticRigidBody*)bm->getRigidBodyObject())->getGeometry(), scene.boundaryModels[body]->color.data());

				meshShader.end();
			}
		}
	}
	else if (renderWalls == 6)
	{
		glDisable(GL_LIGHTING);
		glPointSize(4.0f);
		glBegin(GL_POINTS);
		const float colors[11][3] = {
			{1,1,1},
			{0.8,0.2,0},
			{0.6,0.4,0},
			{0.4,0.6,0},
			{0.2,0.8,0},
			{0,1.0,0},
			{0,0.8,0.2},
			{0,0.6,0.4},
			{0,0.4,0.6},
			{0,0.2,0.8},
			{0,0,1.0},
		};
		for (int body = sim->numberOfControlBoundaryModels() - 1; body >= 0; body--)
		{
			if ((renderWalls  ==6 ))//|| (!scene.controlModels[body]->isWall))
			{
				ControlBoundaryModel* cbm = sim->getControlBoundaryModel(body);
				for (unsigned int i = 0; i < cbm->numberOfParticles(); i++)
				{
					const int face = cbm->getFace(i);
					glColor3fv(colors[face]);
					
					//glColor3fv(scene.controlModels[body]->color.data());
					glVertex3v(&cbm->getPosition(i)[0]);
				}
			}
		}
		glEnd();
		glEnable(GL_LIGHTING);
	}
	else if (renderWalls == 7)
	{
		glDisable(GL_LIGHTING);
		glPointSize(4.0f);

		glBegin(GL_POINTS);
		for (int body = sim->numberOfTargetModels() - 1; body >= 0; body--)
		{
			if ((renderWalls == 6))// || (!scene.controlModels[body]->isWall))
			{
				TargetModel* tm = sim->getTargetModel(body);
				for (unsigned int i = 0; i < tm->numberOfParticles(); i++)
				{
					glColor3fv(scene.controlModels[body]->color.data());
					glVertex3v(&tm->getPosition(i)[0]);
				}
			}
		}
		glEnd();
		glEnable(GL_LIGHTING);
	}
}

void loadObj(const std::string &filename, TriangleMesh &mesh, const Vector3r &scale)
{
	std::vector<OBJLoader::Vec3f> x;
	std::vector<OBJLoader::Vec3f> normals;
	std::vector<MeshFaceIndices> faces;
	OBJLoader::Vec3f s = { (float) scale[0], (float)scale[1], (float)scale[2] };
	OBJLoader::loadObj(filename, &x, &faces, &normals, nullptr, s);

	mesh.release();
	const unsigned int nPoints = (unsigned int)x.size();
	const unsigned int nFaces = (unsigned int)faces.size();
	mesh.initMesh(nPoints, nFaces);
	for (unsigned int i = 0; i < nPoints; i++)
	{
		mesh.addVertex(Vector3r(x[i][0], x[i][1], x[i][2]));
	}
	for (unsigned int i = 0; i < nFaces; i++)
	{
		// Reduce the indices by one
		int posIndices[3];
		for (int j = 0; j < 3; j++)
		{
			posIndices[j] = faces[i].posIndices[j] - 1;
		}

		mesh.addFace(&posIndices[0]);
	}

	LOG_INFO << "Number of triangles: " << nFaces;
	LOG_INFO << "Number of vertices: " << nPoints;
}
void initTargetModel()
{
	LOG_INFO << "Add target particles......from volume sampling(.bgeo)";
	Simulation* sim = Simulation::getCurrent();
	SceneLoader::Scene& m_scene = base->getScene();
	std::map<std::string, unsigned int> targetIDs;
	unsigned int index = 0;
	for (unsigned int i = 0; i < m_scene.targetModels.size(); i++)
	{
		if (targetIDs.find(m_scene.targetModels[i]->id) == targetIDs.end())
			targetIDs[m_scene.targetModels[i]->id] = index++;
	}

	const unsigned int numberOfTargetModels = static_cast<unsigned int>(targetIDs.size());
	std::vector<std::vector<Vector3r>> targetParticles;
	std::vector<std::vector<Vector3r>> targetVelocities;
	targetParticles.resize(numberOfTargetModels);
	targetVelocities.resize(numberOfTargetModels);

	std::string base_path = FileSystem::getFilePath(base->getSceneFile());
	for (unsigned int i = 0; i < m_scene.targetModels.size(); i++)
	{
		const unsigned int targetIndex = targetIDs[m_scene.targetModels[i]->id];

		std::string fileName;
		if (FileSystem::isRelativePath(m_scene.targetModels[i]->samplesFile))
			fileName = base_path + "/" + m_scene.targetModels[i]->samplesFile;
		else
			fileName = m_scene.targetModels[i]->samplesFile;

		PartioReaderWriter::readParticles(fileName, m_scene.targetModels[i]->translation, m_scene.targetModels[i]->rotation,m_scene.targetModels[i]->scale2, targetParticles[targetIndex],targetVelocities[targetIndex]);
	}
	for (auto it = targetIDs.begin(); it != targetIDs.end(); it++)
	{
		const unsigned int index = it->second;
		Simulation::getCurrent()->addTargetModel(it->first, (unsigned int)targetParticles[index].size(), targetParticles[index].data());
	}
	Simulation::getCurrent()->updateTargetSearch();
}
void initControlModel()
{
	LOG_INFO << "Initialize control particles......from volume sampling(.bgeo)";

	Simulation* sim = Simulation::getCurrent();
	SceneLoader::Scene& m_scene = base->getScene();
	//////////////////////////////////////////////////////////////////////////
	// Determine number of different control IDs
	//////////////////////////////////////////////////////////////////////////
	std::map<std::string, unsigned int> controlIDs;
	unsigned int index = 0;
	for (unsigned int i = 0; i < m_scene.controlModels.size(); i++)
	{
		if (controlIDs.find(m_scene.controlModels[i]->id) == controlIDs.end())
			controlIDs[m_scene.controlModels[i]->id] = index++;
	}
	const unsigned int numberOfControlModels = static_cast<unsigned int>(controlIDs.size());
	
	std::vector<std::vector<Vector3r>> controlParticles;
	std::vector<std::vector<Vector3r>> controlVelocities;
	controlParticles.resize(numberOfControlModels);
	controlVelocities.resize(numberOfControlModels);

	std::string base_path = FileSystem::getFilePath(base->getSceneFile());

	for (unsigned int i = 0; i < m_scene.controlModels.size(); i++)
	{
		const unsigned int controlIndex = controlIDs[m_scene.controlModels[i]->id];

		std::string fileName;
		if (FileSystem::isRelativePath(m_scene.controlModels[i]->samplesFile))
			fileName = base_path + "/" + m_scene.controlModels[i]->samplesFile;
		else
			fileName = m_scene.controlModels[i]->samplesFile;

		PartioReaderWriter::readParticles(fileName, m_scene.controlModels[i]->translation, m_scene.controlModels[i]->rotation, m_scene.controlModels[i]->scale2, controlParticles[controlIndex], controlVelocities[controlIndex]);
		//Simulation::getCurrent()->setValue(Simulation::PARTICLE_RADIUS, m_scene.particleRadius);
	}
	for (auto it = controlIDs.begin(); it != controlIDs.end(); it++)
	{
		const unsigned int index = it->second;
		Simulation::getCurrent()->addControlModelf(it->first, (unsigned int)controlParticles[index].size(),controlParticles[index].data());
	}
	Simulation::getCurrent()->updateControlVloume();
}
void initControlData()
{
	std::string scene_path = FileSystem::getFilePath(base->getSceneFile());
	std::string scene_file_name = FileSystem::getFileName(base->getSceneFile());
	SceneLoader::Scene& scene = base->getScene();
	const bool useCache = base->getUseParticleCaching();

	string cachePath = scene_path + "/ControlCache";
	
	for (unsigned int i = 0; i < scene.controlModels.size(); i++)
	{
		string meshFileName = FileSystem::normalizePath(scene_path + "/" + scene.controlModels[i]->meshFile);

		std::string md5FileName = FileSystem::normalizePath(cachePath + "/" + FileSystem::getFileNameWithExt(meshFileName) + ".md5");
		bool md5 = false;
		if (useCache)
		{
			string md5Str = FileSystem::getFileMD5(meshFileName);
			if (FileSystem::fileExists(md5FileName))
				md5 = FileSystem::checkMD5(md5Str, md5FileName);
		}

		std::vector<Vector3r> controlParticles;
		if (scene.controlModels[i]->samplesFile != "") //�����ļ������Բ���
		{
			string particleFileName = scene_path + "/" + scene.controlModels[i]->samplesFile;
			PartioReaderWriter::readParticles(particleFileName, scene.controlModels[i]->translation, scene.controlModels[i]->rotation, scene.controlModels[i]->scale[0], controlParticles);
		}
		ControlBody* cb = new ControlBody();
		TriangleMesh& geo = cb->getGeometry();
		loadObj(meshFileName, geo, scene.controlModels[i]->scale);

		if (scene.controlModels[i]->samplesFile == "")
		{
			// Cache sampling
			std::string mesh_base_path = FileSystem::getFilePath(scene.controlModels[i]->meshFile);
			std::string mesh_file_name = FileSystem::getFileName(scene.controlModels[i]->meshFile);

			const string resStr = to_string(scene.controlModels[i]->scale[0]) + "_" + to_string(scene.controlModels[i]->scale[1]) + "_" + to_string(scene.controlModels[i]->scale[2]);
			const string particleFileName = FileSystem::normalizePath(cachePath + "/" + mesh_file_name + "_sbd_" + std::to_string(scene.particleRadius) + "_" + resStr + ".bgeo");

			// check MD5 if cache file is available
			bool foundCacheFile = false;

			if (useCache)
				foundCacheFile = FileSystem::fileExists(particleFileName);

			if (useCache && foundCacheFile && md5)
			{
				PartioReaderWriter::readParticles(particleFileName, scene.controlModels[i]->translation, scene.controlModels[i]->rotation, 1.0, controlParticles);
				//(particleFileName, scene.boundaryModels[i]->translation, scene.boundaryModels[i]->rotation, 1.0, boundaryParticles);
				LOG_INFO << "Loaded cached control sampling: " << particleFileName;
			}

			if (!useCache || !foundCacheFile || !md5)
			{
				LOG_INFO << "Surface sampling of " << meshFileName;
				START_TIMING("Poisson disk sampling");
				PoissonDiskSampling sampling;
				sampling.sampleMesh(geo.numVertices(), geo.getVertices().data(), geo.numFaces(), geo.getFaces().data(), scene.particleRadius, 10, 1, controlParticles);
				STOP_TIMING_AVG;

				// Cache sampling
				if (useCache && (FileSystem::makeDir(cachePath) == 0))
				{
					LOG_INFO << "Save particle sampling: " << particleFileName;
					PartioReaderWriter::writeParticles(particleFileName, (unsigned int)controlParticles.size(), controlParticles.data(), nullptr, scene.particleRadius);
					FileSystem::writeMD5File(meshFileName, md5FileName);
				}

				// transform particles
				for (unsigned int j = 0; j < (unsigned int)controlParticles.size(); j++)
					controlParticles[j] = scene.controlModels[i]->rotation * controlParticles[j] + scene.controlModels[i]->translation;
			}
		}
		for (unsigned int j = 0; j < geo.numVertices(); j++)
			geo.getVertices()[j] = scene.controlModels[i]->rotation * geo.getVertices()[j] + scene.controlModels[i]->translation;

		geo.updateNormals();
		geo.updateVertexNormals();

		Simulation::getCurrent()->addControlModel(cb, i, static_cast<unsigned int>(controlParticles.size()), &controlParticles[0]);
	}
	Simulation::getCurrent()->updateControlVloume();
}
void initBoundaryData()
{
	std::string scene_path = FileSystem::getFilePath(base->getSceneFile());
	std::string scene_file_name = FileSystem::getFileName(base->getSceneFile());
	SceneLoader::Scene &scene = base->getScene();
	const bool useCache = base->getUseParticleCaching();

	string cachePath = scene_path + "/Cache";

	for (unsigned int i = 0; i < scene.boundaryModels.size(); i++)
	{
		string meshFileName = FileSystem::normalizePath(scene_path + "/" + scene.boundaryModels[i]->meshFile);

		std::string md5FileName = FileSystem::normalizePath(cachePath + "/" + FileSystem::getFileNameWithExt(meshFileName) + ".md5");
		bool md5 = false;
		if (useCache)
		{
			string md5Str = FileSystem::getFileMD5(meshFileName);
			if (FileSystem::fileExists(md5FileName))
				md5 = FileSystem::checkMD5(md5Str, md5FileName);
		}

		std::vector<Vector3r> boundaryParticles;
		if (scene.boundaryModels[i]->samplesFile != "")
		{
			string particleFileName = scene_path + "/" + scene.boundaryModels[i]->samplesFile;
			PartioReaderWriter::readParticles(particleFileName, scene.boundaryModels[i]->translation, scene.boundaryModels[i]->rotation, scene.boundaryModels[i]->scale[0], boundaryParticles);
		}

		StaticRigidBody *rb = new StaticRigidBody();
		TriangleMesh &geo = rb->getGeometry();
		loadObj(meshFileName, geo, scene.boundaryModels[i]->scale);

		if (scene.boundaryModels[i]->samplesFile == "")
		{
			// Cache sampling
			std::string mesh_base_path = FileSystem::getFilePath(scene.boundaryModels[i]->meshFile);
			std::string mesh_file_name = FileSystem::getFileName(scene.boundaryModels[i]->meshFile);

			const string resStr = to_string(scene.boundaryModels[i]->scale[0]) + "_" + to_string(scene.boundaryModels[i]->scale[1]) + "_" + to_string(scene.boundaryModels[i]->scale[2]);
			const string particleFileName = FileSystem::normalizePath(cachePath + "/" + mesh_file_name + "_sbd_" + std::to_string(scene.particleRadius) + "_" + resStr + ".bgeo");

			// check MD5 if cache file is available
			bool foundCacheFile = false;

			if (useCache)
				foundCacheFile = FileSystem::fileExists(particleFileName);

			if (useCache && foundCacheFile && md5)
			{
				PartioReaderWriter::readParticles(particleFileName, scene.boundaryModels[i]->translation, scene.boundaryModels[i]->rotation, 1.0, boundaryParticles);
				LOG_INFO << "Loaded cached boundary sampling: " << particleFileName;
			}

			if (!useCache || !foundCacheFile || !md5)
			{
				LOG_INFO << "Surface sampling of " << meshFileName;
				START_TIMING("Poisson disk sampling");
				PoissonDiskSampling sampling;
				sampling.sampleMesh(geo.numVertices(), geo.getVertices().data(), geo.numFaces(), geo.getFaces().data(), scene.particleRadius, 10, 1, boundaryParticles);
				STOP_TIMING_AVG;

				// Cache sampling
				if (useCache && (FileSystem::makeDir(cachePath) == 0))
				{
					LOG_INFO << "Save particle sampling: " << particleFileName;
					PartioReaderWriter::writeParticles(particleFileName, (unsigned int)boundaryParticles.size(), boundaryParticles.data(), nullptr, scene.particleRadius);
					FileSystem::writeMD5File(meshFileName, md5FileName);
				}

				// transform particles
				for (unsigned int j = 0; j < (unsigned int)boundaryParticles.size(); j++)
					boundaryParticles[j] = scene.boundaryModels[i]->rotation * boundaryParticles[j] + scene.boundaryModels[i]->translation;
			}
		}
		for (unsigned int j = 0; j < geo.numVertices(); j++)
			geo.getVertices()[j] = scene.boundaryModels[i]->rotation * geo.getVertices()[j] + scene.boundaryModels[i]->translation;

		geo.updateNormals();
		geo.updateVertexNormals();

		Simulation::getCurrent()->addBoundaryModel(rb, static_cast<unsigned int>(boundaryParticles.size()), &boundaryParticles[0]);
	}
	Simulation::getCurrent()->updateBoundaryVolume();
	

#ifdef GPU_NEIGHBORHOOD_SEARCH
	// copy the particle data to the GPU
	Simulation::getCurrent()->getNeighborhoodSearch()->update_point_sets();
#endif 
}

void TW_CALL setCurrentFluidModel(const void *value, void *clientData)
{
	const unsigned int val = *(const unsigned int *)(value);
	*((unsigned int*)clientData) = val;

	initParameters();
}

void TW_CALL getCurrentFluidModel(void *value, void *clientData)
{
	*(unsigned int *)(value) = *((unsigned int*)clientData);
}

void TW_CALL setColorField(const void *value, void *clientData)
{
	const unsigned int val = *(const unsigned int *)(value);
	base->setColorField(currentFluidModel, colorFieldNames[val]);
}

void TW_CALL getColorField(void *value, void *clientData)
{
	const std::string &fieldName = base->getColorField(currentFluidModel);
	unsigned int index = 0;
	for (auto i = 0; i < colorFieldNames.size(); i++)
	{
		if (colorFieldNames[i] == fieldName)
		{
			index = i;
			break;
		}
	}
	*(unsigned int *)(value) = index;
}

void TW_CALL setRenderMaxValue(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);
	base->setRenderMaxValue(currentFluidModel, val);
}

void TW_CALL getRenderMaxValue(void *value, void *clientData)
{
	*(Real *)(value) = base->getRenderMaxValue(currentFluidModel);
}

void TW_CALL setRenderMinValue(const void *value, void *clientData)
{
	const Real val = *(const Real *)(value);
	base->setRenderMinValue(currentFluidModel, val);
}

void TW_CALL getRenderMinValue(void *value, void *clientData)
{
	*(Real *)(value) = base->getRenderMinValue(currentFluidModel);
}

void TW_CALL setColorMapType(const void *value, void *clientData)
{
	const unsigned int val = *(const unsigned int *)(value);
	base->setColorMapType(currentFluidModel, val);
}

void TW_CALL getColorMapType(void *value, void *clientData)
{
	*(unsigned int *)(value) = base->getColorMapType(currentFluidModel);
}

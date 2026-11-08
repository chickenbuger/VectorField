#pragma once

#include "Vec3.h"
#include "Particle.h"
#include <vector>
#include "VectorField2D.h"
#include "GL/glut.h"

using namespace std;

#define G 0.06674 //6.674x10^-11

class GPlanet
{
public:
	double	_r;
	double	_mass;
	double _gravity;
	Vec3<double> _pos;
public:
	GPlanet();
	GPlanet(double x, double y);
	GPlanet(double x, double y, double z);
	GPlanet(Vec3<double> pos);
	~GPlanet();
public:
	void Draw();
};

class Planet
{
public:
	VectorField2D* _vF;
	vector<GPlanet*> _gP;
public:
	Planet();
	~Planet();
public:
	void Update();
	void CreatePlanet();
	void CreatePlanet(double x, double y);
	void CreatePlanet(Vec3<double> pos);
public:
	void Draw();
};


#ifndef __VECTOR_FIELD_2D_H__
#define	__VECTOR_FIELD_2D_H__

#pragma once

#include "Vec3.h"
#include "Particle.h"
#include <vector>
#include "GL\glut.h"

using namespace std;

class VectorField2D
{
public:
	vector<Particle*>	m_Particles; //모든 입자들
	double						m_diff;
public:
	VectorField2D();
	VectorField2D(int x, int y);
	~VectorField2D();
public:
	void Init(int x, int y);
	void Update();
	void	SettingForce();
public:
	void Draw();
	void PlanetForceDraw();
};

#endif
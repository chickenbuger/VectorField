#ifndef __VECTOR_FIELD_3D_H__
#define	__VECTOR_FIELD_3D_H__

#pragma once

#include "Vec3.h"
#include "Particle.h"
#include <vector>
#include "GL\glut.h"

using namespace std;

class VectorField3D
{
public:
	vector<Particle*> m_Particles; //모든 입자들
public:
	VectorField3D();
	VectorField3D(int x, int y, int z);
	~VectorField3D();
public:
	void Init(int x, int y, int z);
	void Update();
	void	SettingForce();
public:
	void Draw();
};

#endif
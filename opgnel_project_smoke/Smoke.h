#ifndef __SMOKE_H__
#define __SMOKE_H__

#pragma once

#include "Vec3.h"
#include "Particle.h"
#include <vector>
#include "GL/glut.h"

using namespace std;

class Smoke
{
public:
	vector<Particle*>	m_Particle; //���ڵ�
public:
	double						m_Viscosity; //����
	double						m_Pressure; //�з�
	double						m_Density; //�е�
public:
	Smoke(void) {};
	~Smoke(void) {};
public:
	void Start(void); //���� ���� ����
	void Update(void);
	void Init(void); //�ʱⰪ ����
public:
	void SetPressure(double newP);
public:
	void Draw(void);
};

#endif
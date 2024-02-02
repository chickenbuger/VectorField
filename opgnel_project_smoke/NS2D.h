#ifndef __NS_2D_H__
#define __NS_2D_H__

#pragma once

#include "Vec3.h"
#include "Particle.h"
#include <vector>
#include "GL\glut.h"

#define MINPOS -1
#define MAXPOS 1

using namespace std;

enum AirType
{
	UNUSED = 0,
	Smoke
};

class Air
{
public:
	Particle		m_Particle; //�ϳ��� ����
	AirType		m_Type;
public:
	double			m_Viscosity; //����
	double			m_Pressure; //�з�
	double			m_Density; //�е�
public:
	Air(void) {};
	Air(AirType);
	~Air(void) {};
public:
	//�����������
	void Init(void); //�ʱⰪ ����
public:
	//���
	void SetPressure(double newP);
};

class NS2D
{
public:
	vector<Air>		m_Air;
	AirType			m_CurrentType;
	int					m_GridX;
	int					m_GridY;
public:
	NS2D();
	NS2D(int,int);
	~NS2D();
public:
	void Init(void);
	void Update(void);
public:
	void Draw(void);
};

#endif
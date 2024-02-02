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
	Particle		m_Particle; //하나의 입자
	AirType		m_Type;
public:
	double			m_Viscosity; //점성
	double			m_Pressure; //압력
	double			m_Density; //밀도
public:
	Air(void) {};
	Air(AirType);
	~Air(void) {};
public:
	//사용하지않음
	void Init(void); //초기값 설정
public:
	//사용
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
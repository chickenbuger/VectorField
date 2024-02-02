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
	vector<Particle*>	m_Particle; //입자들
public:
	double						m_Viscosity; //점성
	double						m_Pressure; //압력
	double						m_Density; //밀도
public:
	Smoke(void) {};
	~Smoke(void) {};
public:
	void Start(void); //연기 생성 시작
	void Update(void);
	void Init(void); //초기값 설정
public:
	void SetPressure(double newP);
public:
	void Draw(void);
};

#endif
#ifndef __WATERBALL2_H__
#define __WATERBALL2_H__

#pragma once

#include "Vec3.h"
#include "Particle.h"
#include <vector>
#include "GL\glut.h"

using namespace std;

class WaterBall2
{
public:
	Particle					Center;
	vector<Particle*>	DrawParticle;
	double						R;
	double						InvDensity;
	double						HowMuchIn;
public:
	WaterBall2();
	~WaterBall2();
public:
	void Init();
	void Update(double dt);
	void ApplyExternalForce(double dt);
	void CheckInWater();
	void DrawCircle();
	void Bending();
public:
	void Draw();
};

#endif
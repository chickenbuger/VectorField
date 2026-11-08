#ifndef __WATERBALL_H__
#define __WATERBALL_H__

#pragma once

#include "Vec3.h"
#include "Particle.h"
#include <vector>
#include "GL\glut.h"

using namespace std;

enum class Mode
{
	Water = 0,
	Floor
};

class WaterBall
{
public:
	Particle					Center;
	vector<Particle*>	DrawParticle;
	double						R;
	double						InvDensity;
	double						HowMuchIn;
	Mode						CurrentMode;
public:
	WaterBall();
	~WaterBall();
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
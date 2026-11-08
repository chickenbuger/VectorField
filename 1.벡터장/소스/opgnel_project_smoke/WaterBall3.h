#pragma once

#include "Vec3.h"
#include "Particle.h"
#include <vector>
#include "VectorField2D.h"
#include "GL/glut.h"

using namespace std;

class WaterBall3
{
public:
	Particle	Center;
	Vec3<double> newPos;
	double		springR;
	double		y0;
	double		diff;	
public:
	VectorField2D* vF;
public:
	WaterBall3();
	~WaterBall3();
public:
	void Init();
	void Update(double dt);
public:
	void Draw();
};
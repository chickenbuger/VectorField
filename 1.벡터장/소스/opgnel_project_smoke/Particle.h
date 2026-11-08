#ifndef __PARTICLE_H__
#define __PARTICLE_H__

#pragma once

#include "Vec3.h"
#include <math.h>

//#define VF1
//#define VF2
//#define VF3
//#define VF4
#define VF5

class Particle
{
public:
	double				m_InverseMass; // mass, F=ma -> a = F/m
	double				m_Damping; // -velocity)
	Vec3<double>	m_Position; // pos 
	Vec3<double>	m_Velocity; // vel
	Vec3<double>	m_Acceleration; // acc
	Vec3<double>	m_Force; // force, a = f/m
	Vec3<double>	m_VField;
public:
	Particle();
	~Particle();
public:
	inline void			SetMass(const double mass) { m_InverseMass = 1.0 / mass; }
	inline void			SetInverseMass(const double inverseMass) { m_InverseMass = inverseMass; }
	inline void			SetDamping(const double damping) { m_Damping = damping; }
	inline void			SetPos(const Vec3<double> pos) { m_Position = pos; }
	inline void			SetPos(const double x, const double y, const double z) { m_Position.Set(x, y, z); }
	inline void			SetVel(const Vec3<double> vel) { m_Velocity = vel; }
	inline void			SetVel(const double x, const double y, const double z) { m_Velocity.Set(x, y, z); }
	inline void			SetAcc(const Vec3<double> acc) { m_Acceleration = acc; }
	inline void			SetAcc(const double x, const double y, const double z) { m_Acceleration.Set(x, y, z); }
	inline void			SetVField(const Vec3<double> vf) { m_VField = vf; }
	inline void			SetVField(const double x, const double y, const double z) { m_VField.Set(x, y, z); }
public:
	inline double				GetMass(void) { return 1.0 / m_InverseMass; }
	inline double				GetInverseMass(void) { return m_InverseMass; }
	inline double				GetDamping(void) { return m_Damping; }
	inline Vec3<double>	GetPos(void) { return m_Position; }
	inline Vec3<double>	GetVel(void) { return m_Velocity; }
	inline Vec3<double>	GetAcc(void) { return m_Acceleration; }
	inline Vec3<double>	GetVField(void) { return m_VField; }
public:
	void	ClearForce(void) { m_Force.Clear(); }
	void	AddForce(const Vec3<double> force) { m_Force += force; }
	bool	HasFiniteMass(void) { return m_InverseMass >= 0.0; }
	void	Integrate(double dt);
	void	ComputeVector2D();
	void	ComputeVector3D();
};

#endif
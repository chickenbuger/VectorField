#include "Particle.h"

Particle::Particle()
{
}

Particle::~Particle()
{
}

void Particle::Integrate(double dt) //ХыЧе
{
	if (m_InverseMass <= 0.0f)
	{
		return;
	}

	// accumulating acc : a = f/m (Newton's 2)
	Vec3<double> acc = m_Acceleration;
	acc += m_Force * m_InverseMass; // a=F/m

	m_Velocity += acc * dt; // intgrate by time
	m_Velocity *= pow(m_Damping, dt); // resistance force
	m_Position += m_Velocity * dt; // integrate by time

	ClearForce();
}

void Particle::ComputeVector2D()
{
	double PosX = m_Position.x();
	double PosY = m_Position.y();

#ifdef VF1 
	//F(x,y) = <-y,x>
	double Fx = -PosY;
	double Fy = PosX;
#endif

#ifdef VF2 
	//F(x,y) = <xx,yy>
	double Fx = powf(PosX, 2);
	double Fy = powf(PosY, 2);
#endif

#ifdef VF3 
	//F(x,y) = <xxx,yyy>
	double Fx = powf(PosX, 3);
	double Fy = powf(PosY, 3);
#endif

#ifdef VF4
	double Fx = PosY;
	double Fy = PosX;
#endif

#ifdef VF5
	double Fx = PosY;
	double Fy = sin(Fx);
#endif

	SetVField(Fx, Fy, 0);
}

void Particle::ComputeVector3D()
{
	double PosX = m_Position.x();
	double PosY = m_Position.y();
	double PosZ = m_Position.z();

#ifdef VF1 
	double Fx = PosY;
	double Fy = PosZ;
	double Fz = PosX;
#endif

#ifdef VF2
	double Fx = powf(PosX, 2);
	double Fy = powf(PosY, 2);;
	double Fz = powf(PosZ, 2);;
#endif

#ifdef VF3 
	double Fx = PosY;
	double Fy = -2;
	double Fz = PosX;
#endif

#ifdef VF4 
	double Fx = PosY / PosZ;
	double Fy = -PosX / PosZ;
	double Fz = PosZ / 4;
#endif

#ifdef VF5 
	double Fx = PosX / sqrt(powf(PosX, 2) + powf(PosY, 2) + powf(PosZ, 2));
	double Fy = PosY / sqrt(powf(PosX, 2) + powf(PosY, 2) + powf(PosZ, 2));
	double Fz = PosZ / sqrt(powf(PosX, 2) + powf(PosY, 2) + powf(PosZ, 2));
#endif

	SetVField(Fx, Fy, Fz);
}

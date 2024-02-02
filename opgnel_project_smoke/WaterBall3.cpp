#include "WaterBall3.h"

WaterBall3::WaterBall3()
{
	Center.SetPos(0.0, 1.0, 0.0);
	Center.SetInverseMass(1.0f);
	Center.SetDamping(0.99);

	springR = 0.3f;
	y0 = Center.GetPos().y();
	diff = 0.0;
	Init();
}

WaterBall3::~WaterBall3()
{

}

void WaterBall3::Init()
{
	vF = new VectorField2D(10, 10);
	vF->m_diff = 2.0;
	vF->Init(15,15);
}

void WaterBall3::Update(double dt)
{
	//gravity
	Vec3<double> gravity(0.0, -9.8, 0.0);
	double damping = 0.99;

	//Apply gravity
	Center.SetVel(Center.GetVel() + gravity * dt * Center.GetInverseMass());
	Center.SetVel(Center.GetVel() * damping);
	newPos = Center.GetPos() + Center.GetVel() * dt;

	//Bending
	double currentDiff = (Vec3<double>(0, 0, 0) - Center.GetPos()).Length() - diff;
	Vec3<double> dp = Vec3<double>(0, 0, 0) - Center.GetPos();
	dp.Normalize();
	dp *= Center.GetInverseMass() / (Center.GetInverseMass() + 1/10000000) * currentDiff;
	newPos += dp;

	//Integrate
	Center.SetVel((newPos-Center.GetPos())/dt);
	Center.SetPos(newPos);

	//Test
	//printf("Center Pos %.2f %.2f %.2f\n", Center.GetPos().x(), Center.GetPos().y(), Center.GetPos().z());
	//printf("Center ACC %.2f %.2f %.2f\n", Center.GetAcc().x(), Center.GetAcc().y(), Center.GetAcc().z());
	//printf("dp : %.2f %.2f %.2f", dp.GetX(), dp.GetY(), dp.GetZ());

	/*
	Center.SetAcc(Center.GetAcc() + gravity); //중력가속도가 쌓임

	//Integrate
	Center.Integrate(dt);
	*/
}

void WaterBall3::Draw()
{
	//printf("WB3 Draw Call \n");
	glColor3f(1, 0, 0);
	glPointSize(5.0f);
	glBegin(GL_POINTS);
	for (double i = 0; i < 360; i++)
	{
		glVertex3f(cos(i) * springR + Center.GetPos().x(), sin(i) * springR + Center.GetPos().y(), 0);
	}
	//glVertex3f(Center.GetPos().x(), Center.GetPos().y(), Center.GetPos().z());
	glEnd();

	//Direction Point
	glColor3f(0, 1, 0);
	glPointSize(3.0f);
	glBegin(GL_POINTS);
	for (auto& p : vF->m_Particles)
	{
		double df = p->GetPos().y() - 0.0f;
		df /= 10.f;
		glVertex3f(p->GetPos().x(), p->GetPos().y() - df, p->GetPos().z());
	}
	glEnd();
	glPointSize(1.0f);

	//Vector Field
	glColor3f(1, 0, 0);
	glBegin(GL_LINES);
	for (auto& p : vF->m_Particles)
	{
		glVertex3f(p->GetPos().x(), p->GetPos().y(), p->GetPos().z());
		double df = p->GetPos().y() - 0.0f;
		df /= 10.f;
		glVertex3f(p->GetPos().x(), p->GetPos().y() - df, p->GetPos().z());
	}
	glEnd();
}

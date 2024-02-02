#include "VectorField3D.h"

VectorField3D::VectorField3D()
{
}

VectorField3D::VectorField3D(int x, int y, int z)
{
	printf("VF3D 생성자 x,y 호출 x : %d y : %d z : %d\n", x, y, z);
	double total_size = x * y * z;

	if (x <= 1 || y <= 1 || z <= 1)
	{
		return;
	}

	for (int i = 0; i < total_size; i++)
	{
		Particle* p = new Particle;
		double dx = 1.6 / (x - 1);
		double dy = 1.6 / (y - 1);
		double dz = 1.6 / (z - 1);
		//printf("dx %f dy %f\n", dx, dy);

		double _x = (i % (x * y)) % x;
		double _y = (i % (x * y)) / x;
		double _z = i/(x*y);

		p->SetMass(1);
		p->SetDamping(1);
		p->SetPos(dx * _x - 0.8, dy * _y - 0.8, dz * _z - 0.8);
		//printf("%d %.2f %.2f %.2f\n", i, p->GetPos().x(), p->GetPos().y(), p->GetPos().z());
		p->SetVel(0, 0, 0);

		m_Particles.push_back(p);
	}
}

VectorField3D::~VectorField3D()
{
}

void VectorField3D::Init(int x, int y, int z)
{
	printf("VF3D Init 함수 호출\n");
	double total_size = x * y * z;

	if (x <= 1 || y <= 1 || z <= 1)
	{
		return;
	}

	for (int i = 0; i < total_size; i++)
	{
		Particle* p = new Particle;
		double dx = 2.0 / (x - 1);
		double dy = 2.0 / (y - 1);
		double dz = 2.0 / (z - 1);
		//printf("dx %f dy %f\n", dx, dy);

		double _x = (i % (x * y)) % x;
		double _y = (i % (x * y)) / x;
		double _z = i / (x * y);

		p->SetMass(1);
		p->SetDamping(1);
		p->SetPos(dx * _x - 1, dy * _y - 1, dz * _z - 1);
		p->SetVel(0, 0, 0);

		m_Particles.push_back(p);
	}
}

void VectorField3D::Update()
{
	printf("Call Update Func");
}

void VectorField3D::SettingForce()
{
	printf("Call SettingForce Func\n");
	//F(x,y) = <-y,x>
	for (auto* p : m_Particles)
	{
		p->ComputeVector3D();
		/*printf("Cur Pos X : %.2f Y : %.2f Z : %.2f ||  VPos X : %.2f Y : %.2f Z : %.2f \n",
			p->GetPos().x(), p->GetPos().y(), p->GetPos().z(), 
			p->GetVField().x(), p->GetVField().y(), p->GetVField().z());*/
	}
}

void VectorField3D::Draw()
{
	printf("VF3D Draw 함수 호출 \n");

	SettingForce();

	int i = 0;
	for (auto& p : m_Particles)
	{
		auto pos = p->GetPos();
		auto Vpos = p->GetVField();
		float Color = Vpos.Length();
		Color = 1 / Color;
		//printf("Color : %.2f\n", Color);
		/*
		Vec3<double> Dir = Vpos - pos;
		Dir.Normalize();
		Dir /= 15.0f;
		*/
		//printf("%d : 반복문 pos : %f, %f \n", i++, pos.x(), pos.y());
		glColor3f(Color - 0.5, 1 - Color, 0.5);
		glBegin(GL_LINES);
			glVertex3f(pos.x(), pos.y(), pos.z());
			glVertex3f(pos.x() + Vpos.x() , pos.y() + Vpos.y() , pos.z() + Vpos.z() );
			//glVertex3f(pos.x() + Vpos.x() / 3.0f, pos.y() + Vpos.y() / 3.0f, pos.y() + Vpos.z() / 3.0f);
		glEnd();
		glPointSize(2.0f);
		glBegin(GL_POINTS);
			//glVertex3f(pos.x(), pos.y(), pos.z());
			glColor3f(1, 0, 0);
			glVertex3f(pos.x() + Vpos.x(), pos.y() + Vpos.y(), pos.z() + Vpos.z());
			//glVertex3f(pos.x() + Vpos.x() / 3.0f, pos.y() + Vpos.y() / 3.0f, pos.z() + Vpos.z() / 3.0f);
		glEnd();
	}
}
#include "VectorField2D.h"

VectorField2D::VectorField2D()
{
	m_diff = 1.0f;
}

VectorField2D::VectorField2D(int x, int y)
{
	printf("생성자 x,y 호출 x : %d y : %d\n", x, y);
	double total_size = x * y;
	m_diff = 1.0f;
	if (x <= 1 || y <= 1)
	{
		return;
	}

	for (int i = 0; i < total_size; i++)
	{
		Particle* p = new Particle;
		double dx = (2.0 * m_diff) / (x - 1);
		double dy = (2.0 * m_diff) / (y - 1);
		//printf("dx %f dy %f\n", dx, dy);
		p->SetMass(1);
		p->SetDamping(1);
		p->SetPos(dx * (i % x) - m_diff, dy * (i / x) - m_diff, 0); //-1~1사이 점을 분포시킴
		//printf("%d %.2f %.2f\n", i,p->GetPos().x(),p->GetPos().y());
		p->SetVel(0, 0, 0);
		
		m_Particles.push_back(p);
	}
}

VectorField2D::~VectorField2D()
{
}

void VectorField2D::Init(int x, int y)
{
	printf("VF Init 함수 호출\n");
	m_Particles.clear();
	double total_size = x * y;
	if (x <= 1 || y <= 1)
	{
		return;
	}

	for (int i = 0; i < total_size; i++)
	{
		Particle* p = new Particle;
		double dx = (2.0 * m_diff) / (x - 1);
		double dy = (2.0 * m_diff) / (y - 1);
		//printf("dx %f dy %f\n", dx, dy);
		p->SetMass(1);
		p->SetDamping(1);
		p->SetPos(dx * (i % x) - m_diff, dy * (i / x) - m_diff, 0); //-1~1사이 점을 분포시킴
		//printf("%d %.2f %.2f\n", i,p->GetPos().x(),p->GetPos().y());
		p->SetVel(0, 0, 0);

		m_Particles.push_back(p);
	}
}

void VectorField2D::Update()
{
	printf("Call Update Func");
}

void VectorField2D::SettingForce()
{
	printf("Call SettingForce Func\n");
	//F(x,y) = <-y,x>
	for (auto* p : m_Particles)
	{
		p->ComputeVector2D();
	}
}

void VectorField2D::Draw()
{
	//printf("VF2D Draw 함수 호출 \n");

	SettingForce();

	int i = 0;
	for (auto& p : m_Particles)
	{
		auto pos = p->GetPos();
		auto Vpos = p->GetVField();
		float Color = Vpos.Length();
		Color = 1 / Color;
		/*
		Vec3<double> Dir = Vpos - pos;
		Dir.Normalize();
		Dir /= 15.0f;
		*/
		//printf("%d : 반복문 pos : %f, %f \n", i++, pos.x(), pos.y());
		glColor3f(Color - 0.5, 1 - Color, 0.5);
		glBegin(GL_LINES);
			glVertex2f(pos.x(), pos.y());	
			glVertex2f(pos.x() + Vpos.x(), pos.y() + Vpos.y());
			//glVertex2f(pos.x() + Vpos.x() / 3.0f, pos.y() + Vpos.y() / 3.0f);
		glEnd();
		glPointSize(4.0f);
		glBegin(GL_POINTS);
			glColor3f(1.0, 1.0, 0.0);
			glVertex2f(pos.x() + Vpos.x(), pos.y() + Vpos.y());
			//glVertex2f(pos.x() + Vpos.x() / 3.0f, pos.y() + Vpos.y() / 3.0f);
		glEnd();
	}
}

void VectorField2D::PlanetForceDraw()
{
	int i = 0;
	for (auto& p : m_Particles)
	{
		auto pos = p->GetPos();
		auto Vpos = p->GetVField();
		float Color = Vpos.Length();
		Color = 1 / Color;

		glColor3f(Color - 0.5, 1 - Color, 0.5);
		glBegin(GL_LINES);
		glVertex2f(pos.x(), pos.y());
		glVertex2f(pos.x() + Vpos.x(), pos.y() + Vpos.y());

		glEnd();
		glPointSize(4.0f);
		glBegin(GL_POINTS);
		glColor3f(1.0, 1.0, 0.0);
		glVertex2f(pos.x() + Vpos.x(), pos.y() + Vpos.y());

		glEnd();
	}
}

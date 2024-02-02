#include "Planet.h"

Planet::Planet()
{
	_vF = new VectorField2D(15, 15);
	_vF->m_diff = 2.0f;
	_vF->Init(15, 15);
}

Planet::~Planet()
{
}

void Planet::Update()
{
	//각 입자에 대하여 외력 계산  => VField에 저장
	for (auto& _p : _vF->m_Particles)
	{
		for (auto& _g : _gP)
		{
			Vec3<double> pos1 =	_g->_pos;
			Vec3<double> pos2 = _p->GetPos();
			Vec3<double> dir = pos1 - pos2;
			double distance = dir.Length();
			dir.Normalize();

			double gF = G * (_p->GetMass() * _g->_mass) / distance;

			_p->SetVField((_p->GetVField() + (dir * gF)) / (_gP.size()/2.0f));
			/*
			printf("pos1 %.2f %.2f %.2f pos2 %.2f %.2f %.2f dir %.2f %.2f %.2f dis %.2f gF %.2f",
				pos1.GetX(), pos1.GetY(), pos1.GetZ(), pos2.GetX(), pos2.GetY(), pos2.GetZ(),
				dir.GetX(), dir.GetY(), dir.GetZ(), distance, gF);
				*/
		}
	}
}

void Planet::CreatePlanet()
{
	GPlanet* g = new GPlanet();

	_gP.push_back(g);
}

void Planet::CreatePlanet(double x, double y)
{
	GPlanet* g = new GPlanet();
	
	g->_pos = Vec3<double>(x, y, 0);

	_gP.push_back(g);
}

void Planet::CreatePlanet(Vec3<double> pos)
{
	GPlanet* g = new GPlanet();

	g->_pos = pos;

	_gP.push_back(g);
}

void Planet::Draw()
{
	//Draw Planet
	for (auto& gP : _gP)
	{
		gP->Draw();
	}

	//VF Point
	_vF->PlanetForceDraw();
}

GPlanet::GPlanet()
{
	_r = 0.2f;
	_mass = 1.0f;
	_gravity = 9.8f;
	_pos.Set(0, 0, 0);
}

GPlanet::GPlanet(double x, double y)
{
	_r = 0.2f;
	_mass = 1.0f;
	_gravity = 9.8f;
	_pos.Set(x, y, 0);
}

GPlanet::GPlanet(double x, double y, double z)
{
	_r = 0.2f;
	_mass = 1.0f;
	_gravity = 9.8f;
	_pos.Set(x, y, z);
}

GPlanet::GPlanet(Vec3<double> pos)
{
	_r = 0.2f;
	_mass = 1.0f;
	_gravity = 9.8f;
	_pos.Set(pos);
}

GPlanet::~GPlanet()
{
}

void GPlanet::Draw()
{
	//printf("GP draw함수 호출\n");
	glColor3f(0, 1, 0);
	glBegin(GL_POLYGON);
	for (int i = 0; i < 360; i++)
	{
		double angle = i * 3.141592 / 180;
		double x = _r * cos(angle);
		double y = _r * sin(angle);
		glVertex3f(x + _pos.x(), y + _pos.y(), 0);
	}
	glEnd();
}

#include "NS2D.h"

void Air::SetPressure(double newP)
{
	m_Density = m_Pressure = newP;
}

Air::Air(AirType Air)
{
	m_Type = Air;
}

void Air::Init()
{
	switch (m_Type)
	{
	case Smoke:
		//연기의 점성 (섭씨 20도 : 1.86x10^-6) -> 0에 수렴 또는 0.00000186
		m_Viscosity = 0;
		//압력은 default를 1로 설정 (대기압)
		m_Pressure = 1;
		//PV=nRT 기체의 상태방정식
		//T=20, P=1, V=`1, R=1, n=1/20
		//P=n -> V=1/n, D=n
		//V*D=1
		m_Density = m_Pressure;
		break;
	}
	m_Particle.SetPos(0, 0, 0);
	m_Particle.ClearForce();
}

NS2D::NS2D()
{
	m_CurrentType = Smoke;
	m_Air.clear();

	for (auto& p : m_Air)
	{
		p.m_Type = m_CurrentType;
	}
}

NS2D::NS2D(int x, int y)
{
	m_GridX = x;
	m_GridY = y;
	Init();
}

void NS2D::Init(void)
{
	switch (m_CurrentType)
	{
	case Smoke:
		for (auto& p : m_Air)
		{
			//연기의 점성 (섭씨 20도 : 1.86x10^-6) -> 0에 수렴 또는 0.00000186
			p.m_Viscosity = 0;
			//밀도 = 질량 / 부피 (부피는 1로 가정)
			p.m_Density= p.m_Particle.GetMass();
			//밀도= 질량*압력 / 기체 상수*온도 => 20R(기체상수) = 압력 => 기체상수는 1/20으로 가정
			//공기의 특별 기체상수 = 287.058J / kg*K
			//대기의 밀도 : T= 20인경우 1.204
		}
		break;
	}
}

void NS2D::Update(void)
{
	double dt = 0.001;

	for (auto& p : m_Air)
	{
		p.m_Particle.Integrate(dt); //업데이트
	}
}

void NS2D::Draw(void)
{
}

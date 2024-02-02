#include "Smoke.h"

void Smoke::Update(void)
{
	double dt = 0.001;

	for (auto& p : m_Particle)
	{
		p->Integrate(dt);
	}
}

void Smoke::Init(void)
{
	//연기의 점성 (섭씨 20도 : 1.86x10^-6) -> 0에 수렴 또는 0.00000186
	m_Viscosity = 0;
	//압력은 default를 1로 설정 (대기압)
	m_Pressure = 1;
	//PV=nRT 기체의 상태방정식
	//T=20, P=1, V=`1, R=1, n=1/20
	//P=n -> V=1/n, D=n
	//V*D=1
	m_Density = m_Pressure;
}

void Smoke::SetPressure(double newP)
{
	m_Density = m_Pressure = newP;
}

void Smoke::Draw(void)
{
	Vec3<double> pos;
}

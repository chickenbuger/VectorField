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
		//������ ���� (���� 20�� : 1.86x10^-6) -> 0�� ���� �Ǵ� 0.00000186
		m_Viscosity = 0;
		//�з��� default�� 1�� ���� (����)
		m_Pressure = 1;
		//PV=nRT ��ü�� ���¹�����
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
			//������ ���� (���� 20�� : 1.86x10^-6) -> 0�� ���� �Ǵ� 0.00000186
			p.m_Viscosity = 0;
			//�е� = ���� / ���� (���Ǵ� 1�� ����)
			p.m_Density= p.m_Particle.GetMass();
			//�е�= ����*�з� / ��ü ���*�µ� => 20R(��ü���) = �з� => ��ü����� 1/20���� ����
			//������ Ư�� ��ü��� = 287.058J / kg*K
			//����� �е� : T= 20�ΰ�� 1.204
		}
		break;
	}
}

void NS2D::Update(void)
{
	double dt = 0.001;

	for (auto& p : m_Air)
	{
		p.m_Particle.Integrate(dt); //������Ʈ
	}
}

void NS2D::Draw(void)
{
}

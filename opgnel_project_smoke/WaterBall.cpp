#include "WaterBall.h"

WaterBall::WaterBall()
{
	//printf("������ WaterBall\n");
	Center.SetPos(0, 1.0, 0);
	R = 0.3;
	InvDensity = 1;
	HowMuchIn = 0;
	CurrentMode = Mode::Floor;
	//printf("ȣ��\n");
	Init();
}

WaterBall::~WaterBall()
{
}

void WaterBall::Init()
{
	//printf("Init�Լ� ȣ�� Water Ball\n");
	for (double i = 0; i < 360 / InvDensity; i++)
	{
		Particle* p = new Particle();
		
		//�� �׸���(������ �������� ������ R)
		//x=cos(x)*R, y=sin(x)*R
		p->SetPos(cos(i) * R + Center.GetPos().x(), sin(i) * R + Center.GetPos().y(), 0);
		p->SetInverseMass(1.0f);
		p->SetDamping(0.1);

		DrawParticle.push_back(p);
	}
	//printf("Init �Լ�\n DrawParticle size : %d\n", DrawParticle.size());
}

void WaterBall::Update(double dt)
{
	Bending();
	ApplyExternalForce(dt);
	
	for (auto& p : DrawParticle)
	{
		p->Integrate(dt);
	}
}

void WaterBall::ApplyExternalForce(double dt)
{
	
	//Gravity
	Vec3<double> gravity(0.0, -9.8, 0.0);

	for (auto& p : DrawParticle)
	{
		p->SetAcc(p->GetAcc() + gravity);
	}
	Center.SetAcc(Center.GetAcc() + gravity);
	

	//�þ ����
	double StretchY = Center.GetPos().y() - 0.0;
	Vec3<double> dp = Vec3<double>(0, 0, 0) - Center.GetPos();
	dp.Normalize();
	dp *= 1 / 10000 * StretchY;
	Center.SetPos(Center.GetPos() + dp);

	/*
	//printf("���� �߷����� ACC : %f \n", DrawParticle[0]->GetAcc().y());
	//Buoyancy
	// �е� = ����/����  ����=>1 �е�=>1 ���� =>1
	// ���ǰ� �����ϸ� �е��� �����Ѵ�.
	//���Ǵ� R
	//�η� = ���� �е� * ��� ���� * �߷°��ӵ�
	if (CurrentMode == Mode::Water)
	{
		if (HowMuchIn != 0)
		{
			printf("�η½��� ���� H : %.2f  Center�� Pos x : %.2f, y : %.2f, z : %.2f\n", HowMuchIn, Center.GetPos().x(), Center.GetPos().y(), Center.GetPos().z());
			double den = DrawParticle[0]->GetMass() / R;
			Vec3<double> buoyancy(-gravity * den * HowMuchIn / DrawParticle.size());
			printf("���� �η����� ACC : %f\n", DrawParticle[0]->GetAcc().y());

			for (auto& p : DrawParticle)
			{
				p->SetAcc(p->GetAcc() + buoyancy);
			}
			Center.SetAcc(Center.GetAcc() + buoyancy);
		}
	}
	*/
}

void WaterBall::CheckInWater()
{
	int CountValue = 0;
	for (auto p : DrawParticle)
	{
		if (p->GetPos().y() <= 0)
		{
			CountValue++;
		}
	}
	HowMuchIn = CountValue;
}


void WaterBall::DrawCircle()
{
	printf("Draw Circle ȣ��\n");
	glColor3f(1, 0, 0);
	glPointSize(3.0f);
	glBegin(GL_POINTS);
	printf("%.2f %.2f %.2f\n", Center.GetPos().x(), Center.GetPos().y(), Center.GetPos().z());
	glVertex3f(Center.GetPos().x(), Center.GetPos().y(), 0);
	glVertex3f(0, 0, 0);
	/*for (double i = 0; i < 360 / InvDensity; i++)
	{
		glVertex3f(cos(i) * R + Center.GetPos().x(), sin(i) * R + Center.GetPos().y(), 0);
	}*/
	glEnd();
	glPointSize(1.0f);
}

void WaterBall::Bending()
{
	//printf("Bending ȣ��, ���� center�� %.2f %.2f %.2f\n", Center.GetPos().x(), Center.GetPos().y(), Center.GetPos().z());
	double StretchY = Center.GetPos().y() - 0.0;
	//printf("Stretch �� : %.2f\n", StretchY);
	Vec3<double> dp = Vec3<double>(0, 0, 0) - Center.GetPos();
	//printf("dp �� %.2f \n", dp.GetY());
	dp.Normalize();
	//printf("dp �� %.2f \n", dp.GetY());
	dp *= 1.0 / 100 * StretchY;
	//printf("dp �� %.4f \n", dp.GetY());
	Center.SetPos(Center.GetPos() + dp);
	//DrawCircle();
	Draw();
}

void WaterBall::Draw()
{
	//printf("WaterBall Draw �Լ� ȣ��");
	glColor3f(1, 0, 0);
	glPointSize(3.0f);
	glBegin(GL_POINTS);
	for (double i = 0; i < 360 / InvDensity; i++)
	{
		//glVertex3f(Center.GetPos().x(), Center.GetPos().y(), 0);
		glVertex3f(cos(i) * R + Center.GetPos().x(), sin(i) * R + Center.GetPos().y(), 0);
	}
	/*for (auto& p : DrawParticle)
	{
		glVertex3f(p->GetPos().x(), p->GetPos().y(), 0);
	}*/
	glEnd();
	glPointSize(1.0f);

	/*
	glPointSize(3.0f);
	glBegin(GL_POINTS);
	glVertex3f(0, 1.5, 0);
	glEnd();
	glPointSize(1.0f);
	*/
}

#include <iostream>
#include "VectorField2D.h"
#include "VectorField3D.h"
#include "WaterBall.h"
#include "WaterBall2.h"
#include "WaterBall3.h"
#include "Planet.h"
#include "GL/glut.h"

//#define VectorField_2D
//#define VectorField_3D
//#define Buoyancy
//#define Buoyancy2
//#define Buoyancy3
#define PlanetGravity

VectorField2D*	_VF2D; 
VectorField3D*	_VF3D;
WaterBall*			_Water;
WaterBall2*		_Water2;
WaterBall3*		_Water3;
Planet*				_Planet;

int width = 800;
int height = 800;

float _zoom = 15.0f; // 화면 확대,축소
float _rot_x = 0.0f; // x축 회전
float _rot_y = 0.001f; // y축 회전
float _trans_x = 0.0f; // x축 이동
float _trans_y = 0.0f; // y축 이동
int _last_x = 0; // 이전 마우스 클릭 x위치
int _last_y = 0; // 이전 마우스 클릭 y위치
unsigned char _buttons[3] = { 0 }; // 마우스 상태(왼쪽,오른쪽,휠 버튼)

void Init(void)
{
	glEnable(GL_DEPTH_TEST);

	//printf("Init 함수 호출\n");

#ifdef VectorField_2D
	_VF2D = new VectorField2D(10, 10);
#endif

#ifdef VectorField_3D
	//_VF3D = new VectorField3D(10, 10, 10);
	_VF3D = new VectorField3D(15, 15, 15);
#endif

#ifdef Buoyancy
	_Water = new WaterBall();
#endif

#ifdef Buoyancy2
	_Water2 = new WaterBall2();
#endif

#ifdef Buoyancy3
	_Water3 = new WaterBall3();
#endif

#ifdef PlanetGravity
	_Planet = new Planet();
#endif
}

void DrawXY()
{
	//x축
	glColor3f(0.0, 0.0, 0.0);
	glLineWidth(3.0f);
	glBegin(GL_LINES);	
		glVertex2f(-2.0f, 0.0f);
		glVertex2f(2.0f, 0.0f);
	glEnd();

	//y축
	glColor3f(0.0, 0.0, 0.0);
	glBegin(GL_LINES);
		glVertex2f(0.0f, -2.0f);
		glVertex2f(0.0f, 2.0f);
	glEnd();
	glLineWidth(1.0f);
}

void DrawXYZ()
{
	//x축
	glColor3f(1.0, 0.0, 0.0);
	glLineWidth(3.0f);
	glBegin(GL_LINES);
	glVertex3f(-2.0f, 0.0f, 0.0f);
	glVertex3f(2.0f, 0.0f, 0.0f);
	glEnd();

	//y축
	glColor3f(0.0, 1.0, 0.0);
	glBegin(GL_LINES);
	glVertex3f(0.0f, -2.0f, 0.0f);
	glVertex3f(0.0f, 2.0f, 0.0f);
	glEnd();

	//z축
	glColor3f(0.0, 0.0, 1.0);
	glBegin(GL_LINES);
	glVertex3f(0.0f, 0.0f, -2.0f);
	glVertex3f(0.0f, 0.0f, 2.0f);
	glEnd();
	glLineWidth(1.0f);
}

void DrawSurface()
{
	glColor3f(0.0, 0.0, 0.501961);
	glLineWidth(3.0f);
	glBegin(GL_LINES);
	glVertex3f(-2.0f, 0.0f, 0.0f);
	glVertex3f(2.0f, 0.0f, 0.0f);
	glEnd();
}

void Update(void)
{
#ifdef Buoyancy
	_Water->CheckInWater();
	_Water->Update(0.0001);
#endif

#ifdef Buoyancy3
	_Water3->Update(0.1f);
#endif

#ifdef PlanetGravity
	//_Planet->Update();
#endif

	glutPostRedisplay();
}

void Draw(void)
{
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);

	glDisable(GL_LIGHTING);

#ifdef VectorField_2D
	_VF2D->Draw();
	DrawXY();
#endif

#ifdef VectorField_3D
	_VF3D->Draw();
	DrawXYZ();
#endif

#ifdef Buoyancy
	_Water->Draw();
	DrawSurface();
#endif

#ifdef Buoyancy2
	_Water2->Draw();
	DrawSurface();
#endif

#ifdef Buoyancy3
	_Water3->Draw();
	DrawSurface();
#endif

#ifdef PlanetGravity
	_Planet->Draw();
#endif
}

void Display(void)
{
	glClearColor(0.3f, 0.3f, 0.3f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();

#ifdef VectorField_3D
	glTranslatef(0, 0, -_zoom);
	glTranslatef(_trans_x, _trans_y, 0);
	glRotatef(_rot_x, 1, 0, 0);
	glRotatef(_rot_y, 0, 1, 0);

	//camera setting
	//glTranslatef(0.0, 1.0, 0.0);
	glRotatef(-45.0f, 1, 0, 0);
	glRotatef(-45.0f, 0, 1, 0);
	glRotatef(90.0f, 0, 0, 1);
#endif

	Draw();
	glutSwapBuffers(); // GLUT_DOUBLE 
}

void Mouse(int button, int state, int x, int y)
{
#ifdef VectorField_3D
	// 이전 마우스 클릭 위치
	_last_x = x;
	_last_y = y;

	switch (button)
	{
	case GLUT_LEFT_BUTTON:
		_buttons[0] = state == GLUT_DOWN ? 1 : 0;
		break;
	case GLUT_MIDDLE_BUTTON:
		_buttons[1] = state == GLUT_DOWN ? 1 : 0;
		break;
	case GLUT_RIGHT_BUTTON:
		_buttons[2] = state == GLUT_DOWN ? 1 : 0;
		break;
	default:
		break;
	}
	glutPostRedisplay(); // 화면갱신
#endif

#ifdef PlanetGravity
	printf("%f %f\n", x / 200.0f, y / 200.0f);
	y = 800 - y;
	_Planet->CreatePlanet(x / 200.0f - 2, y / 200.0f - 2);

	_Planet->Update();

	glutPostRedisplay();
#endif
}

void Motion(int x, int y)
{
#ifdef VectorField_3D
	int diff_x = x - _last_x;
	int diff_y = y - _last_y;
	_last_x = x;
	_last_y = y;

	if (_buttons[2]) {
		_zoom -= (float)0.02f * diff_x;
	}
	else if (_buttons[1]) {
		_trans_x += (float)0.02f * diff_x;
		_trans_y -= (float)0.02f * diff_y;
	}
	else if (_buttons[0]) {
		_rot_x += (float)0.2f * diff_y;
		_rot_y += (float)0.2f * diff_x;
	}
	glutPostRedisplay(); // 화면갱신
#endif
}


void Reshape(int w, int h)
{
	if (w == 0) 
	{
		h = 1;
	}
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-2.0, 2.0, -2.0, 2.0, -20, 20);
	glMatrixMode(GL_MODELVIEW); 
	glLoadIdentity();
}

void main(int argc, char** argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(width, height);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("Smoke");
	glutDisplayFunc(Display);
	glutReshapeFunc(Reshape);
	glutMouseFunc(Mouse); // 마우스 클릭 이벤트	
	glutMotionFunc(Motion); // 마우스 이동 이벤트
	glutIdleFunc(Update);
	Init();
	glutMainLoop();
}
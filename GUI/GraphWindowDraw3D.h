
/*************************************************************
	GraphWindowDraw3D.h
*************************************************************/

#include "GraphWindow3D.h"
#include "scandata.h"

#ifndef GRAPHWINDOWCF
#define GRAPHWINDOWCF

class GraphWindowDraw3D : public GraphWindow3D
{
protected:
	void draw();
	int handle(int);

	void SetRTS2();
	void draw3DCurveFold();

public:
	int disp_axis, disp_axis2;
	int disp_LIN_SMOOTH;
	int disp_POLY_OFFSET;

	int sts_alt2, sts_sft2, sts_ctrl2;
	int push_button2;
	int push_x2, push_y2;

	scandata *sd;

	CQuaternion quat2;
	float prevquat2[4], currquat2[4], dispquat2[4];
	float prevtrans2[4], currtrans2[4], disptrans2[4];
	float mObject[16];

	GraphWindowDraw3D(int X, int Y, int W, int H, const char *L);
	GraphWindowDraw3D(int X, int Y, int W, int H);
	~GraphWindowDraw3D();

	void init();
	void initTexture();
	void initObject();
	void resetObjTrans();
	void resetObjRot();

	void exportImage(char *filename);
};

#endif	// GRAPHWINDOWCF

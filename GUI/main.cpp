
#if 1
#pragma comment(lib,"opengl32.lib")
//#pragma comment(lib,"glu32.lib")
//#pragma comment(lib,"glut32.lib")
#pragma comment(lib,"freeglut.lib")
#pragma comment(lib,"fltk.lib")
#pragma comment(lib,"fltkforms.lib")
#pragma comment(lib,"fltkgl.lib")
#pragma comment(lib,"fltkimages.lib")
#pragma comment(lib,"fltkjpeg.lib")
#pragma comment(lib,"fltkpng.lib")
#pragma comment(lib,"fltkzlib.lib")
#endif

#ifdef _DEBUG
#pragma comment(lib,"opencv_world420d.lib")
#else
#pragma comment(lib,"opencv_world420.lib")
#endif

#include "GraphWindowDraw3D.h"
#include "GraphWindowDraw2D.h"
#include "ControlPanel.h"
#include "scandata.h"

#define GCF_W 640
#define GCF_H 480
#define GCP_W 320
#define CP_W 280

GraphWindowDraw3D *gwin=NULL;
ControlPanel *cwin=NULL;
GraphWindowDraw2D *gwin_cp=NULL;

int main(int argc, char* argv[])
{
	Fl_Window *window = new Fl_Window( 5, 50, GCF_W + GCP_W + CP_W, GCF_H);
	
	scandata* sdata = new scandata();

	gwin = new GraphWindowDraw3D(0, 0, GCF_W, GCF_H);
	gwin->end();

	gwin_cp = new GraphWindowDraw2D(GCF_W, 0, GCP_W, GCF_H);
	gwin_cp->end();

	cwin = new ControlPanel(GCF_W + GCP_W, 0, CP_W, GCF_H, "aaa", gwin);
	cwin->end();

	gwin_cp->disp_axis = gwin->disp_axis;
	cwin->gwin = gwin;
	cwin->gwin_cp = gwin_cp;
	gwin_cp->cwin = (void*)cwin;
	cwin->sd = gwin->sd = gwin_cp->sd = sdata;
	cwin->vi_axis_y->value(sdata->axis_y);
	cwin->vi_axis_z->value(sdata->axis_z);

	window->add(gwin);
	window->add(gwin_cp);
	window->add(cwin);
	window->end();
	window->show();

	Fl::run();

	return 0;
}



/*************************************************************
	GraphWindowDraw2D.h
*************************************************************/

#include <FL/Fl.h>
#include <FL/Fl_Double_Window.h>
#include <FL/Fl_JPEG_Image.h>
#include <FL/Fl_Box.h>
#include "scandata.h"

#ifndef GRAPHWINDOWCP
#define GRAPHWINDOWCP

class GraphWindowDraw2D : public Fl_Double_Window
{
protected:
	void draw();
	int handle(int);
	void init();

public:
	GraphWindowDraw2D(int X, int Y, int W, int H, const char *L);
	GraphWindowDraw2D(int X, int Y, int W, int H);

	Fl_JPEG_Image *jpg;
	int disp_axis;
	void *cwin;
	scandata* sd;

	double mt[9];
};

#endif // GRAPHWINDOWCP



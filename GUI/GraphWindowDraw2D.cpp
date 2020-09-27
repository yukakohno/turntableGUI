#include <FL/fl_draw.H>
#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include "GraphWindowDraw2D.h"
#include "ControlPanel.h"
#include "util.h"

GraphWindowDraw2D::GraphWindowDraw2D(int X, int Y, int W, int H, const char *L) : Fl_Double_Window(X, Y, W, H, L)
{
	box(FL_FLAT_BOX);
	jpg = NULL;
	init();
	draw();
}

GraphWindowDraw2D::GraphWindowDraw2D(int X, int Y, int W, int H) : Fl_Double_Window(X, Y, W, H)
{
	box(FL_FLAT_BOX);
	jpg = NULL;
	init();
	draw();
}

void GraphWindowDraw2D::init()
{
	if( jpg ){ delete jpg; }
	//jpg = new Fl_JPEG_Image( TEXFNAME );
	sd = NULL;
}

void GraphWindowDraw2D::draw()
{
	/* Draw 2D objects with fltk */

	// background
	fl_color(FL_WHITE);
	fl_rectf(0,0,w(),h());

	int cx = w() / 2, cy = h() / 2;

	if (disp_axis) {
		fl_color(0, 0, 255);
		fl_line_style(FL_SOLID, 3);
		fl_line(cx, cy, cx - 100, cy);
		fl_color(0, 255, 0);
		fl_line_style(FL_SOLID, 3);
		fl_line(cx, cy, cx, cy - 100);
	}

	fl_color(127, 127, 127);
	fl_line_style(FL_SOLID, 3);
	fl_circle(cx, cy, 100);

	if (sd) {
		double sina[MAX_LINE_CNT], cosa[MAX_LINE_CNT], rad[MAX_LINE_CNT];
		for (int i = 0; i < sd->linecnt; i++) {
			sina[i] = sin(sd->lineangle[i]);
			cosa[i] = cos(sd->lineangle[i]);
			rad[i] = 36.0 + sd->lineangv[i] * 10000.0;
		}

		for (int i = 1; i < sd->linecnt; i++) {

			int r, g, b;
			float sec = (float)i / (float)sd->linecnt;

			if (sec < 0.33) {
				r = (int)((0.33 - sec) / 0.33 * 255);
				g = (int)(sec / 0.33 * 255);
				b = 0;
			}
			else if (sec < 0.66) {
				r = 0;
				g = (int)((0.66 - sec) / 0.33 * 255);
				b = (int)((sec - 0.33) / 0.33 * 255);
			}
			else {
				r = (int)((sec - 0.66) / 0.34 * 255);
				g = 0;
				b = (int)((1.0 - sec) / 0.34 * 255);
			}
			fl_color(r, g, b);
			fl_line_style(FL_SOLID, 2);
			fl_line(cx + rad[i - 1] * cosa[i - 1], cy - rad[i - 1] * sina[i - 1],
				cx + rad[i] * cosa[i], cy - rad[i] * sina[i]);
		}
#if 0
		fl_color(255, 0, 0);
		fl_line_style(FL_SOLID, 2);
		fl_line(cx + rad[sd->linecnt - 1] * cosa[sd->linecnt - 1],
			cy - rad[sd->linecnt - 1] * sina[sd->linecnt - 1],
			cx + rad[0] * cosa[0], cy -  rad[0] * sina[0]);
#endif
		fl_color(0, 0, 0);
		fl_line_style(FL_SOLID, 3);
		for (int i = 0; i < sd->cpcnt; i++) {
			int idx = sd->cpidx[i];
			double sina = sin(sd->lineangle[idx]);
			double cosa = cos(sd->lineangle[idx]);
			double rad = 36.0 + sd->cpangv[i] * 10000.0;
			fl_circle(cx + rad * cosa, cy - rad * sina, 5);
		}
	}
}

int GraphWindowDraw2D::handle(int event)
{
	int key, push_x, push_y;

	switch(event){
	case FL_FOCUS:
		return 1; // to detect FL_KEYDOWN & FL_KEYUP
		break;
	case FL_UNFOCUS:
		return 1; // to detect FL_KEYDOWN & FL_KEYUP
		break;
	case FL_KEYDOWN:
		key = Fl::event_key();
		if(key == 'a'){
			return 1;
		}
		break;
	case FL_KEYUP:
		break;
	case FL_PUSH:
		push_x = Fl::event_x();
		push_y = Fl::event_y();
		if( Fl::event_button() == FL_LEFT_MOUSE ){
			std::cout<< "FL_PUSH FL_LEFT_MOUSE x: "<< push_x << " y: " << push_y << std::endl;
		} else if (Fl::event_button() == FL_RIGHT_MOUSE) {
			std::cout << "FL_PUSH FL_RIGHT_MOUSE x: " << push_x << " y: " << push_y << std::endl;
		}
		return 1; // activate FL_DRAG event
		break;
	case FL_DRAG:
		push_x = Fl::event_x();
		push_y = Fl::event_y();
		if(Fl::event_button() == FL_LEFT_MOUSE || Fl::event_button() == FL_RIGHT_MOUSE ){
			std::cout << "FL_DRAG x: " << push_x << " y: " << push_y << std::endl;
			redraw();
		}
		break;
	case FL_RELEASE:
		push_x = Fl::event_x();
		push_y = Fl::event_y();
		if (Fl::event_button() == FL_LEFT_MOUSE) {
			std::cout << "FL_RELEASE FL_LEFT_MOUSE x: " << push_x << " y: " << push_y << std::endl;
		}
		else if (Fl::event_button() == FL_RIGHT_MOUSE) {
			std::cout << "FL_RELEASE FL_RIGHT_MOUSE x: " << push_x << " y: " << push_y << std::endl;
		}
		return 1; // activate FL_DRAG event
		break;
	case FL_MOVE:
		break;
	default:
		break;
	}
	return Fl_Double_Window::handle(event);
}

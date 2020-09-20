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
	jpg = new Fl_JPEG_Image( TEXFNAME );
}

void GraphWindowDraw2D::draw()
{
	/* Draw 2D objects with fltk */

	// background
	fl_color(FL_WHITE);
	fl_rectf(0,0,w(),h());

	if (disp_axis) {
		fl_color(255, 0, 0);
		fl_line_style(FL_SOLID, 3);
		fl_line(0, 0, w(), 0);
		fl_color(0, 255, 0);
		fl_line_style(FL_SOLID, 3);
		fl_line(0, 0, 0, h());
	}

	// sample 2D objects

	if( jpg ){
		jpg->draw( 50, 80, jpg->w(), jpg->h()/2 );
	}

	fl_color(127, 127, 127);
	fl_rectf( 200, 250, 180, 120);

	fl_color(0, 255, 255);
	fl_line_style(FL_SOLID, 2);
	fl_line(50, 80, 200, 250);

	fl_color(255, 0, 255);
	fl_line_style(FL_DOT, 2);
	fl_line(50+jpg->w(), 80 + jpg->h()/2, 200+180, 250+120);

	fl_color(255, 255, 0);
	fl_line_style( FL_SOLID, 3);
	fl_circle( 150, 250, 200 );

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

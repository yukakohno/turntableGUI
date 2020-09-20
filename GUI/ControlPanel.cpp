#include <time.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <Windows.h>
#include "ControlPanel.h"

ControlPanel::ControlPanel(int X, int Y, int W, int H, GraphWindowDraw3D *_gwin )
: Fl_Window(X, Y, W, H)
{
	ControlPanel(X, Y, W, H, "", _gwin );
}

ControlPanel::ControlPanel(int X, int Y, int W, int H, const char *L, GraphWindowDraw3D *_gwin )
: Fl_Window(X, Y, W, H, L)
{
	gwin = _gwin;
	sd = NULL;
	createPanel();
}

ControlPanel::~ControlPanel()
{
	gwin = NULL;
}

void ControlPanel::cb_rb_gwin(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	GraphWindowDraw3D* gwin = This->gwin;
	char wgtlbl[256];	strcpy(wgtlbl, wgt->label());

	if (!strcmp(wgtlbl, "Y-DN") || !strcmp(wgtlbl, "Z-UP") || !strcmp(wgtlbl, "X-UP"))
	{
		// reset rotation: up = 0:y-up, 1:z-up, 2:x-up
		if (!strcmp(wgtlbl, "Y-DN")) {
			wgt->label("Z-UP");
			gwin->resetRotation(1);
		}
		else if (!strcmp(wgtlbl, "Z-UP")) {
			wgt->label("X-UP");
			gwin->resetRotation(2);
		}
		else if (!strcmp(wgtlbl, "X-UP")) {
			//gwin->resetRotation(0);
			wgt->label("Y-DN");
			gwin->resetRotation(-1);
		}
	}
	else if (!strcmp(wgtlbl, "SCL")) {	// reset scale
		gwin->resetScale();
	}
	else if (!strcmp(wgtlbl, "TRNS")) {	// reset translation
		gwin->resetTranslation();
	}
	gwin->redraw();
}

void ControlPanel::createPanel()
{
	int wgt_x=0, wgt_y=0;

	Fl_Tabs* tab = new Fl_Tabs(0, 0, this->w(), this->h());
	{
		Fl_Group* g = new Fl_Group(0, 20, this->w(), this->h()-20, "D0");
		//g->hide();
		{ 
			wgt_x = 10;
			wgt_y = 20;

			// ------------------------- fILE_IO -------------------------------------------
			Fl_Box *bx_FILE_I = new Fl_Box(0, wgt_y, g->w(), 20, "--- LOAD/SAVE ---");

			wgt_y += 20;

			btn_load = new Fl_Button(wgt_x, wgt_y, 80, 20, "load");
			btn_load->callback(cb_btn_load, (void*)this);
			fc = new Fl_File_Chooser(DIR_IN, "(*.{txt})", Fl_File_Chooser::SINGLE, "Fl_File_Chooser");
			g->add(btn_load);

			wgt_x += 90;

			btn_savescreen = new Fl_Button(wgt_x, wgt_y, 80, 20, "save screen");
			btn_savescreen->callback(cb_btn_savescreen, (void*)this);
			g->add(btn_savescreen);

			wgt_x = 10;
			wgt_y += 25;

			// ------------------------- DISPLAY -------------------------------------------
			Fl_Box *bx_DISPLAY= new Fl_Box(0, wgt_y, g->w(), 20, "--- DISPLAY ---");
			wgt_y += 20;

			btn_resetRot = new Fl_Button(wgt_x, wgt_y, 45, 20, "Y-DN");
			btn_resetRot->callback(cb_rb_gwin, (void*)this);
			wgt_x += 45;
			btn_resetScale = new Fl_Button(wgt_x, wgt_y, 45, 20, "SCL");
			btn_resetScale->callback(cb_rb_gwin, (void*)this);
			wgt_x += 45;
			btn_resetTrans = new Fl_Button(wgt_x, wgt_y, 45, 20, "TRNS");
			btn_resetTrans->callback(cb_rb_gwin, (void*)this);
			g->add(btn_resetRot);
			g->add(btn_resetScale);
			g->add(btn_resetTrans);

			wgt_x = 10;
			wgt_y += 25;

			cb_dispaxis = new Fl_Check_Button(wgt_x, wgt_y, 45, 20, "AXIS (OBJECT)");
			cb_dispaxis->value(this->gwin->disp_axis);
			cb_dispaxis->callback(cb_cb_dispaxis, (void*)this);
			wgt_y += 25;
			cb_dispaxis2 = new Fl_Check_Button(wgt_x, wgt_y, 45, 20, "AXIS (LRF)");
			cb_dispaxis2->value(this->gwin->disp_axis2);
			cb_dispaxis2->callback(cb_cb_dispaxis2, (void*)this);
			wgt_y += 25;
			//cb_dispLSMT = new Fl_Check_Button(wgt_x, wgt_y, 45, 20, "Line smooth");
			//cb_dispLSMT->value(this->gwin->disp_LIN_SMOOTH);
			//cb_dispLSMT->callback(cb_cb_dispLSMT, (void*)this);
			//wgt_y += 25;
			//cb_dispOFST = new Fl_Check_Button(wgt_x, wgt_y, 45, 20, "Line offset");
			//cb_dispOFST->value(this->gwin->disp_POLY_OFFSET);
			//cb_dispOFST->callback(cb_cb_dispOFST, (void*)this);

			g->add(cb_dispaxis);
			g->add(cb_dispaxis2);
			//g->add(cb_dispLSMT);
			//g->add(cb_dispOFST);

			wgt_x = 10;
			wgt_y += 25;

			// ------------------------- CONTROL -------------------------------------------
			Fl_Box* bx_CONTROL = new Fl_Box(0, wgt_y, g->w(), 20, "--- CONTROL ---");

			wgt_x = 60;
			wgt_y += 20;

			vi_axis_y = new Fl_Value_Input(wgt_x, wgt_y, 75, 20, "axis_y");
			vi_axis_y->value(0.0);
			vi_axis_y->callback(cb_vi_axis, (void*)this);
			g->add(vi_axis_y);

			wgt_x += 125;

			vi_axis_z = new Fl_Value_Input(wgt_x, wgt_y, 75, 20, "axis_z");
			vi_axis_z->value(0.0);
			vi_axis_z->callback(cb_vi_axis, (void*)this);
			g->add(vi_axis_z);

			wgt_x = 60;
			wgt_y += 25;

			vs_cpidx = new Fl_Value_Slider(wgt_x, wgt_y, 200, 20, "cp #");
			vs_cpidx->bounds(0, 8);	vs_cpidx->step(1);	vs_cpidx->value(0);
			vs_cpidx->align(FL_ALIGN_LEFT);
			vs_cpidx->type(FL_HORIZONTAL);
			vs_cpidx->callback(cb_vs_cpidx, (void*)this);
			g->add(vs_cpidx);

			wgt_y += 25;

			vs_cpangv = new Fl_Value_Slider(wgt_x, wgt_y, 200, 20, "cp angv");
			vs_cpangv->bounds(-1, 1);	vs_cpangv->step(0.0001);	vs_cpangv->value(0);
			vs_cpangv->align(FL_ALIGN_LEFT);
			vs_cpangv->type(FL_HORIZONTAL);
			vs_cpangv->callback(cb_vs_cpangv, (void*)this);
			g->add(vs_cpangv);


		}
		g->end();
	}
	{
		Fl_Group* g = new Fl_Group(0, 20, this->w(), this->h()-20, "D1");
		g->hide();
		{ 
			wgt_x = 10;	wgt_y = 20;

		}
		g->end();
	}
}

int ControlPanel::handle(int event)
{
	int i, key, pos=-1, prm=-1;

	switch(event){
	case FL_FOCUS:
		//return 1; // to detect FL_KEYDOWN & FL_KEYUP
		break;
	case FL_UNFOCUS:
		//return 1; // to detect FL_KEYDOWN & FL_KEYUP
		break;
	case FL_KEYDOWN:
		//printf("FL_KEYDOWN\n");
		key = Fl::event_key();
		if(key=='a'){

		}
		break;
	default:
		//return Fl_Window::handle(event);
		break;
	}
	return Fl_Window::handle(event);
}

#define _USE_MATH_DEFINES
#include <math.h>
#include <sys\stat.h>
#include <direct.h>
#include <time.h>

#include "ControlPanel.h"
#include "util.h"

void ControlPanel::cb_vi_axis(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	This->sd->axis_y = (double)This->vi_axis_y->value();
	This->sd->axis_z = (double)This->vi_axis_z->value();
	This->gwin->redraw();
	This->gwin_cp->redraw();
}

void ControlPanel::cb_vs_cpidx(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	int cpidx = (int)This->vs_cpidx->value();
	This->vs_cpangv->value(This->sd->cpangv[cpidx]);
}
void ControlPanel::cb_vs_cpangv(Fl_Widget* wgt, void* idx)
{
	ControlPanel* This = (ControlPanel*)idx;
	int cpidx = (int)This->vs_cpidx->value();
	double cpval = (double)This->vs_cpangv->value();
	This->sd->cpangv[cpidx] = cpval;
	This->sd->cpangv2lineangle();
	This->gwin->redraw();
	This->gwin_cp->redraw();
}

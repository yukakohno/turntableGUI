#include "ControlPanel.h"

void ControlPanel::cb_cb_dispaxis(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	GraphWindowDraw3D *gwin = This->gwin;
	GraphWindowDraw2D *gwin_cp = This->gwin_cp;
	gwin->disp_axis = gwin_cp->disp_axis = ((Fl_Check_Button*)wgt)->value();
	gwin->redraw();
	gwin_cp->redraw();
}

void ControlPanel::cb_cb_dispaxis2(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	GraphWindowDraw3D *gwin = This->gwin;
	//GraphWindowDraw2D *gwin_cp = This->gwin_cp;
	gwin->disp_axis2 = ((Fl_Check_Button*)wgt)->value();
	gwin->redraw();
	//gwin_cp->redraw();
}

void ControlPanel::cb_cb_dispLSMT(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	GraphWindowDraw3D *gwin = This->gwin;
	gwin->disp_LIN_SMOOTH = ((Fl_Check_Button*)wgt)->value();
	gwin->redraw();
}

void ControlPanel::cb_cb_dispOFST(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	GraphWindowDraw3D *gwin = This->gwin;
	gwin->disp_POLY_OFFSET = ((Fl_Check_Button*)wgt)->value();
	gwin->redraw();
}


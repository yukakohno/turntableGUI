
/*************************************************************
	ControlPanel.h
*************************************************************/

#include <FL/Fl.h>
#include <FL/Fl_Window.h>
#include <FL/Fl_Box.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_File_Chooser.H>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_Group.h>
#include <FL/Fl_Int_Input.H>
#include <FL/Fl_Round_Button.h>
#include <FL/Fl_Value_Slider.h>
#include <FL/Fl_Value_Input.H>
#include <FL/Fl_Tabs.H>
#include "GraphWindowDraw3D.h"
//#include "GraphWindowDraw2D.h"

#define DIR_IN "./input/"
#define DIR_OUT "./output/"
#define TEXFNAME "./texture/grid_bw3.jpg"
#define TEXWIDTH  256	// texture width
#define TEXHEIGHT 256	// texture height

#ifndef CPANEL
#define CPANEL

class ControlPanel : public Fl_Window
{
public:
	GraphWindowDraw3D *gwin;
	//GraphWindowDraw2D *gwin_cp;
	scandata *sd;

	ControlPanel(int X, int Y, int W, int H, GraphWindowDraw3D *_gwin );
	ControlPanel(int X, int Y, int W, int H, const char *L, GraphWindowDraw3D *_gwin );
	~ControlPanel();
	void initPanel();
	void createPanel();
	void refresh(int init);
	void setorg();

private:
	int handle(int event);
	int value_grpfix();
	int value_grpparam();

	// ------------------------- fILE_IO -------------------------------------------
public:
	Fl_Button *btn_load;
	Fl_File_Chooser *fc;
	Fl_Button *btn_loadtex;
	Fl_Button *btn_savescreen;

private:
	static void cb_btn_load( Fl_Widget *wgt, void *idx);
	static void cb_btn_savescreen( Fl_Widget *wgt, void *idx);

	// ------------------------- DISPLAY -------------------------------------------
public:
	Fl_Button *btn_resetRot;
	Fl_Button *btn_resetScale;
	Fl_Button *btn_resetTrans;
	Fl_Check_Button *cb_dispaxis;
	Fl_Check_Button* cb_dispaxis2;
	Fl_Check_Button* cb_dispLSMT;
	Fl_Check_Button* cb_dispOFST;
private:
	static void cb_rb_gwin(Fl_Widget* wgt, void* idx);
	static void cb_cb_dispaxis(Fl_Widget *wgt, void *idx);
	static void cb_cb_dispaxis2(Fl_Widget* wgt, void* idx);
	static void cb_cb_dispLSMT(Fl_Widget *wgt, void *idx);
	static void cb_cb_dispOFST(Fl_Widget *wgt, void *idx);

	// ------------------------- CONTROL -------------------------------------------
public:
	Fl_Value_Input* vi_axis_y;
	Fl_Value_Input* vi_axis_z;
	Fl_Value_Slider* vs_cpidx;
	Fl_Value_Slider* vs_cpangv;

private:
	static void cb_vi_axis(Fl_Widget* wgt, void* idx);
	static void cb_vs_cpidx(Fl_Widget* wgt, void* idx);
	static void cb_vs_cpangv(Fl_Widget* wgt, void* idx);

};

#endif	// CPANEL

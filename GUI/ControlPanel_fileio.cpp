#define _USE_MATH_DEFINES
#include <math.h>
#include <sys\stat.h>
#include <direct.h>
#include <time.h>

#include "ControlPanel.h"
#include "util.h"

void ControlPanel::cb_btn_load(Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;
	GraphWindowDraw3D *gwin = This->gwin;
	Fl_File_Chooser *fc = This->fc;
	char wgtlbl[256];	strcpy(wgtlbl, wgt->label());
	int count;

	if(!strcmp(wgtlbl,"load"))
	{
		fc->show();

		while (fc->visible())
			Fl::wait();

		count = fc->count();

		if( count > 0 )
		{
			int sts = 0;
			char fname[1024], fname_m2m3[1024];
			strcpy(fname, fc->value());

			if( strstr( fname, ".txt") )
			{
				// load fname
				sts = This->sd->load(fname);
				This->vs_cpidx->bounds(0, This->sd->cpcnt - 1);
				This->vs_cpangv->value( This->sd->cpangv[(int)This->vs_cpidx->value()] );
				gwin->redraw();
			}
			if( sts==0 )
			{
				wgt->label("clr");
			}
		}
	}
	else if(!strcmp(wgtlbl,"clr"))
	{
		// clear data
		This->sd->clear();
		gwin->clear();
		gwin->initObject();
		gwin->redraw();
		wgt->label("load");
	}
}

void ControlPanel::cb_btn_savescreen( Fl_Widget *wgt, void *idx)
{
	ControlPanel *This = (ControlPanel *)idx;

	time_t t0 = time(NULL);
	struct tm *t1 = localtime(&t0);
	char fname[128];

	sprintf(fname, "%sscreen_%02d%02d%02d.bmp", DIR_OUT, t1->tm_hour, t1->tm_min, t1->tm_sec );
	This->gwin->exportImage( fname );
}

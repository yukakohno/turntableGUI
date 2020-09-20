#define _USE_MATH_DEFINES
#include <math.h>
#include "scandata.h"
#include "util.h"
#include "Spline.h"

scandata::scandata()
{
	linecnt = 0;
	memset(ptcnt, 0, sizeof(int) * MAX_LINE_CNT);
	cpcnt = 0;

	axis_y = 2.0;
	axis_z = 0.1;
}

scandata::~scandata()
{
}

int scandata::load(std::string fname)
{
	// for now
	int ret = makeSampleData();

	cpcnt = 9;
	for (int i = 0; i < cpcnt-1; i++) {
		cpidx[i] = (double)i * (double)linecnt / (double)(cpcnt-1);
	}
	cpidx[cpcnt-1] = linecnt-1;

	for (int i = 0; i < linecnt; i++) {
		lineangle[i] = 2.0 * M_PI * (double)i / (double)linecnt;
	}
	lineangv[0] = 0;
	for (int i = 1; i < linecnt; i++) {
		lineangv[i] = lineangle[i] - lineangle[i-1];
	}
	for (int i = 0; i < cpcnt; i++) {
		cpangv[i] = lineangv[cpidx[i]];
	}
	return ret;
}

int scandata::makeSampleData()
{
	double sample_axis_y = 2.0;
	linecnt = 1000;
	for (int i = 0; i < linecnt; i++) {
		ptcnt[i] = 500;
		for (int j = 0; j < ptcnt[i]; j++) {
			double a = M_PI * ((double)j / (double)ptcnt[i] -0.5);
			x[i][j] = sin( a );
			y[i][j] = sample_axis_y - cos( a );
		}
	}

	return 0;
}

void scandata::clear()
{
	linecnt = 0;
	memset(ptcnt, 0, sizeof(int) * MAX_LINE_CNT);
	memset(ptcnt, 0, sizeof(int) * MAX_LINE_CNT);
}

int scandata::cpangv2lineangle()
{
	int ret = 0;

	Spline sp;
	sp.init(cpangv, cpcnt);
	for (int i = 0; i < linecnt; i++) {
		double t = (double)i / (double)linecnt * (double)(cpcnt - 1);
		lineangv[i] = sp.culc(t);
	}

	lineangle[0] = 0;
	for (int i = 1; i < linecnt; i++) {
		lineangle[i] = lineangle[i - 1] + lineangv[i];
	}

	return ret;
}

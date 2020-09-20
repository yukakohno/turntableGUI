/*************************************************************
	scandata.h
*************************************************************/

#ifndef SCANDATA
#define SCANDATA

#include <string>

#define MAX_LINE_CNT 1000
#define MAX_PT_CNT 1000
#define MAX_CP_CNT 10

class scandata
{
public:
	scandata();
	~scandata();

	int load(std::string fname);
	int makeSampleData();
	void clear();
	int cpangv2lineangle();

	double axis_y, axis_z;

	int linecnt;
	int ptcnt[MAX_LINE_CNT];
	double x[MAX_LINE_CNT][MAX_PT_CNT];
	double y[MAX_LINE_CNT][MAX_PT_CNT];
	double lineangle[MAX_LINE_CNT];
	double lineangv[MAX_LINE_CNT];
	int cpcnt;
	int cpidx[MAX_CP_CNT];
	//double cpangle[MAX_CP_CNT];
	double cpangv[MAX_CP_CNT];

private:

};

#endif	// SCANDATA

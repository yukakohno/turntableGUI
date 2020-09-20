/***************************************
	Quatanian.cpp
***************************************/

#include "Quaternion.h"
#include <math.h>

#define RENORMCOUNT 97
#define TRACKBALLSIZE (0.8f)

CQuaternion::CQuaternion()
{

}

CQuaternion::~CQuaternion()
{

}

void CQuaternion::vzero(float *v)
{
	v[0] = 0.0f;
	v[1] = 0.0f;
	v[2] = 0.0f;
}

void CQuaternion::vset(float *v, float x, float y, float z)
{
	v[0] = x;
	v[1] = y;
	v[2] = z;
}

void CQuaternion::vsub(const float *src1, const float *src2, float *dst)
{
	dst[0] = src1[0] - src2[0];
	dst[1] = src1[1] - src2[1];
	dst[2] = src1[2] - src2[2];
}

void CQuaternion::vcopy(const float *v1, float *v2)
{
	register int i;
	for (i=0; i<3; i++)
		v2[i] = v1[i];
}

void CQuaternion::vcross(const float *v1, const float *v2, float *cross)
{
	float temp[3];

	temp[0] = (v1[1] * v2[2]) - (v1[2] * v2[1]);
	temp[1] = (v1[2] * v2[0]) - (v1[0] * v2[2]);
	temp[2] = (v1[0] * v2[1]) - (v1[1] * v2[0]);
	vcopy(temp, cross);
}

float CQuaternion::vlength(const float *v)
{
	return (float)sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

void CQuaternion::vscale(float *v, float div)
{
	v[0] *= div;
	v[1] *= div;
	v[2] *= div;
}

void CQuaternion::vnormal(float *v)
{
	vscale(v, 1.0f / vlength(v));
}

float CQuaternion::vdot(const float *v1, const float *v2)
{
	return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

void CQuaternion::vadd(const float *src1, const float *src2, float *dst)
{
	dst[0] = src1[0] + src2[0];
	dst[1] = src1[1] + src2[1];
	dst[2] = src1[2] + src2[2];
}

void CQuaternion::trackball(float q[4], float p1x, float p1y, float p2x, float p2y)
{
	float a[3];		/* Axis of rotation */
	float phi;		/* How much to rotate about axis. */
	float p1[3], p2[3], d[3];
	float t;

	if (p1x == p2x && p1y == p2y) {
		/* Zero rotation */
		vzero(q);
		q[3] = 1.0;
		return;
	}

	vset(p1, p1x, p1y, tb_project_to_sphere(TRACKBALLSIZE, p1x, p1y));
	vset(p2, p2x, p2y, tb_project_to_sphere(TRACKBALLSIZE, p2x, p2y));

	vcross(p2, p1, a);

	vsub(p1, p2, d);
	t = vlength(d) / (2.0f * TRACKBALLSIZE);

	/* Avoid problems with out-of-control values */
	if (t>1.0)
		t = 1.0;
	if (t < -1.0)
		t = -1.0;
	phi = 2.0f * (float)asin(t);

	axis_to_quat(a, phi, q);
}

/* Given an axis and angle, compute quatenion */
void CQuaternion::axis_to_quat(float a[3], float phi, float q[4])
{
	vnormal(a);
	vcopy(a, q);
	vscale(q, (float)sin(phi / 2.0f));
	q[3] = (float)cos(phi / 2.0f);
}

float CQuaternion::tb_project_to_sphere(float r, float x, float y)
{
	float d, t, z;

	d = (float)sqrt(x * x + y * y);
	if (d < r * 0.70710678118654752440) { /* Inside sphere */
		z = (float)sqrt(r * r - d * d);
	} else { /* On hyperbola. */
		t = r / 1.41421356237309504880f;
		z = t * t / d;
	}
	return z;
}


void CQuaternion::add_quat(float q1[4], float q2[4], float dest[4])
{
	static int count = 0;
	float t1[4], t2[4], t3[4];
	float tf[4];

	vcopy(q1, t1);
	vscale(t1, q2[3]);

	vcopy(q2, t2);
	vscale(t2, q1[3]);

	vcross(q2, q1, t3);
	vadd(t1, t2, tf);
	vadd(t3, tf, tf);
	tf[3] = q1[3] * q2[3] - vdot(q1, q2);

	dest[0] = tf[0];
	dest[1] = tf[1];
	dest[2] = tf[2];
	dest[3] = tf[3];

	//if (++count > RENORMCOUNT) 
	{
		count = 0;
		normalize_quat(dest);
	}
}

void CQuaternion::normalize_quat(float q[4])
{
	int i;
	float mag;

	//mag = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
	mag = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
	for (i=0; i<4; i++)
		q[i] /= mag;
}

/* Build a rotation matrix, given a quatenion rotation. */
void CQuaternion::build_rotmatrix(float m[4][4], float q[4])
{
	m[0][0] = 1.0f - 2.0f * (q[1] * q[1] + q[2] * q[2]);
	m[0][1] = 2.0f * (q[0] * q[1] - q[2] * q[3]);
	m[0][2] = 2.0f * (q[2] * q[0] + q[1] * q[3]);
	m[0][3] = 0.0f;

	m[1][0] = 2.0f * (q[0] * q[1] + q[2] * q[3]);
	m[1][1] = 1.0f - 2.0f * (q[2] * q[2] + q[0] * q[0]);
	m[1][2] = 2.0f * (q[1] * q[2] - q[0] * q[3]);
	m[1][3] = 0.0f;

	m[2][0] = 2.0f * (q[2] * q[0] - q[1] * q[3]);
	m[2][1] = 2.0f * (q[1] * q[2] + q[0] * q[3]);
	m[2][2] = 1.0f - 2.0f * (q[1] * q[1] + q[0] * q[0]);
	m[2][3] = 0.0f;

	m[3][0] = 0.0f;
	m[3][1] = 0.0f;
	m[3][2] = 0.0f;	
	m[3][3] = 1.0f;
}

void CQuaternion::build_rotmatrix(float m[16], float q[4])
{
	m[0] = 1.0f - 2.0f * (q[1] * q[1] + q[2] * q[2]);
	m[1] = 2.0f * (q[0] * q[1] - q[2] * q[3]);
	m[2] = 2.0f * (q[2] * q[0] + q[1] * q[3]);
	m[3] = 0.0f;

	m[4] = 2.0f * (q[0] * q[1] + q[2] * q[3]);
	m[5] = 1.0f - 2.0f * (q[2] * q[2] + q[0] * q[0]);
	m[6] = 2.0f * (q[1] * q[2] - q[0] * q[3]);
	m[7] = 0.0f;

	m[8] = 2.0f * (q[2] * q[0] - q[1] * q[3]);
	m[9] = 2.0f * (q[1] * q[2] + q[0] * q[3]);
	m[10] = 1.0f - 2.0f * (q[1] * q[1] + q[0] * q[0]);
	m[11] = 0.0f;

	m[12] = 0.0f;
	m[13] = 0.0f;
	m[14] = 0.0f;	
	m[15] = 1.0f;
}

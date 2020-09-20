#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>

//#include <GL/glu.h>
#include <GL/glut.h>
//#include <glut.h>
#include <FL/Fl.h>
#include <FL/gl.h>
#include <FL/Enumerations.H>
#include <FL/fl_draw.H>
#include "GraphWindowDraw3D.h"
#include "ControlPanel.h"
#include "util.h"

/* OpenCV */
#include <opencv2/opencv.hpp>
//using namespace std;
//using namespace cv;

GraphWindowDraw3D::GraphWindowDraw3D(int X, int Y, int W, int H, const char *L)
: GraphWindow3D(X, Y, W, H, L)
{
	GraphWindowDraw3D(X, Y, W, H);
}

GraphWindowDraw3D::GraphWindowDraw3D(int X, int Y, int W, int H)
: GraphWindow3D(X, Y, W, H)
{
	init();
}

GraphWindowDraw3D::~GraphWindowDraw3D()
{
	init();
}

void GraphWindowDraw3D::initObject()
{
	prevquat2[0] = prevquat2[1] = prevquat2[2] = 0; prevquat2[3] = 1;
	currquat2[0] = currquat2[1] = currquat2[2] = 0; currquat2[3] = 1;
	dispquat2[0] = dispquat2[1] = dispquat2[2] = 0; dispquat2[3] = 1;
	prevtrans2[0] = prevtrans2[1] = prevtrans2[2] = 0; prevtrans2[3] = 1;
	currtrans2[0] = currtrans2[1] = currtrans2[2] = 0; currtrans2[3] = 1;
	disptrans2[0] = disptrans2[1] = disptrans2[2] = 0; disptrans2[3] = 1;
	unit_m44( mObject );
}

void GraphWindowDraw3D::initTexture()
{
	/* テクスチャの読み込みに使う配列 */
	GLubyte texture[TEXHEIGHT][TEXWIDTH][4];

	/* テクスチャ画像の読み込み */
	cv::Mat tex = cv::imread(TEXFNAME, 1);
	for(int j=0 ; j<tex.rows; j++){
		for(int i=0 ; i<tex.cols; i++){
			texture[j][i][0] = (GLubyte)((cv::Vec3b &)tex.at<cv::Vec3b>(j,i))[0];
			texture[j][i][1] = (GLubyte)((cv::Vec3b &)tex.at<cv::Vec3b>(j,i))[1];
			texture[j][i][2] = (GLubyte)((cv::Vec3b &)tex.at<cv::Vec3b>(j,i))[2];
			texture[j][i][3] = (GLubyte)255;
		}
	}

	/* テクスチャ画像はバイト単位に詰め込まれている */
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

	/* テクスチャを拡大・縮小する方法の指定 */
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

	/* テクスチャの繰り返し方法の指定 */
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

	/* テクスチャ環境 */
	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	/* テクスチャの割り当て */
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, TEXWIDTH, TEXHEIGHT,
		0, GL_RGBA, GL_UNSIGNED_BYTE, texture);

	/* アルファテストの判別関数 */
	//glAlphaFunc(GL_GREATER, 0.5);
}

void GraphWindowDraw3D::init()
{
	GraphWindow3D::init();

	disp_axis = 1;
	disp_axis2 = 1;
	disp_LIN_SMOOTH = 0;
	disp_POLY_OFFSET = 0;

	push_button2 = 0;
	push_x2 = push_y2 = 0;
	initObject();

	resetRotation(-1);

	sd = NULL;
}

void GraphWindowDraw3D::resetObjTrans()
{
	prevtrans2[0] = prevtrans2[1] = prevtrans2[2] = 0; prevtrans2[3] = 1;
	currtrans2[0] = currtrans2[1] = currtrans2[2] = 0; currtrans2[3] = 1;
	disptrans2[0] = disptrans2[1] = disptrans2[2] = 0; disptrans2[3] = 1;
	mObject[12] = disptrans2[0];
	mObject[13] = disptrans2[1];
	mObject[14] = disptrans2[2];
	//redraw();
}

void GraphWindowDraw3D::resetObjRot()
{
	quat2.trackball(prevquat2, 0.0, 0.0, 0.0, 0.0);
	quat2.trackball(currquat2, 0.0, 0.0, 0.0, 0.0);
	quat2.trackball(dispquat2, 0.0, 0.0, 0.0, 0.0); 
	quat2.build_rotmatrix(mObject, dispquat2);
	mObject[12] = disptrans2[0];
	mObject[13] = disptrans2[1];
	mObject[14] = disptrans2[2];
	//redraw();
}

extern double norm( double *x, double *y, double *z );

void GraphWindowDraw3D::draw3DCurveFold()
{
	/* Draw 3D objects with OpenGL */

	// set camera position & rotation 
	// which is same as to set position & rotation of the whole world
	glPushMatrix();
	glTranslated(0.,0.,-n_plane-1000);
	SetRTS();

	if (disp_axis) {
		glScaled(0.1, 0.1, 0.1);
		drawAxis();
	}

	// set object position and rotation
	//glPushMatrix();
	//SetRTS2();

	// lrf origin & axis
	if( disp_axis2 )
	{
		glPushMatrix();
		glTranslatef(0.0, -sd->axis_y, -sd->axis_z);
		glScaled( 0.1, 0.1, 0.1 );
		drawAxis();
		glPopMatrix();
	}

	glPushMatrix();
	if (sd) {
		glPushMatrix();
		//glScaled(100.0, 100.0, 100.0);
		for (int i = 0; i < sd->linecnt; i++) {
			float r, g, b;
			float sec = (float)i / (float)sd->linecnt;
			glPointSize(1.0);

			for (int j = 0; j < sd->cpcnt; j++) {
				if (i == sd->cpidx[j]) {
					sec = -1;
					break;
				}
			}
			if (sec < 0) {
				r = g = b = 0.0;
				glPointSize(2.0);
			}
			else if (sec < 0.33) {
				r = (0.33-sec)/0.33;
				g = sec/0.33;
				b = 0.0;
			}
			else if (sec < 0.66) {
				r = 0.0;
				g = (0.66-sec)/0.33;
				b = (sec-0.33)/0.33;
			}
			else {
				r = (sec-0.66)/0.34;
				g = 0.0;
				b = (1.0-sec)/0.34;
			}
			glColor3f(r, g, b);
			glPushMatrix();
			glRotated(sd->lineangle[i]*180.0/M_PI, 1.0, 0.0, 0.0);
			glTranslatef(0.0, -sd->axis_y, -sd->axis_z);
			glBegin(GL_POINTS);
			for (int j = 0; j < sd->ptcnt[i]; j++) {
				glVertex3f(sd->x[i][j], sd->y[i][j], 0.0);
			}
			glEnd();
			glPopMatrix();
		}
		glPopMatrix();
	}
	//glPopMatrix();
	glPopMatrix();
}

void GraphWindowDraw3D::SetRTS2()
{
	glMultMatrixf(mObject);
}

void GraphWindowDraw3D::draw()
{
	if (!valid())
	{
		initTexture();
		SetLightMat();	// 光源環境設定

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();

		if(projection)
			glFrustum(-this->w()/2, this->w()/2, -this->h()/2, this->h()/2, n_plane, f_plane);// l,r, b,t, n,f 
		else
			glOrtho(-this->w()/2, this->w()/2, -this->h()/2, this->h()/2, n_plane, f_plane);

		glClearColor (1.f, 1.f, 1.f, 1.f);
	}
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_TEXTURE);
	glLoadIdentity();
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glViewport (0, 0, this->w(), this->h()); // 書かなくてもOK

	if( disp_LIN_SMOOTH ){
		glDepthFunc(GL_LEQUAL);
		glEnable(GL_LINE_SMOOTH);
		//glEnable(GL_POLYGON_SMOOTH);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		//glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	}

	draw3DCurveFold();

	glDisable(GL_LINE_SMOOTH);
	glDisable(GL_BLEND);

	glFlush();
	//	glutPostRedisplay();
}

int GraphWindowDraw3D::handle(int event)
{
	int i, key;

	switch(event){
	case FL_FOCUS:
		//return 1; // to detect FL_KEYDOWN & FL_KEYUP
		break;
	case FL_UNFOCUS:
		//return 1; // to detect FL_KEYDOWN & FL_KEYUP
		break;
	case FL_KEYDOWN:
#ifdef DEBUG_MOUSE_EVENT
		printf("FL_KEYDOWN\n");
#endif
		key = Fl::event_key();
		if(key=='r' || key=='R' ){ // reset rotation
			quat2.trackball(prevquat2, 0.0, 0.0, 0.0, 0.0);
			quat2.trackball(currquat2, 0.0, 0.0, 0.0, 0.0);
			quat2.trackball(dispquat2, 0.0, 0.0, 0.0, 0.0); 
			quat2.build_rotmatrix(mObject, dispquat2);
		} else if( key=='t' || key=='T' ){ // reset translation
			prevtrans2[0] = prevtrans2[1] = prevtrans2[2] = 0; prevtrans2[3] = 1;
			currtrans2[0] = currtrans2[1] = currtrans2[2] = 0; currtrans2[3] = 1;
			disptrans2[0] = disptrans2[1] = disptrans2[2] = 0; disptrans2[3] = 1;
		}
		redraw();
		break;
	case FL_KEYUP:
#ifdef DEBUG_MOUSE_EVENT
		printf("FL_KEYUP\n");
#endif
		break;
	case FL_PUSH:
		push_button2 = Fl::event_button();
		sts_sft2 = Fl::event_shift();
		sts_alt2 = Fl::event_alt();
		sts_ctrl2 = Fl::event_ctrl();
		push_x2 = Fl::event_x();
		push_y2 = Fl::event_y();

#ifdef DEBUG_MOUSE_EVENT
		printf("FL_PUSH\n");
		printf("sts_sft:%d, sts_alt:%d, sts_ctrl:%d, (%d,%d)\n",
			sts_sft, sts_alt, sts_ctrl, push_x, push_y);
#endif
		if(push_button2 == FL_RIGHT_MOUSE && (sts_sft2 || sts_alt2 || sts_ctrl2)){
			return 1; // ドラッグイベントを有効にする
		}else{
			return GraphWindow3D::handle(event);
			//return GraphWindow::handle(event);
		}
		break;
	case FL_DRAG:
#ifdef DEBUG_MOUSE_EVENT
		printf("FL_DRAG\n");
		printf("mouse dragged: (%d,%d)\n", Fl::event_x(), Fl::event_y());
#endif
		if( push_button2 == FL_RIGHT_MOUSE ){
			if(sts_sft2){	// rotation
				// ドラッグした量をcurrquatに保存
				quat2.trackball(currquat2, (2.0*push_x2 - w_w)/w_w, (w_h - 2.0*push_y2)/w_h,
					(2.0*Fl::event_x() - w_w)/w_w, (w_h - 2.0*Fl::event_y())/w_h);
				// prevquat（前回分）とcurrquat（今回分）を足してdispquat（表示用）に格納
				float tmpq0[4], tmpq1[4], camq[4], objq[4], curq[4];
				mat_quat( mCamera, camq );
				mat_quat( mObject, objq );
				camq[3] = -camq[3];
				objq[3] = -objq[3];
				quat2.add_quat(currquat2, camq, tmpq0);
				quat2.add_quat(tmpq0, objq, tmpq1);
				camq[3] = -camq[3];
				objq[3] = -objq[3];
				quat2.add_quat(camq, tmpq1, tmpq0);
				quat2.add_quat(objq, tmpq0, curq);
				quat2.add_quat(prevquat2, curq, dispquat2);

				quat2.build_rotmatrix(mObject, dispquat2);
				mObject[12] = disptrans2[0];
				mObject[13] = disptrans2[1];
				mObject[14] = disptrans2[2];

			} else if(sts_ctrl2){ // translation
				// マウス移動量 -> カメラ座標での物体移動量
				currtrans2[0] = Fl::event_x() - push_x2;
				currtrans2[1] = -(Fl::event_y() - push_y2);
				currtrans2[2] = 0;

				// 回転に応じて移動方向を修正
				mult_m44_v4(mCamera, currtrans2);

				// スケールに応じて移動量修正
				currtrans2[0] /= dispscale;
				currtrans2[1] /= dispscale;
				currtrans2[2] /= dispscale;

				// 前の移動量に加算
				disptrans2[0] = prevtrans2[0] + currtrans2[0];
				disptrans2[1] = prevtrans2[1] + currtrans2[1];
				disptrans2[2] = prevtrans2[2] + currtrans2[2];
				mObject[12] = disptrans2[0];
				mObject[13] = disptrans2[1];
				mObject[14] = disptrans2[2];
			}
			redraw();
		}
		break;
	case FL_RELEASE:
#ifdef DEBUG_MOUSE_EVENT
		printf("FL_RELEASE\n");
		printf("mouse released: (%d,%d)\n", Fl::event_x(), Fl::event_y());
#endif
		if( push_button2 == FL_RIGHT_MOUSE && (sts_sft2 || sts_alt2 || sts_ctrl2) )
		{
			if(sts_sft2){
				for(i=0;i<4;i++) prevquat2[i] = dispquat2[i];
			} else if(sts_ctrl2){
				for(i=0;i<3;i++) prevtrans2[i] = disptrans2[i];
			}

			redraw();

			// キー状態を解除
			sts_sft2 = 0;
			sts_ctrl2 = 0;
			push_button2 = 0;
			push_x2 = 0;
			push_y2 = 0;

			return 1;
		}

		// キー状態を解除
		sts_sft2 = 0;
		sts_ctrl2 = 0;
		push_button2 = 0;
		push_x2 = 0;
		push_y2 = 0;

		break;
	case FL_MOVE:
		break;
	default:
		break;
	}

	return GraphWindow3D::handle(event);
}

void GraphWindowDraw3D::exportImage(char *filename)
{
	int imgw, imgh, imgstep;
	unsigned char *img=NULL;

	draw(); // redraw immediately

	imgw = this->w();	// may need to be 4*X
	imgh = this->h();
	img = new unsigned char[imgw*imgh*3];
	imgstep = imgw*3;

	glReadPixels(0,0, imgw, imgh, GL_RGB, GL_UNSIGNED_BYTE, img);

	cv::Mat image( cv::Size(imgw, imgh), CV_8UC3 );
	for (int j = 0; j < imgh; j++) {
		for (int i = 0; i < imgw; i++) {
			image.at<cv::Vec3b>(j, i)[0] = img[(imgh-1-j)*imgstep + i*3+2];
			image.at<cv::Vec3b>(j, i)[1] = img[(imgh-1-j)*imgstep + i*3+1];
			image.at<cv::Vec3b>(j, i)[2] = img[(imgh-1-j)*imgstep + i*3];
		}
	}
	imwrite( filename, image);

	delete[] img; img = NULL;
}

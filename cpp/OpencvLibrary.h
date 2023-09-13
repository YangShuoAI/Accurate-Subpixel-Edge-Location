#include <stdio.h>
#include <tchar.h>
#include <iostream>
#include <fstream>

#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui_c.h>

#define CV_VERSION_ID CVAUX_STR(CV_MAJOR_VERSION) \
CVAUX_STR(CV_MINOR_VERSION) CVAUX_STR(CV_SUBMINOR_VERSION)

#ifdef _DEBUG
#define CVLIB(name) "opencv_" name CV_VERSION_ID "d"
#else
#define CVLIB(name) "opencv_" name CV_VERSION_ID
#endif

#pragma comment(lib, CVLIB("core"))
#pragma comment(lib, CVLIB("highgui"))
#pragma comment(lib, CVLIB("imgcodecs"))
#pragma comment(lib, CVLIB("imgproc"))
//#pragma comment(lib, CVLIB("alphamat"))
//#pragma comment(lib, CVLIB("aruco"))
//#pragma comment(lib, CVLIB("barcode"))
//#pragma comment(lib, CVLIB("bgsegm"))
//#pragma comment(lib, CVLIB("calib3d"))
//#pragma comment(lib, CVLIB("ccalib"))
//#pragma comment(lib, CVLIB("dnn"))
//#pragma comment(lib, CVLIB("features2d"))
//#pragma comment(lib, CVLIB("flann"))
//#pragma comment(lib, CVLIB("gapi"))
//#pragma comment(lib, CVLIB("ml"))
//#pragma comment(lib, CVLIB("objdetect"))
//#pragma comment(lib, CVLIB("photo"))
//#pragma comment(lib, CVLIB("rgbd"))
//#pragma comment(lib, CVLIB("saliency"))
//#pragma comment(lib, CVLIB("shape"))
//#pragma comment(lib, CVLIB("stitching"))
//#pragma comment(lib, CVLIB("surface_matching"))
//#pragma comment(lib, CVLIB("text"))
//#pragma comment(lib, CVLIB("video"))
//#pragma comment(lib, CVLIB("videoio"))
//#pragma comment(lib, CVLIB("videostab"))
//#pragma comment(lib, CVLIB("viz"))
//#pragma comment(lib, CVLIB("wechat_qrcode"))
//#pragma comment(lib, CVLIB("xfeatures2d"))
//#pragma comment(lib, CVLIB("ximgproc"))
//#pragma comment(lib, CVLIB("xobjdetect"))
//#pragma comment(lib, CVLIB("xphoto"))
//#pragma comment(lib, CVLIB("world"))





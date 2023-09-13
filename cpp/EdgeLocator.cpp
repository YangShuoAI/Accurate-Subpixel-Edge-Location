// Accurate subpixel edge location based on partial area effect
// YangShuo
// 2023-9-12
#include<stdio.h>
#include<direct.h>
#include<stdlib.h>
#include<io.h>
#include <minmax.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include "EdgeLocation.h"

#define TEST_RESULTS 1

#if TEST_RESULTS
////////////////////////////////////////////////////////////////////////////////////////////////////
// 注释：mkdirs
// mkdirs创建递归文件夹
//      dir      -I      文件夹路径
//      
////////////////////////////////////////////////////////////////////////////////////////////////////
int mkdirs(char* dir)
{
    int i, len;
    char str[100];
    strcpy(str, dir);//缓存文件路径
    len = (int)strlen(str);
    for (i = 0; i < len; i++)
    {
        if (str[i] == '\\')
        {
            str[i] = '\0';
            if (_access(str, 0) != 0)
            {
                int sts = _mkdir(str);
            }
            str[i] = '\\';
        }
    }
    if (len > 0 && _access(str, 0) != 0) //检测是否创建成功
    {
        int sts = _mkdir(str);
    }

    return 0;
}
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////
// 注释：Edge构造函数初始化
// EdgeLocator::Edge构造函数接口变量定义
//      cur_position      -I      图像一维索引值
//      cur_x             -I      图像亚像素X位置
//      cur_y             -I      图像亚像素Y位置
//      cur_nx            -I      沿x方向法向量
//      cur_ny            -I      沿y方向法向量
//      cur_curv          -I      曲率
//      cur_i0            -I      当前位置两侧强度
//      cur_i1            -I      当前位置两侧强度
//      
////////////////////////////////////////////////////////////////////////////////////////////////////
EdgeLocator::Edge::Edge(int    cur_position,
                        double cur_x,
                        double cur_y,
                        double cur_nx,
                        double cur_ny,
                        double cur_curv,
                        double cur_i0,
                        double cur_i1)
{
    position = cur_position;
    x = cur_x;
    y = cur_y;
    nx = cur_nx;
    ny = cur_ny;
    curv = cur_curv;
    i0 = cur_i0;
    i1 = cur_i1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// 注释：get 1D index inside image position
// EdgeLocator::Edge函数接口变量定义
//
////////////////////////////////////////////////////////////////////////////////////////////////////
int EdgeLocator::Edge::get_position()
{
    return position;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// 注释：get x
// EdgeLocator::Edge函数接口变量定义
//
////////////////////////////////////////////////////////////////////////////////////////////////////
double EdgeLocator::Edge::get_x()
{
    return x;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// 注释：get y
// EdgeLocator::Edge函数接口变量定义
//
////////////////////////////////////////////////////////////////////////////////////////////////////
double EdgeLocator::Edge::get_y()
{
    return y;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// 注释：get nx
// EdgeLocator::Edge函数接口变量定义
//
////////////////////////////////////////////////////////////////////////////////////////////////////
double EdgeLocator::Edge::get_nx()
{
    return nx;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// 注释：get ny
// EdgeLocator::Edge函数接口变量定义
//
////////////////////////////////////////////////////////////////////////////////////////////////////
double EdgeLocator::Edge::get_ny()
{
    return ny;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// 注释：get curv
// EdgeLocator::Edge函数接口变量定义
//
////////////////////////////////////////////////////////////////////////////////////////////////////
double EdgeLocator::Edge::get_curv()
{
    return curv;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// 注释：get i0
// EdgeLocator::Edge函数接口变量定义
//
////////////////////////////////////////////////////////////////////////////////////////////////////
double EdgeLocator::Edge::get_i0()
{
    return i0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// 注释：get i1
// EdgeLocator::Edge函数接口变量定义
//
////////////////////////////////////////////////////////////////////////////////////////////////////
double EdgeLocator::Edge::get_i1()
{
    return i1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// 注释：meshgrid初始化
// 创建meshgrid索引表
//      xgv               -I      xgv范围
//      ygv               -I      ygv范围
//      X                 -O      X索引表
//      Y                 -O      Y索引表
//
////////////////////////////////////////////////////////////////////////////////////////////////////
void  meshgrid(const cv::Range& xgv, const cv::Range& ygv, cv::Mat& X, cv::Mat& Y)
{
    std::vector<int> t_x, t_y;
    for (int i = xgv.start; i <= xgv.end; i++) t_x.push_back(i);
    for (int i = ygv.start; i <= ygv.end; i++) t_y.push_back(i);
    cv::Mat xgv_mat = cv::Mat(t_x);
    cv::Mat ygv_mat = cv::Mat(t_y);
    cv::repeat(xgv_mat.reshape(1, 1), (int)ygv_mat.total(), 1, X);
    cv::repeat(ygv_mat.reshape(1, 1).t(), 1, (int)xgv_mat.total(), Y);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// 注释：meshgrid2double初始化
// 创建meshgrid2double索引表(双精度浮点型)
//      xgv               -I      xgv范围
//      ygv               -I      ygv范围
//      X                 -O      X索引表
//      Y                 -O      Y索引表
//
////////////////////////////////////////////////////////////////////////////////////////////////////
void  meshgrid2double(const cv::Range& xgv, const cv::Range& ygv, cv::Mat& X, cv::Mat& Y)
{
    std::vector<double> t_x, t_y;
    for (int i = xgv.start; i <= xgv.end; i++) t_x.push_back(i);
    for (int i = ygv.start; i <= ygv.end; i++) t_y.push_back(i);
    cv::Mat xgv_mat = cv::Mat(t_x);
    cv::Mat ygv_mat = cv::Mat(t_y);
    cv::repeat(xgv_mat.reshape(1, 1), (int)ygv_mat.total(), 1, X);
    cv::repeat(ygv_mat.reshape(1, 1).t(), 1, (int)xgv_mat.total(), Y);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// 注释：circleGrid初始化
// 创建circleGrid函数
//      x               -I      x坐标向量
//      y               -I      y坐标向量
//      radius2         -I      圆半径
//      xCenter         -I      圆中心x坐标
//      yCenter         -I      圆中心x坐标
//      dx              -I      dx索引表
//      dy              -I      dy索引表
//      p               -O      p索引向量
//
////////////////////////////////////////////////////////////////////////////////////////////////////
void circleGrid(std::vector<double>  x,
                std::vector<double>  y,
                double               radius2,
                double               xCenter,
                double               yCenter,
                cv::Mat              dx,
                cv::Mat              dy,
                std::vector<double> &p)
{
    int numPixels = (int)x.size();
    double* dx_data = (double*)dx.data;
    double* dy_data = (double*)dy.data;

    if (numPixels > 0)
    {
        for (int n = 0; n < numPixels; n++)
        {
            cv::Mat grid_x = cv::Mat::zeros(dx.size(), CV_64FC1);
            cv::Mat grid_y = cv::Mat::zeros(dy.size(), CV_64FC1);
            cv::Mat grid = cv::Mat::zeros(dy.size(), CV_64FC1);
            grid_x = x[n] + dx - xCenter;
            grid_y = y[n] + dy - yCenter;
            grid.setTo(1, (grid_x.mul(grid_x) + (grid_y.mul(grid_y)) < radius2));
            p.push_back(cv::mean(grid)[0]);
        }
    }
    else
    {
        for (int n = 0; n < numPixels; n++)
        {
            p.push_back(0.0);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// 注释：circleVerticalWindow初始化
// 创建circleVerticalWindow函数
//      x               -I      x坐标向量
//      y               -I      y坐标向量
//      radius2         -I      圆半径
//      xCenter         -I      圆中心x坐标
//      yCenter         -I      圆中心x坐标
//      dx              -I      dx索引表
//      dy              -I      dy索引表
//      p               -O      p索引向量
//
////////////////////////////////////////////////////////////////////////////////////////////////////
void circleVerticalWindow(double  xWindowCenter,
                          double  yWindowCenter,
                          double  xCenter,
                          double  yCenter,
                          double  radius,
                          double  innerIntensity,
                          double  outerIntensity,
                          double  gridResolution,
                          cv::Mat &image)
{
    // computes a synthetic window of 3 columns with a circle inside
    // compute pixels completely outside or inside
    double r2 = radius * radius;
    cv::Mat x;
    cv::Mat y;
    cv::Range xgv = cv::Range((int)xWindowCenter - 1, (int)xWindowCenter + 1);
    cv::Range ygv = cv::Range((int)yWindowCenter - 4, (int)yWindowCenter + 4);
    meshgrid2double(xgv, ygv, x, y);

    cv::Mat c = cv::Mat::zeros(x.size(), CV_64FC1);
    cv::Mat cx1 = (x - 0.5 - xCenter);
    cv::Mat cy1 = (y - 0.5 - yCenter);
    cv::Mat cx2 = (x + 0.5 - xCenter);
    cv::Mat cy2 = (y + 0.5 - yCenter);

    c = ((cx1.mul(cx1) + cy1.mul(cy1)) < r2) / 255.0;
    c = c + ((cx1.mul(cx1) + cy2.mul(cy2)) < r2) / 255.0;
    c = c + ((cx2.mul(cx2) + cy1.mul(cy1)) < r2) / 255.0;
    c = c + ((cx2.mul(cx2) + cy2.mul(cy2)) < r2) / 255.0;
    c.convertTo(c, CV_64FC1);
    if (image.empty())
    {
        image = cv::Mat::zeros(x.size(), CV_64FC1);
    }
    c.copyTo(image);
    image.setTo(outerIntensity, c == 0);
    image.setTo(innerIntensity, c == 4);

    // compute contour pixels
    double delta = 1.0 / (gridResolution - 1);

    int dx_rows = (int)(ceil((0.5 - (-0.5)) / delta));
    int dx_cols = (int)(ceil((0.5 - (-0.5)) / delta));

    cv::Mat dx = cv::Mat::zeros(cv::Size(dx_rows, dx_cols), CV_64FC1);
    cv::Mat dy = cv::Mat::zeros(cv::Size(dx_rows, dx_cols), CV_64FC1);
    std::vector<double> d_x, d_y;
    double delta_step = (0.5 - (-0.5)) / (dx_rows * 1.0);
    for (double i = -0.5; i <= 0.5 + 1e-8; i = i + delta)
    {
        d_x.push_back(i);
    }
    for (double i = -0.5; i <= 0.5 + 1e-8; i = i + delta)
    {
        d_y.push_back(i);
    }
    cv::Mat xdx_mat = cv::Mat(d_x);
    cv::Mat ydy_mat = cv::Mat(d_y);
    cv::repeat(xdx_mat.reshape(1, 1), (int)ydy_mat.total(), 1, dx);
    cv::repeat(ydy_mat.reshape(1, 1).t(), 1, (int)xdx_mat.total(), dy);

    cv::Mat c_mat = cv::Mat::zeros(x.size(), CV_64FC1);
    c_mat.setTo(1, c > 0 & c < 4);
    cv::Mat c_mat_transpose;
    cv::Mat x_mat_transpose;
    cv::Mat y_mat_transpose;
    cv::transpose(c_mat, c_mat_transpose);
    cv::transpose(x, x_mat_transpose);
    cv::transpose(y, y_mat_transpose);
    cv::Mat x_flat = x.reshape(1, (int)(x.total() * x.channels()));
    std::vector<double> x_mat_vec = x.isContinuous() ? x_flat : x_flat.clone();
    cv::Mat y_flat = y.reshape(1, (int)(y.total() * y.channels()));
    std::vector<double> y_mat_vec = y.isContinuous() ? y_flat : y_flat.clone();
    std::vector<double> x_vec;
    std::vector<double> y_vec;

    for (int i = 0; i < c_mat_transpose.rows; i++)
    {
        for (int j = 0; j < c_mat_transpose.cols; j++)
        {
            int index = i * c_mat_transpose.cols + j;
            double c_mat_data = c_mat_transpose.at<double>(i, j);
            if (c_mat_data > 0)
            {
                x_vec.push_back(x_mat_transpose.at<double>(i, j));
                y_vec.push_back(y_mat_transpose.at<double>(i, j));
            }
        }
    }

    std::vector<double> p, i_val;
    circleGrid(x_vec, y_vec, r2, xCenter, yCenter, dx, dy, p);

    for (size_t i = 0; i < p.size(); i++)
    {
        i_val.push_back(outerIntensity + (innerIntensity - outerIntensity) * p[i]);
    }

    int i_val_index = 0;
    for (int i = 0; i < image.cols; i++)
    {
        for (int j = 0; j < image.rows; j++)
        {
            double c_mat_data = c_mat.at<double>(j, i);
            if (c_mat_data > 0)
            {
                image.at<double>(j, i) = i_val[i_val_index];
                i_val_index = i_val_index + 1;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// 注释：circleHorizontalWindow初始化
// 创建circleHorizontalWindow函数
//      x               -I      x坐标向量
//      y               -I      y坐标向量
//      radius          -I      圆半径
//      xCenter         -I      圆中心x坐标
//      yCenter         -I      圆中心x坐标
//      dx              -I      dx索引表
//      dy              -I      dy索引表
//      p               -O      p索引向量
//
////////////////////////////////////////////////////////////////////////////////////////////////////
void circleHorizontalWindow(double  xWindowCenter,
                            double  yWindowCenter,
                            double  xCenter,
                            double  yCenter,
                            double  radius,
                            double  innerIntensity,
                            double  outerIntensity,
                            double  gridResolution,
                            cv::Mat &image)
{
    // computes a synthetic window of 3 columns with a circle inside
    // compute pixels completely outside or inside
    double r2 = radius * radius;
    cv::Mat x;
    cv::Mat y;
    cv::Range xgv = cv::Range((int)xWindowCenter - 4, (int)xWindowCenter + 4);
    cv::Range ygv = cv::Range((int)yWindowCenter - 1, (int)yWindowCenter + 1);
    meshgrid2double(xgv, ygv, x, y);

    cv::Mat c = cv::Mat::zeros(x.size(), CV_64FC1);
    cv::Mat cx1 = (x - 0.5 - xCenter);
    cv::Mat cy1 = (y - 0.5 - yCenter);
    cv::Mat cx2 = (x + 0.5 - xCenter);
    cv::Mat cy2 = (y + 0.5 - yCenter);

    c = ((cx1.mul(cx1) + cy1.mul(cy1)) < r2) / 255.0;
    c = c + ((cx1.mul(cx1) + cy2.mul(cy2)) < r2) / 255.0;
    c = c + ((cx2.mul(cx2) + cy1.mul(cy1)) < r2) / 255.0;
    c = c + ((cx2.mul(cx2) + cy2.mul(cy2)) < r2) / 255.0;
    c.convertTo(c, CV_64FC1);
    if (image.empty())
    {
        image = cv::Mat::zeros(x.size(), CV_64FC1);
    }
    c.copyTo(image);
    image.setTo(outerIntensity, c == 0);
    image.setTo(innerIntensity, c == 4);

    // compute contour pixels
    double delta = 1.0 / (gridResolution - 1);

    int dx_rows = (int)(ceil((0.5 - (-0.5)) / delta));
    int dx_cols = (int)(ceil((0.5 - (-0.5)) / delta));

    cv::Mat dx = cv::Mat::zeros(cv::Size(dx_rows, dx_cols), CV_64FC1);
    cv::Mat dy = cv::Mat::zeros(cv::Size(dx_rows, dx_cols), CV_64FC1);
    std::vector<double> d_x, d_y;
    double delta_step = (0.5 - (-0.5)) / (dx_rows * 1.0);
    for (double i = -0.5; i <= 0.5 + 1e-8; i = i + delta)
    {
        d_x.push_back(i);
    }
    for (double i = -0.5; i <= 0.5 + 1e-8; i = i + delta)
    {
        d_y.push_back(i);
    }
    cv::Mat xdx_mat = cv::Mat(d_x);
    cv::Mat ydy_mat = cv::Mat(d_y);
    cv::repeat(xdx_mat.reshape(1, 1), (int)ydy_mat.total(), 1, dx);
    cv::repeat(ydy_mat.reshape(1, 1).t(), 1, (int)xdx_mat.total(), dy);

    cv::Mat c_mat = cv::Mat::zeros(x.size(), CV_64FC1);
    c_mat.setTo(1, c > 0 & c < 4);
    cv::Mat c_mat_transpose;
    cv::Mat x_mat_transpose;
    cv::Mat y_mat_transpose;
    cv::transpose(c_mat, c_mat_transpose);
    cv::transpose(x, x_mat_transpose);
    cv::transpose(y, y_mat_transpose);
    cv::Mat x_flat = x.reshape(1, (int)(x.total() * x.channels()));
    std::vector<double> x_mat_vec = x.isContinuous() ? x_flat : x_flat.clone();
    cv::Mat y_flat = y.reshape(1, (int)(y.total() * y.channels()));
    std::vector<double> y_mat_vec = y.isContinuous() ? y_flat : y_flat.clone();
    std::vector<double> x_vec;
    std::vector<double> y_vec;
    for (int i = 0; i < c_mat_transpose.rows; i++)
    {
        for (int j = 0; j < c_mat_transpose.cols; j++)
        {
            int index = i * c_mat_transpose.cols + j;
            double c_mat_data = c_mat_transpose.at<double>(i, j);
            if (c_mat_data > 0)
            {
                x_vec.push_back(x_mat_transpose.at<double>(i, j));
                y_vec.push_back(y_mat_transpose.at<double>(i, j));
            }
        }
    }

    std::vector<double> p, i_val;
    circleGrid(x_vec, y_vec, r2, xCenter, yCenter, dx, dy, p);

    for (size_t i = 0; i < p.size(); i++)
    {
        i_val.push_back(outerIntensity + (innerIntensity - outerIntensity) * p[i]);
    }

    int i_val_index = 0;
    for (int i = 0; i < image.cols; i++)
    {
        for (int j = 0; j < image.rows; j++)
        {
            double c_mat_data = c_mat.at<double>(j, i);
            if (c_mat_data > 0)
            {
                image.at<double>(j, i) = i_val[i_val_index];
                i_val_index = i_val_index + 1;
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// 注释：finalDetectorIterN内容
// 创建finalDetectorIterN函数
//      F               -I      输入图像
//      threshold       -I      输入阈值
//      order           -I      输入阶次
//      ep              -O      输出亚像素边缘信息
//      I               -O      输出边缘计算增强后图像
//
////////////////////////////////////////////////////////////////////////////////////////////////////
int EdgeLocator::finalDetectorIterN(cv::Mat            F,
                                    double             threshold,
                                    int                order,
                                    std::vector<Edge> &ep,
                                    cv::Mat&           I)
{
    // initialization
    int rows = F.rows;
    int cols = F.cols;
    int pixelGridResol = 50;
    std::vector<Edge> edgepixels;
    cv::Mat x;
    cv::Mat y;
    meshgrid(cv::Range(1, cols), cv::Range(1, rows), x, y);

    cv::Mat C = cv::Mat::zeros(x.size(), CV_64FC1);

    // smooth image
    double w = 0.75;// case for smooth mask 'average', H(1:3, 1 : 3) = 1 / 9
    F.convertTo(F, CV_64FC1);
    cv::Mat G = F.clone();
    cv::Rect roi_region(1, 1, cols - 2, rows - 2);
    cv::blur(G(roi_region), G(roi_region), cv::Size(3, 3), cv::Point(-1, -1));

    // compute partial derivatives
    cv::Mat Gx;
    cv::Mat Gy;
    cv::Mat grad;
    cv::Sobel(G, Gx, CV_64FC1, 1, 0, 1);
    cv::Sobel(G, Gy, CV_64FC1, 0, 1, 1);
    Gx = 0.5 * Gx;
    Gy = 0.5 * Gy;
    cv::sqrt(Gx.mul(Gx) + Gy.mul(Gy), grad);

    // detect edge pixels with maximum Gy (not including margins)
    cv::Rect gy_region(2, 5, cols - 2 - 2, rows - 5 - 5);
    cv::Mat absGyInner = cv::abs(Gy(gy_region).clone());
    std::vector<int> edges;
    cv::Mat E = cv::Mat::zeros(F.size(), CV_32SC1);
    cv::Mat F_transpose;
    cv::Mat G_transpose;
    cv::Mat Gx_transpose;
    cv::Mat Gy_transpose;
    cv::Mat x_transpose;
    cv::Mat y_transpose;
    cv::transpose(F, F_transpose);
    cv::transpose(G, G_transpose);
    cv::transpose(Gx, Gx_transpose);
    cv::transpose(Gy, Gy_transpose);
    cv::transpose(x, x_transpose);
    cv::transpose(y, y_transpose);
    double* F_data = (double*)F_transpose.data;
    double* G_data = (double*)G_transpose.data;
    double* Gx_data = (double*)Gx_transpose.data;
    double* Gy_data = (double*)Gy_transpose.data;
    int* x_data = (int*)x_transpose.data;
    int* y_data = (int*)y_transpose.data;
    //int *E_data = (int *)E.data;

    cv::Rect E_roi_region(2, 5, cols - 2 - 2, rows - 5 - 5);
    cv::Rect grad_roi_region(2, 5, cols - 2 - 2, rows - 5 - 5);
    cv::Rect Gx_roi_region(2, 5, cols - 2 - 2, rows - 5 - 5);
    cv::Rect Gy_roi_region1(2, 4, cols - 2 - 2, rows - 6 - 4);
    cv::Rect Gy_roi_region2(2, 6, cols - 2 - 2, rows - 4 - 6);

    cv::Mat E_roi = E(E_roi_region).clone();
    cv::Mat grad_roi = grad(grad_roi_region).clone();
    cv::Mat Gx_roi1 = Gx(Gx_roi_region).clone();
    cv::Mat Gy_roi1 = Gy(Gy_roi_region1).clone();
    cv::Mat Gy_roi2 = Gy(Gy_roi_region2).clone();

    int* E_roi_data = (int*)E_roi.data;
    double* absGyInner_data = (double*)absGyInner.data;
    double* grad_roi_data = (double*)grad_roi.data;
    double* Gx_roi1_data = (double*)Gx_roi1.data;
    double* Gy_roi1_data = (double*)Gy_roi1.data;
    double* Gy_roi2_data = (double*)Gy_roi2.data;

    int roi_rows = absGyInner.rows;
    int roi_cols = absGyInner.cols;
    for (int i = 0; i < roi_rows; i++)
    {
        for (int j = 0; j < roi_cols; j++)
        {
            int index = i * roi_cols + j;
            if ((grad_roi_data[index] > threshold)
                & (absGyInner_data[index] >= abs(Gx_roi1_data[index] + 1e-8))
                & (absGyInner_data[index] >= abs(Gy_roi1_data[index] + 1e-8))
                & (absGyInner_data[index] > abs(Gy_roi2_data[index] + 1e-8)))
            {
                E_roi_data[index] = 1;
            }

        }
    }
    E_roi.copyTo(E(E_roi_region));

    cv::Mat xE;
    cv::Mat yE;
    cv::multiply(E, x, xE);
    cv::multiply(E, y, yE);
    int* xE_data = (int*)xE.data;
    int* yE_data = (int*)yE.data;
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            int index = i * cols + j;
            if ((xE_data[index] > 0) && (yE_data[index] > 0))
            {
                int edge = (xE_data[index] - 1) * rows + yE_data[index] - 1;
                edges.push_back(edge);
            }
        }
    }
    std::sort(edges.begin(), edges.end());

    std::vector<double> A;
    std::vector<double> B;
    std::vector<double> a;
    std::vector<double> b;
    std::vector<double> c;
    std::vector<double> nx;
    std::vector<double> ny;
    std::vector<double> curv;
    std::vector<int> valid;

    A.resize(edges.size());
    B.resize(edges.size());
    a.resize(edges.size());
    b.resize(edges.size());
    c.resize(edges.size());
    nx.resize(edges.size());
    ny.resize(edges.size());
    curv.resize(edges.size());
    valid.resize(edges.size());

    // compute all horizontal edges
    for (int k = 0; k < edges.size(); k++)
    {
        // compute window floating limits
        int edge = edges[k];
        int m1 = -1;
        int m2 = 1;

        int m = 0;
        int l1 = 0;
        int r2 = 0;
        int minl1 = 0;
        int maxr2 = 0;
        int l2 = 0;
        int r1 = 0;
        int maxl2 = 0;
        int minr1 = 0;
        if ((Gx_data[edge] * Gy_data[edge]) >= 0.0)
        {
            m = 1;
            l1 = -1;
            r2 = 1;
            minl1 = -3;
            maxr2 = 3;
            l2 = 1;
            r1 = -1;
            maxl2 = 4;
            minr1 = -4;
        }
        else
        {
            m = -1;
            l1 = -1;
            r2 = 1;
            minl1 = -4;
            maxr2 = 4;
            l2 = 1;
            r1 = -1;
            maxl2 = 3;
            minr1 = -3;
        }

        while ((l1 > minl1) && (abs(Gy_data[edge - rows + l1] + 1e-6) >=
            abs(Gy_data[edge - rows + l1 - 1] + 1e-8)))
        {
            l1 = l1 - 1;
        }

        while ((l2 < maxl2) && (abs(Gy_data[edge - rows + l2] + 1e-6) >=
            abs(Gy_data[edge - rows + l2 + 1] + 1e-8)))
        {
            l2 = l2 + 1;
        }
        while ((m1 > -4) && (abs(Gy_data[edge + m1] + 1e-6) >=
            abs(Gy_data[edge + m1 - 1] + 1e-8)))
        {
            m1 = m1 - 1;
        }
        while ((m2 < 4.0) && (abs(Gy_data[edge + m2] + 1e-6) >=
            abs(Gy_data[edge + m2 + 1] + 1e-8)))
        {
            m2 = m2 + 1;
        }
        while ((r1 > minr1) && (abs(Gy_data[edge + rows + r1] + 1e-6) >=
            abs(Gy_data[edge + rows + r1 - 1] + 1e-8)))
        {
            r1 = r1 - 1;
        }
        while ((r2 < maxr2) && (abs(Gy_data[edge + rows + r2] + 1e-6) >=
            abs(Gy_data[edge + rows + r2 + 1] + 1e-8)))
        {
            r2 = r2 + 1;
        }

        // window
        cv::Mat window = cv::Mat::zeros(cv::Size(3, 9), CV_64FC1);
        window(cv::Rect(0, l1 + 5 - 1, 1, (l2 + 5 - l1 - 5 + 1))) = 1.0;
        window(cv::Rect(1, m1 + 5 - 1, 1, (m2 + 5 - m1 - 5 + 1))) = 100.0;
        window(cv::Rect(2, r1 + 5 - 1, 1, (r2 + 5 - r1 - 5 + 1))) = 1.0;

        double AA = 0.0;
        double BB = 0.0;
        // compute intensities
        if (m > 0)
        {
            AA = (G_data[edge + m2] + G_data[edge + rows + r2]) * 0.5;
            BB = (G_data[edge - rows + l1] + G_data[edge + m1]) * 0.5;
        }
        else
        {
            AA = (G_data[edge - rows + l2] + G_data[edge + m2]) * 0.5;
            BB = (G_data[edge + m1] + G_data[edge + rows + r1]) * 0.5;
        }

        // search for a second close edge
        bool uBorder = false;
        bool dBorder = false;
        double partial = 0.0;
        if (m1 > -4)
        {
            partial = abs(G_data[edge + m1 - 2] + 1e-8);
            if ((partial > abs(G_data[edge] * 0.25 + 1e-8)) && (partial > (threshold * 0.5)))
            {
                uBorder = true;
            }
        }
        if (m2 < 4)
        {
            partial = abs(G_data[edge + m2 + 2] + 1e-8);
            if ((partial > abs(G_data[edge] * 0.25 + 1e-8)) && (partial > (threshold * 0.5)))
            {
                dBorder = true;
            }
        }

        double SL = 0.0;
        double SM = 0.0;
        double SR = 0.0;

        int j = (int)(floor((edges[k] - 1.0) / (rows * 1.0))) + 1;
        int i = edges[k] + 1 - rows * (j - 1);
        int p = 0;
        int ll = 0;
        int rr = 0;

        if (uBorder || dBorder)
        {
            cv::Rect image_rect = cv::Rect(j - 3, i - 6, 5, 11);
            cv::Mat rimvt = F(image_rect).clone();
            if (uBorder)
            {
                if (m > 0)
                {
                    BB = (F_data[edge + m1] + F_data[edge - rows + l1]) * 0.5;
                    p = 1;
                }
                else
                {
                    BB = (F_data[edge + m1] + F_data[edge + rows + r1]) * 0.5;
                    p = 0;
                }

                if (Gy_data[edge - 2 * rows + l1 + p] * Gy_data[edge] > 0.0)
                {
                    ll = l1 + p - 1;
                }
                else
                {
                    ll = l1 + p;
                }
                if (Gy_data[edge + 2 * rows + r1 + 1 - p] * Gy_data[edge] > 0.0)
                {
                    rr = r1 - p;
                }
                else
                {
                    rr = r1 + 1 - p;
                }

                for (int xx = 0; xx < (ll + 6); xx++)
                {
                    rimvt.at<double>(xx, 0) = BB;
                }
                for (int xx = 0; xx < (l1 + 6); xx++)
                {
                    rimvt.at<double>(xx, 1) = BB;
                }
                for (int xx = 0; xx < (m1 + 6); xx++)
                {
                    rimvt.at<double>(xx, 2) = BB;
                }
                for (int xx = 0; xx < (r1 + 6); xx++)
                {
                    rimvt.at<double>(xx, 3) = BB;
                }
                for (int xx = 0; xx < (rr + 6); xx++)
                {
                    rimvt.at<double>(xx, 4) = BB;
                }

                l1 = -3 + m;
                m1 = -3;
                r1 = -3 - m;
            }

            if (dBorder)
            {
                if (m > 0)
                {
                    AA = (F_data[edge + m2] + F_data[edge + rows + r2]) * 0.5;
                    p = 1;
                }
                else
                {
                    AA = (F_data[edge + m2] + F_data[edge - rows + l2]) * 0.5;
                    p = 0;
                }

                if ((Gy_data[edge - 2 * rows + l2 + p - 1] * Gy_data[edge]) > 0.0)
                {
                    ll = l2 + p;
                }
                else
                {
                    ll = l2 + p - 1;
                }

                if (Gy_data[edge + 2 * rows + r2 - p] * Gy_data[edge] > 0.0)
                {
                    rr = r2 + 1 - p;
                }
                else
                {
                    rr = r2 - p;
                }

                for (int xx = ll + 5; xx < rimvt.rows; xx++)
                {
                    rimvt.at<double>(xx, 0) = AA;
                }
                for (int xx = l2 + 5; xx < rimvt.rows; xx++)
                {
                    rimvt.at<double>(xx, 1) = AA;
                }
                for (int xx = m2 + 5; xx < rimvt.rows; xx++)
                {
                    rimvt.at<double>(xx, 2) = AA;
                }
                for (int xx = r2 + 5; xx < rimvt.rows; xx++)
                {
                    rimvt.at<double>(xx, 3) = AA;
                }
                for (int xx = rr + 5; xx < rimvt.rows; xx++)
                {
                    rimvt.at<double>(xx, 4) = AA;
                }

                l2 = 3 + m;
                m2 = 3;
                r2 = 3 - m;
            }

            cv::Mat rimv2 = (rimvt(cv::Rect(0, 0, 3, 9)) + rimvt(cv::Rect(1, 0, 3, 9))
                + rimvt(cv::Rect(2, 0, 3, 9)) + rimvt(cv::Rect(0, 1, 3, 9))
                + rimvt(cv::Rect(1, 1, 3, 9)) + rimvt(cv::Rect(2, 1, 3, 9))
                + rimvt(cv::Rect(0, 2, 3, 9)) + rimvt(cv::Rect(1, 2, 3, 9))
                + rimvt(cv::Rect(2, 2, 3, 9))) / 9.0;

            for (int n = l1 + 4; n < l2 + 5; n++)
            {
                SL = SL + rimv2.at<double>(n, 0);
            }
            for (int n = m1 + 4; n < m2 + 5; n++)
            {
                SM = SM + rimv2.at<double>(n, 1);
            }
            for (int n = r1 + 4; n < r2 + 5; n++)
            {
                SR = SR + rimv2.at<double>(n, 2);
            }
        }
        else
        {
            for (int n = l1; n <= l2; n++)
            {
                SL = SL + G_data[edge - rows + n];
            }
            for (int n = m1; n <= m2; n++)
            {
                SM = SM + G_data[edge + n];
            }
            for (int n = r1; n <= r2; n++)
            {
                SR = SR + G_data[edge + rows + n];
            }
        }

        // compute edge features
        double den = 2.0 * (AA - BB);
        if (order == 2)
        {
            c[k] = (SL + SR - 2.0 * SM + AA * (2.0 * m2 - l2 - r2) -
                    BB * (2.0 * m1 - l1 - r1) + 1e-8) / den;
        }
        else
        {
            c[k] = 0.0;
        }
        a[k] = (2.0 * SM - AA * (1.0 + 2.0 * m2) -
                BB * (1.0 - 2.0 * m1) + 1e-8) / den - w * c[k];
        if (abs(a[k]) > 1.0)
        {
            valid[k] = 0;
            continue;
        }
        valid[k] = 1;
        b[k] = (SR - SL + AA * ((1.0 * l2 - 1.0 * r2)) -
                BB * ((1.0 * l1 - 1.0 * r1)) + 1e-8) / den;
        A[k] = AA;
        B[k] = BB;

        int sign = 0;
        double diff_AA_BB = AA - BB;
        if (diff_AA_BB > 0)
        {
            sign = 1;
        }
        else if (diff_AA_BB == 0.0)
        {
            sign = 0;
        }
        else
        {
            sign = -1;
        }

        nx[k] = sign / sqrt(1.0 + b[k] * b[k]) * b[k];
        ny[k] = sign / sqrt(1.0 + b[k] * b[k]);
        curv[k] = 2 * c[k] / ((1 + b[k] * b[k]) * sqrt((1.0 + b[k] * b[k])));
        if (Gy_data[edge] < 0.0)
        {
            curv[k] = -curv[k];
        }

        // generate circle sub image
        double R = abs(1.0 / curv[k]);
        if (R > 1e4)
        {
            R = 1e4;
        }
        else
        {
            if (R < 4.5)
            {
                R = 4.5;
            }
        }

        double innerIntensity = 0.0;
        double outerIntensity = 0.0;
        if (curv[k] > 0.0)
        {
            sign = -1;
            innerIntensity = std::min(AA, BB);
            outerIntensity = std::max(AA, BB);
        }
        else
        {
            sign = 1;
            innerIntensity = std::max(AA, BB);
            outerIntensity = std::min(AA, BB);
        }
        double cent_x = x_data[edge] + sign * R * nx[k];
        double cent_y = y_data[edge] - a[k] + sign * R * ny[k];
        cv::Point2d center = cv::Point2d(cent_x, cent_y);



        cv::Mat subimage;
        circleVerticalWindow(j, i, center.x, center.y,
            R, innerIntensity, outerIntensity,
            pixelGridResol, subimage);

        I(cv::Rect(j - 2, i - 5, 3, 9)) = I(cv::Rect(j - 2, i - 5, 3, 9)) + window.mul(subimage);
        C(cv::Rect(j - 2, i - 5, 3, 9)) = C(cv::Rect(j - 2, i - 5, 3, 9)) + window;
    }

    // edges = edges * valid
    std::transform(valid.begin(), valid.end(), edges.begin(), edges.begin(), std::multiplies<int>());
    // erase 0 elements in edges
    edges.erase(std::remove_if(edges.begin(), edges.end(), [](int n) { return  n == 0; }), edges.end());

    // A = A * valid
    std::transform(valid.begin(), valid.end(), A.begin(), A.begin(), std::multiplies<double>());
    // erase 0 elements in A
    A.erase(std::remove_if(A.begin(), A.end(), [](double n) { return  n == 0.0; }), A.end());

    // B = B * valid
    std::transform(valid.begin(), valid.end(), B.begin(), B.begin(), std::multiplies<double>());
    // erase 0 elements in B
    B.erase(std::remove_if(B.begin(), B.end(), [](double n) { return  n == 0.0; }), B.end());

    // a = a * valid
    std::transform(valid.begin(), valid.end(), a.begin(), a.begin(), std::multiplies<double>());
    // erase 0 elements in a
    a.erase(std::remove_if(a.begin(), a.end(), [](double n) { return  n == 0.0; }), a.end());

    // nx = nx * valid
    std::transform(valid.begin(), valid.end(), nx.begin(), nx.begin(), std::multiplies<double>());
    // erase 0 elements in nx
    nx.erase(std::remove_if(nx.begin(), nx.end(), [](double n) { return  n == 0.0; }), nx.end());

    // ny = ny * valid
    std::transform(valid.begin(), valid.end(), ny.begin(), ny.begin(), std::multiplies<double>());
    // erase 0 elements in ny
    ny.erase(std::remove_if(ny.begin(), ny.end(), [](double n) { return  n == 0.0; }), ny.end());

    // curv = curv * valid
    std::transform(valid.begin(), valid.end(), curv.begin(), curv.begin(), std::multiplies<double>());
    // erase 0 elements in curv
    curv.erase(std::remove_if(curv.begin(), curv.end(), [](double n) { return  n == 0.0; }), curv.end());

    for (size_t k = 0; k < edges.size(); k++)
    {
        int    cur_position = edges[k];
        double cur_x = x_data[cur_position];
        double cur_y = y_data[cur_position] - a[k];
        double cur_nx = nx[k];
        double cur_ny = ny[k];
        double cur_curv = curv[k];
        double cur_i0 = std::min(A[k], B[k]);
        double cur_i1 = std::max(A[k], B[k]);
        Edge edge(cur_position, cur_x, cur_y, cur_nx, cur_ny, cur_curv, cur_i0, cur_i1);
        ep.push_back(edge);
    }

    // detecte edge pixels with maximum Gx (not including margins)
    cv::Rect gx_region(5, 2, cols - 5 - 5, rows - 2 - 2);
    cv::Mat absGxInner = cv::abs(Gx(gx_region).clone());
    cv::Mat E2 = cv::Mat::zeros(F.size(), CV_32SC1);

    cv::Rect E_roi2_region(5, 2, cols - 5 - 5, rows - 2 - 2);
    cv::Rect grad_roi2_region(5, 2, cols - 5 - 5, rows - 2 - 2);
    cv::Rect Gy_roi2_region(5, 2, cols - 5 - 5, rows - 2 - 2);
    cv::Rect Gx_roi2_region1(4, 2, cols - 6 - 4, rows - 2 - 2);
    cv::Rect Gx_roi2_region2(6, 2, cols - 4 - 6, rows - 2 - 2);

    cv::Mat E_roi2 = E2(E_roi2_region).clone();
    cv::Mat grad_roi2 = grad(grad_roi2_region).clone();
    cv::Mat Gy_roi2_1 = Gy(Gy_roi2_region).clone();
    cv::Mat Gx_roi2_1 = Gx(Gx_roi2_region1).clone();
    cv::Mat Gx_roi2_2 = Gx(Gx_roi2_region2).clone();

    int* E_roi2_data = (int*)E_roi2.data;
    double* absGxInner_data = (double*)absGxInner.data;
    double* grad_roi2_data = (double*)grad_roi2.data;
    double* Gy_roi2_1_data = (double*)Gy_roi2_1.data;
    double* Gx_roi2_1_data = (double*)Gx_roi2_1.data;
    double* Gx_roi2_2_data = (double*)Gx_roi2_2.data;

    int roi2_rows = absGxInner.rows;
    int roi2_cols = absGxInner.cols;
    for (int i = 0; i < roi2_rows; i++)
    {
        for (int j = 0; j < roi2_cols; j++)
        {
            int index = i * roi2_cols + j;
            if ((grad_roi2_data[index] > threshold)
                & (absGxInner_data[index] >= abs(Gy_roi2_1_data[index] + 1e-8))
                & (absGxInner_data[index] >= abs(Gx_roi2_1_data[index] + 1e-8))
                & (absGxInner_data[index] > abs(Gx_roi2_2_data[index] + 1e-8)))
            {
                E_roi2_data[index] = 1;
            }

        }
    }
    E_roi2.copyTo(E2(E_roi2_region));

    cv::Mat xE2;
    cv::Mat yE2;
    cv::multiply(E2, x, xE2);
    cv::multiply(E2, y, yE2);
    int* xE2_data = (int*)xE2.data;
    int* yE2_data = (int*)yE2.data;

    edges.clear();
    A.clear();
    B.clear();
    a.clear();
    b.clear();
    c.clear();
    nx.clear();
    ny.clear();
    curv.clear();
    valid.clear();

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            int index = i * cols + j;
            if ((xE2_data[index] > 0) && (yE2_data[index] > 0))
            {
                int edge = (xE2_data[index] - 1) * rows + yE2_data[index] - 1;
                edges.push_back(edge);
            }
        }
    }
    std::sort(edges.begin(), edges.end());

    A.resize(edges.size());
    B.resize(edges.size());
    a.resize(edges.size());
    b.resize(edges.size());
    c.resize(edges.size());
    nx.resize(edges.size());
    ny.resize(edges.size());
    curv.resize(edges.size());
    valid.resize(edges.size());

    // compute all horizontal edges
    for (int k = 0; k < edges.size(); k++)
    {
        // compute window floating limits
        int edge = edges[k];
        int m1 = -1;
        int m2 = 1;

        int m = 0;
        int l1 = 0;
        int r2 = 0;
        int minl1 = 0;
        int maxr2 = 0;
        int l2 = 0;
        int r1 = 0;
        int maxl2 = 0;
        int minr1 = 0;
        if ((Gx_data[edge] * Gy_data[edge]) >= 0.0)
        {
            m = 1;
            l1 = -1;
            r2 = 1;
            minl1 = -3;
            maxr2 = 3;
            l2 = 1;
            r1 = -1;
            maxl2 = 4;
            minr1 = -4;
        }
        else
        {
            m = -1;
            l1 = -1;
            r2 = 1;
            minl1 = -4;
            maxr2 = 4;
            l2 = 1;
            r1 = -1;
            maxl2 = 3;
            minr1 = -3;
        }

        while ((l1 > minl1) && (abs(Gx_data[edge - 1 + l1 * rows] + 1e-6) >=
            abs(Gx_data[edge - 1 + (l1 - 1) * rows] + 1e-8)))
        {
            l1 = l1 - 1;
        }

        while ((l2 < maxl2) && (abs(Gx_data[edge - 1 + l2 * rows] + 1e-6) >=
            abs(Gx_data[edge - 1 + (l2 + 1) * rows] + 1e-8)))
        {
            l2 = l2 + 1;
        }

        while ((m1 > -4) && (abs(Gx_data[edge + m1 * rows] + 1e-6) >=
            abs(Gx_data[edge + (m1 - 1) * rows] + 1e-8)))
        {
            m1 = m1 - 1;
        }

        while ((m2 < 4.0) && (abs(Gx_data[edge + m2 * rows] + 1e-6) >=
            abs(Gx_data[edge + (m2 + 1) * rows] + 1e-8)))
        {
            m2 = m2 + 1;
        }

        while ((r1 > minr1) && (abs(Gx_data[edge + 1 + r1 * rows] + 1e-6) >=
            abs(Gx_data[edge + 1 + (r1 - 1) * rows] + 1e-8)))
        {
            r1 = r1 - 1;
        }

        while ((r2 < maxr2) && (abs(Gx_data[edge + 1 + r2 * rows] + 1e-6) >=
            abs(Gx_data[edge + 1 + (r2 + 1) * rows] + 1e-8)))
        {
            r2 = r2 + 1;
        }

        // window
        cv::Mat window = cv::Mat::zeros(cv::Size(9, 3), CV_64FC1);
        window(cv::Rect(l1 + 5 - 1, 0, (l2 + 5 - l1 - 5 + 1), 1)) = 1.0;
        window(cv::Rect(m1 + 5 - 1, 1, (m2 + 5 - m1 - 5 + 1), 1)) = 100.0;
        window(cv::Rect(r1 + 5 - 1, 2, (r2 + 5 - r1 - 5 + 1), 1)) = 1.0;

        double AA = 0.0;
        double BB = 0.0;
        // compute intensities
        if (m > 0)
        {
            AA = (G_data[edge + m2 * rows] + G_data[edge + 1 + r2 * rows]) * 0.5;
            BB = (G_data[edge - 1 + l1 * rows] + G_data[edge + m1 * rows]) * 0.5;
        }
        else
        {
            AA = (G_data[edge - 1 + l2 * rows] + G_data[edge + m2 * rows]) * 0.5;
            BB = (G_data[edge + m1 * rows] + G_data[edge + 1 + r1 * rows]) * 0.5;
        }

        // search for a second close edge
        bool uBorder = false;
        bool dBorder = false;
        double partial = 0.0;
        if (m1 > -4)
        {
            partial = abs(G_data[edge + (m1 - 2) * rows] + 1e-8);
            if ((partial > abs(G_data[edge] * 0.25 + 1e-8)) && (partial > (threshold * 0.5)))
            {
                uBorder = true;
            }
        }
        if (m2 < 4)
        {
            partial = abs(G_data[edge + (m2 + 2) * rows] + 1e-8);
            if ((partial > abs(G_data[edge] * 0.25 + 1e-8)) && (partial > (threshold * 0.5)))
            {
                dBorder = true;
            }
        }

        double SL = 0.0;
        double SM = 0.0;
        double SR = 0.0;

        int j = (int)(floor((edges[k] - 1.0) / (rows * 1.0))) + 1;
        int i = edges[k] + 1 - rows * (j - 1);
        int p = 0;
        int ll = 0;
        int rr = 0;

        if (uBorder || dBorder)
        {
            cv::Rect image_rect = cv::Rect(j - 6, i - 3, 11, 5);
            cv::Mat rimvt = F(image_rect).clone();
            if (uBorder)
            {
                if (m > 0)
                {
                    BB = (F_data[edge + m1 * rows] + F_data[edge - 1 + l1 * rows]) * 0.5;
                    p = 1;
                }
                else
                {
                    BB = (F_data[edge + m1 * rows] + F_data[edge + 1 + r1 * rows]) * 0.5;
                    p = 0;
                }

                if (Gx_data[edge - 2 + (l1 + p) * rows] * Gx_data[edge] > 0.0)
                {
                    ll = l1 + p - 1;
                }
                else
                {
                    ll = l1 + p;
                }

                if (Gx_data[edge + 2 + (r1 + 1 - p) * rows] * Gx_data[edge] > 0.0)
                {
                    rr = r1 - p;
                }
                else
                {
                    rr = r1 + 1 - p;
                }

                for (int xx = 0; xx < (ll + 6); xx++)
                {
                    rimvt.at<double>(0, xx) = BB;
                }
                for (int xx = 0; xx < (l1 + 6); xx++)
                {
                    rimvt.at<double>(1, xx) = BB;
                }
                for (int xx = 0; xx < (m1 + 6); xx++)
                {
                    rimvt.at<double>(2, xx) = BB;
                }
                for (int xx = 0; xx < (r1 + 6); xx++)
                {
                    rimvt.at<double>(3, xx) = BB;
                }
                for (int xx = 0; xx < (rr + 6); xx++)
                {
                    rimvt.at<double>(4, xx) = BB;
                }

                l1 = -3 + m;
                m1 = -3;
                r1 = -3 - m;
            }

            if (dBorder)
            {
                if (m > 0)
                {
                    AA = (F_data[edge + m2 * rows] + F_data[edge + 1 + r2 * rows]) * 0.5;
                    p = 1;
                }
                else
                {
                    AA = (F_data[edge + m2 * rows] + F_data[edge - 1 + l2 * rows]) * 0.5;
                    p = 0;
                }

                if ((Gx_data[edge - 2 + (l2 + p - 1) * rows] * Gx_data[edge]) > 0.0)
                {
                    ll = l2 + p;
                }
                else
                {
                    ll = l2 + p - 1;
                }

                if (Gx_data[edge + 2 + (r2 - p) * rows] * Gx_data[edge] > 0.0)
                {
                    rr = r2 + 1 - p;
                }
                else
                {
                    rr = r2 - p;
                }

                for (int xx = ll + 5; xx < rimvt.cols; xx++)
                {
                    rimvt.at<double>(0, xx) = AA;
                }
                for (int xx = l2 + 5; xx < rimvt.cols; xx++)
                {
                    rimvt.at<double>(1, xx) = AA;
                }
                for (int xx = m2 + 5; xx < rimvt.cols; xx++)
                {
                    rimvt.at<double>(2, xx) = AA;
                }
                for (int xx = r2 + 5; xx < rimvt.cols; xx++)
                {
                    rimvt.at<double>(3, xx) = AA;
                }
                for (int xx = rr + 5; xx < rimvt.cols; xx++)
                {
                    rimvt.at<double>(4, xx) = AA;
                }

                l2 = 3 + m;
                m2 = 3;
                r2 = 3 - m;
            }

            cv::Mat rimv2 = (rimvt(cv::Rect(0, 0, 9, 3)) + rimvt(cv::Rect(0, 1, 9, 3))
                + rimvt(cv::Rect(0, 2, 9, 3)) + rimvt(cv::Rect(1, 0, 9, 3))
                + rimvt(cv::Rect(1, 1, 9, 3)) + rimvt(cv::Rect(1, 2, 9, 3))
                + rimvt(cv::Rect(2, 0, 9, 3)) + rimvt(cv::Rect(2, 1, 9, 3))
                + rimvt(cv::Rect(2, 2, 9, 3))) / 9.0;

            for (int n = l1 + 4; n < l2 + 5; n++)
            {
                SL = SL + rimv2.at<double>(0, n);
            }
            for (int n = m1 + 4; n < m2 + 5; n++)
            {
                SM = SM + rimv2.at<double>(1, n);
            }
            for (int n = r1 + 4; n < r2 + 5; n++)
            {
                SR = SR + rimv2.at<double>(2, n);
            }
        }
        else
        {
            for (int n = l1; n <= l2; n++)
            {
                SL = SL + G_data[edge - 1 + n * rows];
            }
            for (int n = m1; n <= m2; n++)
            {
                SM = SM + G_data[edge + n * rows];
            }
            for (int n = r1; n <= r2; n++)
            {
                SR = SR + G_data[edge + 1 + n * rows];
            }
        }

        // compute edge features
        double den = 2.0 * (AA - BB);
        if (order == 2)
        {
            c[k] = (SL + SR - 2.0 * SM + AA * (2.0 * m2 - l2 - r2) -
                    BB * (2.0 * m1 - l1 - r1) + 1e-8) / den;
        }
        else
        {
            c[k] = 0.0;
        }
        a[k] = (2.0 * SM - AA * (1.0 + 2.0 * m2) -
                BB * (1.0 - 2.0 * m1) + 1e-8) / den - w * c[k];
        if (abs(a[k]) > 1.0)
        {
            valid[k] = 0;
            continue;
        }
        valid[k] = 1;
        // 乘以1.0是将int转换为double
        b[k] = (SR - SL + AA * ((1.0 * l2 - 1.0 * r2)) -
                BB * ((1.0 * l1 - 1.0 * r1)) + 1e-8) / den;
        A[k] = AA;
        B[k] = BB;

        int sign = 0;
        double diff_AA_BB = AA - BB;
        if (diff_AA_BB > 0)
        {
            sign = 1;
        }
        else if (diff_AA_BB == 0.0)
        {
            sign = 0;
        }
        else
        {
            sign = -1;
        }

        nx[k] = sign / sqrt(1.0 + b[k] * b[k]);
        ny[k] = sign / sqrt(1.0 + b[k] * b[k]) * b[k];
        curv[k] = 2 * c[k] / ((1 + b[k] * b[k]) * sqrt((1.0 + b[k] * b[k])));
        if (Gx_data[edge] < 0.0)
        {
            curv[k] = -curv[k];
        }

        // generate circle sub image
        double R = abs(1.0 / curv[k]);
        if (R > 1e4)
        {
            R = 1e4;
        }
        else
        {
            if (R < 4.5)
            {
                R = 4.5;
            }
        }

        double innerIntensity = 0.0;
        double outerIntensity = 0.0;
        if (curv[k] > 0.0)
        {
            sign = -1;
            innerIntensity = std::min(AA, BB);
            outerIntensity = std::max(AA, BB);
        }
        else
        {
            sign = 1;
            innerIntensity = std::max(AA, BB);
            outerIntensity = std::min(AA, BB);
        }
        double cent_x = x_data[edge] - a[k] + sign * R * nx[k];
        double cent_y = y_data[edge] + sign * R * ny[k];
        cv::Point2d center = cv::Point2d(cent_x, cent_y);

        cv::Mat subimage;
        circleHorizontalWindow(j, i, center.x, center.y,
            R, innerIntensity, outerIntensity,
            pixelGridResol, subimage);

        I(cv::Rect(j - 5, i - 2, 9, 3)) = I(cv::Rect(j - 5, i - 2, 9, 3)) + window.mul(subimage);
        C(cv::Rect(j - 5, i - 2, 9, 3)) = C(cv::Rect(j - 5, i - 2, 9, 3)) + window;
    }

    // edges = edges * valid
    std::transform(valid.begin(), valid.end(), edges.begin(), edges.begin(), std::multiplies<int>());
    // erase 0 elements in edges
    edges.erase(std::remove_if(edges.begin(), edges.end(), [](int n) { return  n == 0; }), edges.end());

    // A = A * valid
    std::transform(valid.begin(), valid.end(), A.begin(), A.begin(), std::multiplies<double>());
    // erase 0 elements in A
    A.erase(std::remove_if(A.begin(), A.end(), [](double n) { return  n == 0.0; }), A.end());

    // B = B * valid
    std::transform(valid.begin(), valid.end(), B.begin(), B.begin(), std::multiplies<double>());
    // erase 0 elements in B
    B.erase(std::remove_if(B.begin(), B.end(), [](double n) { return  n == 0.0; }), B.end());

    // a = a * valid
    std::transform(valid.begin(), valid.end(), a.begin(), a.begin(), std::multiplies<double>());
    // erase 0 elements in a
    a.erase(std::remove_if(a.begin(), a.end(), [](double n) { return  n == 0.0; }), a.end());

    // nx = nx * valid
    std::transform(valid.begin(), valid.end(), nx.begin(), nx.begin(), std::multiplies<double>());
    // erase 0 elements in nx
    nx.erase(std::remove_if(nx.begin(), nx.end(), [](double n) { return  n == 0.0; }), nx.end());

    // ny = ny * valid
    std::transform(valid.begin(), valid.end(), ny.begin(), ny.begin(), std::multiplies<double>());
    // erase 0 elements in ny
    ny.erase(std::remove_if(ny.begin(), ny.end(), [](double n) { return  n == 0.0; }), ny.end());

    // curv = curv * valid
    std::transform(valid.begin(), valid.end(), curv.begin(), curv.begin(), std::multiplies<double>());
    // erase 0 elements in curv
    curv.erase(std::remove_if(curv.begin(), curv.end(), [](double n) { return  n == 0.0; }), curv.end());

    for (size_t k = 0; k < edges.size(); k++)
    {
        int    cur_position = edges[k];
        double cur_x = x_data[cur_position] - a[k];
        double cur_y = y_data[cur_position];
        double cur_nx = nx[k];
        double cur_ny = ny[k];
        double cur_curv = curv[k];
        double cur_i0 = std::min(A[k], B[k]);
        double cur_i1 = std::max(A[k], B[k]);
        Edge edge(cur_position, cur_x, cur_y, cur_nx, cur_ny, cur_curv, cur_i0, cur_i1);
        ep.push_back(edge);
    }

    // compute final subimage
    cv::Mat I1 = cv::Mat::zeros(I.size(), I.type());
    cv::Mat C1 = cv::Mat::ones(C.size(), C.type());
    cv::Mat G1 = cv::Mat::zeros(G.size(), G.type());
    I.copyTo(I1, C > 0.0);
    C.copyTo(C1, C > 0.0);
    cv::Mat II = I1 / C1;
    II.copyTo(I, C > 0.0);

    G.copyTo(G1, C == 0.0);
    G1.copyTo(I, C == 0.0);

#if TEST_RESULTS
    FILE* fp;
    char savepath[128] = "./results/";
    mkdirs(savepath);
    fp = fopen("./results/position.txt", "w+ ");
    for (int iii = 0; iii < ep.size(); iii++)
    {
        fprintf(fp, "%d\n", ep[iii].get_position() + 1);
    }
    fclose(fp);

    fp = fopen("./results/x.txt", "w+ ");
    for (int iii = 0; iii < ep.size(); iii++)
    {
        fprintf(fp, "%.4f\n", ep[iii].get_x());
    }
    fclose(fp);

    fp = fopen("./results/y.txt", "w+ ");
    for (int iii = 0; iii < ep.size(); iii++)
    {
        fprintf(fp, "%.4f\n", ep[iii].get_y());
    }
    fclose(fp);

    fp = fopen("./results/nx.txt", "w+ ");
    for (int iii = 0; iii < ep.size(); iii++)
    {
        fprintf(fp, "%.4f\n", ep[iii].get_nx());
    }
    fclose(fp);

    fp = fopen("./results/ny.txt", "w+ ");
    for (int iii = 0; iii < ep.size(); iii++)
    {
        fprintf(fp, "%.4f\n", ep[iii].get_ny());
    }
    fclose(fp);

    fp = fopen("./results/curv.txt", "w+ ");
    for (int iii = 0; iii < ep.size(); iii++)
    {
        fprintf(fp, "%.4f\n", ep[iii].get_curv());
    }
    fclose(fp);

    fp = fopen("./results/i0.txt", "w+ ");
    for (int iii = 0; iii < ep.size(); iii++)
    {
        fprintf(fp, "%.4f\n", ep[iii].get_i0());
    }
    fclose(fp);

    //FILE* fp;
    fp = fopen("./results/i1.txt", "w+ ");
    for (int iii = 0; iii < ep.size(); iii++)
    {
        fprintf(fp, "%.4f\n", ep[iii].get_i1());
    }
    fclose(fp);

    //FILE* fp;
    fp = fopen("./results/I.txt", "w + ");
    for (int iii = 0; iii < I.rows; iii++)
    {
        for (int jjj = 0; jjj < I.cols; jjj++)
        {
            //int index = iii * I.cols + jjj;
            fprintf(fp, "%.4f ", I.at<double>(iii, jjj));
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
#endif

    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// 注释：subpixelEdges
// 创建subpixelEdges函数
//      image           -I      输入图像
//      threshold       -I      输入阈值
//      iter            -I      输入迭代次数
//      edges           -O      输出亚像素边缘信息
//
////////////////////////////////////////////////////////////////////////////////////////////////////
int EdgeLocator::subpixelEdges(cv::Mat image, double threshold, int iter, std::vector<Edge>& edges)
{
    int order = 2;
    if (image.channels() > 1)
    {
        cv::cvtColor(image, image, cv::COLOR_BGR2GRAY);
    }

    // smoothing iteration
    int smoothIter = iter;
    for (int i = 0; i < smoothIter; i++)
    {
        cv::Mat I = cv::Mat::zeros(image.size(), CV_64FC1);
        int sts = finalDetectorIterN(image, threshold, order, edges, I);
        edges.clear();
        I.copyTo(image);
    }

    return 0;
}





// Accurate subpixel edge location based on partial area effect
// YangShuo
// 2023-9-12
#ifndef _EDGE_LOCATION_H_
#define _EDGE_LOCATION_H_

#include <opencv2/opencv.hpp>

namespace EdgeLocator
{
    class Edge
    {
    private:
        int position;   // 1D index inside image
        double x;       // sub pixel position x
        double y;       // sub pixel position y
        double nx;      // normal vector x(normalized)
        double ny;      // normal vector y(normalized)
        double curv;    // curvature
        double i0;      // intensities at both sides
        double i1;      // intensities at both sides

    public:
        // EdgeLocator::Edge Constructor
        Edge(int position, double x, double y,
             double nx, double ny, double curv,
             double i0, double i1);
        // get 1D index inside image position
        int get_position();
        // get x
        double get_x();
        // get y
        double get_y();
        // get nx
        double get_nx();
        // get ny
        double get_ny();
        // get curv
        double get_curv();
        // get i0
        double get_i0();
        // get i1
        double get_i1();
    };

    // subpixelEdges(image, threshold, 'SmoothingIter', iter);
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // 注释：subpixelEdges
    // 创建subpixelEdges函数
    //      image           -I      输入图像
    //      threshold       -I      输入阈值
    //      iter            -I      输入迭代次数
    //      edges           -O      输出亚像素边缘信息
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    int subpixelEdges(cv::Mat image, double threshold, int iter, std::vector<Edge>& edges);

    // finalDetectorIterN
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
    int finalDetectorIterN(cv::Mat            F,
                           double             threshold,
                           int                order,
                           std::vector<Edge> &ep,
                           cv::Mat           &I);
}
#endif


// Accurate subpixel edge location based on partial area effect
// YangShuo
// 2023-9-12
#include "OpencvLibrary.h"
#include "EdgeLocation.h"
int main()
{
    cv::Mat src = cv::imread("angio2.png", cv::IMREAD_ANYCOLOR | cv::IMREAD_ANYDEPTH);
    cv::Mat gray;
    cv::cvtColor(src, gray, cv::COLOR_BGR2GRAY);

    double threshold = 4.0 + 1e-8;
    int iter = 3;

    // [edges, RI] = subpixelEdges(image, threshold, 'SmoothingIter', iter);
    std::vector<EdgeLocator::Edge> edges;
    int sts = EdgeLocator::subpixelEdges(gray, threshold, iter, edges);

    return 0;
}









# Accurate-Subpixel-Edge-Location
Accurate Subpixel Edge Location Based on Partial Area Effect
Abstract：  
&emsp;&emsp;The estimation of edge features, such as subpixel position, orientation, curvature and change in intensity at
both sides of the edge, from the computation of the gradient vector in each pixel is usually inexact, even in
ideal images. In this paper, we present a new edge detector based on an edge and acquisition model derived
from the partial area effect, which does not assume continuity in the image values. The main goal of this
method consists in achieving a highly accurate extraction of the position, orientation, curvature and contrast
of the edges, even in difficult conditions, such as noisy images, blurred edges, low contrast areas or very close
contours. For this purpose, we first analyze the influence of perfectly straight or circular edges in the surrounding
region, in such a way that, when these conditions are fulfilled, the features can exactly be determined.
Afterward, we extend it to more realistic situations considering how adverse conditions can be
tackled and presenting an iterative scheme for improving the results. We have tested this method in real
as well as in sets of synthetic images with extremely difficult edges, and in both cases a highly accurate characterization
has been achieved.

# 摘要：
基于局部区域效应的精确亚像素边缘定位方法  
&emsp;&emsp;摘要：从每个像素的梯度矢量来计算边缘特征，比如亚像素位置、方向、曲率和边缘两侧强度的变化等，即使是在理想图像下也通常是不准确的。在本文中，我们提出了一种新的边缘检测方法，该方法是由局部区域效应推导出的边缘采集模型，并且假设其在图像值上是不连续的。该方法的主要目标是实现即使在有噪声，边缘模糊，对比度低或者双边缘接近等困难情况下，也能对边缘的位置、方向、曲率和边缘对比度进行高精确提取。为此，我们首先分析了理想直线或圆形区域周边边缘的影响，当这些条件都满足时，边缘特征便可以准确确定。随后，我们把它扩展到更实际的情况，考虑如何在不利条件下进行计算，并提出了一种改善结果的迭代方案。我们已经用该方法测试了在实际场景以及用困难边缘合成的图像，并且在这两种情况下都实现了对边缘特征的高度精确描述。

C++ implementation of 《Accurate Subpixel Edge Location Based on Partial Area Effect》  
《基于局部区域效应的精确亚像素边缘定位方法》的C++源码实现

C++ vs Matlab
![image](https://github.com/YangShuoAI/Accurate-Subpixel-Edge-Location/assets/5794094/f85eaf72-08fe-402d-8011-749ebe581d76)

![2023-09-13 193241](https://github.com/YangShuoAI/Accurate-Subpixel-Edge-Location/assets/5794094/f871946f-9f97-4564-9f28-297dc034ad55)

![2023-09-13 193312](https://github.com/YangShuoAI/Accurate-Subpixel-Edge-Location/assets/5794094/c716ab73-665e-4db1-ae3c-e5bae74801f6)

![2023-09-13 190454](https://github.com/YangShuoAI/Accurate-Subpixel-Edge-Location/assets/5794094/b9896a59-a85d-4962-b144-82fddec9bfc5)

![2023-09-13 190641](https://github.com/YangShuoAI/Accurate-Subpixel-Edge-Location/assets/5794094/93f59515-73b1-4d72-9518-b3575c2523a4)

![2023-09-13 184150](https://github.com/YangShuoAI/Accurate-Subpixel-Edge-Location/assets/5794094/2c45c8f0-e1fe-4a06-9e98-ae8c8b6ce428)

![2023-09-13 184219](https://github.com/YangShuoAI/Accurate-Subpixel-Edge-Location/assets/5794094/65227d11-018d-4064-8537-7481122ba597)

![2023-09-13 185648](https://github.com/YangShuoAI/Accurate-Subpixel-Edge-Location/assets/5794094/b24d3f3a-cd04-4089-8cc5-bb4181a21586)

![2023-09-13 185714](https://github.com/YangShuoAI/Accurate-Subpixel-Edge-Location/assets/5794094/d50c494a-3345-4126-a8fb-d1fb5ea304b0)







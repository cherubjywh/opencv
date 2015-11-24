//
//  main.cpp
//  image_warping
//
//  Created by Yuyin Sun on 15-11-23.
//  Copyright (c) 2015年 Yuyin Sun. All rights reserved.
//

#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

using namespace cv;
using namespace std;

int main(int argc, const char * argv[]) {
    
    Point2f srcTri[3];
    Point2f dstTri[3];
    

    
    
    Mat im = imread("/Users/sunyuyin/Desktop/test_image.jpg");
    
    if (im.empty()) {
        exit(-1);
    }
    
    namedWindow("Image");
    imshow("Image", im);
    
    Mat& src = im;
    
    srcTri[0] = Point2f(0, 0);
    srcTri[1] = Point2f(src.cols - 1, 0);
    srcTri[2] = Point2f(0, src.rows - 1);
    
    dstTri[0] = Point2f(src.cols * 0.0, src.rows*0.33);
    dstTri[1] = Point2f(src.cols * 0.85, src.rows * 0.25);
    dstTri[2] = Point2f(src.cols * 0.15, src.rows * 0.7);
    
    Mat warp_mat = getAffineTransform(srcTri, dstTri);
    
    
    Mat dst = Mat::zeros(src.rows, src.cols, src.type());
    
    warpAffine(src, dst, warp_mat, dst.size());
    

    Point center = Point(dst.cols/2, dst.rows/2);
    double angle = -50.0;
    double scale = 0.6;
    
    Mat rot_mat = getRotationMatrix2D(center, angle, scale);
    
    warpAffine(dst, dst, rot_mat, dst.size());
 
    namedWindow("Warped Image");
    imshow("Warped Image", dst);
    
    
    waitKey(0);
    
    destroyAllWindows();
    
    return 0;
}

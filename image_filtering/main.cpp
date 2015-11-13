//
//  main.cpp
//  image_filtering
//
//  Created by Yuyin Sun on 15-11-12.
//  Copyright (c) 2015年 Yuyin Sun. All rights reserved.
//

#include <iostream>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv/cv.h>
#include <opencv/highgui.h>

using namespace std;
using namespace cv;

bool erode_image_opencv() {
    IplImage* img = cvLoadImage("/Users/sunyuyin/Desktop/test_image.jpg");
    cvNamedWindow("My Image");
    cvShowImage("My Image", img);
    
    cvErode(img, img, 0, 2);
    cvNamedWindow("Enroded");
    cvShowImage("Enroded", img);
    
    cvWaitKey(0);
    
    cvDestroyWindow("My Image");
    cvDestroyWindow("Enroded");
    cvReleaseImage(&img);
    
    return true;
}

bool bilateral_image_opencv2() {
    Mat src = imread("/Users/sunyuyin/Desktop/test_image.jpg", CV_LOAD_IMAGE_UNCHANGED);
    if (src.empty()) {
        return false;
    }
    
    namedWindow("My Image", WINDOW_AUTOSIZE);
    imshow("My Image", src);
    
    Mat dst(src.rows, src.cols, src.type());
    
    bilateralFilter(src, dst, 5, 10.0, 1.0);
    
    namedWindow("Bilateral Filter", WINDOW_AUTOSIZE);
    imshow("Bilateral Filter", dst);
    
    waitKey(0);
    destroyAllWindows();
    
    return true;
}

bool blur_image_opencv2() {
    
    Mat src = imread("/Users/sunyuyin/Desktop/test_image.jpg", CV_LOAD_IMAGE_UNCHANGED);
    if (src.empty()) {
        return false;
    }
    
    namedWindow("My Image", WINDOW_AUTOSIZE);
    imshow("My Image", src);
    
    Mat dst(src.rows, src.cols, src.type());
    
    blur(src, dst, Size(5, 5));
    
    namedWindow("Blurred Image", WINDOW_AUTOSIZE);
    imshow("Blurred Image", dst);
    
    waitKey(0);
    destroyAllWindows();

    return true;
}

bool call_box_filter() {
    
    Mat src = imread("/Users/sunyuyin/Desktop/test_image.jpg", CV_LOAD_IMAGE_UNCHANGED);
    if (src.empty()) {
        return false;
    }
    
    namedWindow("My Image", WINDOW_AUTOSIZE);
    imshow("My Image", src);
    
    Mat dst(src.rows, src.cols, src.type());
    
    boxFilter(src, dst, -1, Size(15, 15));
    namedWindow("Box Filter Image", WINDOW_AUTOSIZE);
    imshow("Box Fileter Image", dst);
    
    waitKey(0);
    destroyAllWindows();
    
    return true;
}

bool call_build_pyramid() {
    Mat src = imread("/Users/sunyuyin/Desktop/test_image.jpg", CV_LOAD_IMAGE_UNCHANGED);
    if (src.empty()) {
        return false;
    }
    
    namedWindow("My Image", WINDOW_AUTOSIZE);
    imshow("My Image", src);
    vector<Mat> dst;
    
    buildPyramid(src, dst, 4);
    
    namedWindow("Pyramid", WINDOW_AUTOSIZE);
    
    for (auto &i : dst) {
        imshow("Pyramid", i);
        waitKey(0);
    }
    
    destroyAllWindows();
    
    return true;
}




bool call_copy_make_border() {
    
    Mat src = imread("/Users/sunyuyin/Desktop/test_image.jpg", CV_LOAD_IMAGE_UNCHANGED);
    if (src.empty()) {
        return false;
    }
    
    namedWindow("My Image", WINDOW_AUTOSIZE);
    imshow("My Image", src);

    int top = 10, bottom = 20, left = 30, right = 40;
    Mat dst(src.cols + left + right, src.rows + top + bottom, src.type());
    
    copyMakeBorder(src, dst, top, bottom, left, right, BORDER_WRAP);
    
    namedWindow("Copy make border", WINDOW_AUTOSIZE);
    imshow("Copy make border", dst);
    
    waitKey(0);
    destroyAllWindows();
    
    return true;
}

bool call_gaussian_filter() {
    
    Mat src = imread("/Users/sunyuyin/Desktop/test_image.jpg", CV_LOAD_IMAGE_UNCHANGED);
    if (src.empty()) {
        return false;
    }
    
    namedWindow("My Image", WINDOW_AUTOSIZE);
    imshow("My Image", src);
    
    Mat dst(src.rows, src.cols, src.type());
    
    GaussianBlur(src, dst, Size(15, 15), 1);
    namedWindow("Gaussian Filter Image", WINDOW_AUTOSIZE);
    imshow("Gaussian Fileter Image", dst);
    
    waitKey(0);
    destroyAllWindows();
    
    return true;
}

bool call_erode() {
    
    Mat src = imread("/Users/sunyuyin/Desktop/test_image.jpg", CV_LOAD_IMAGE_UNCHANGED);
    if (src.empty()) {
        return false;
    }
    
    namedWindow("My Image", WINDOW_AUTOSIZE);
    imshow("My Image", src);
    
    Mat dst(src.rows, src.cols, src.type());
    
    erode(src, dst, Mat(5, 5, CV_8UC1, Scalar(1)));
    
    namedWindow("erode Image", WINDOW_AUTOSIZE);
    imshow("erode Image", dst);
    
    waitKey(0);
    destroyAllWindows();
    
    return true;
}

int main(int argc, const char * argv[]) {
    
//    bilateral_image_opencv2();
//    blur_image_opencv2();
//    call_box_filter();
//    call_build_pyramid();
//    call_copy_make_border();
//    call_gaussian_filter();
    call_erode();
    
    return 0;
}

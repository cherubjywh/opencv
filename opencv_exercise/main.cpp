//
//  main.cpp
//  opencv_exercise
//
//  Created by Yuyin Sun on 15-11-11.
//  Copyright (c) 2015年 Yuyin Sun. All rights reserved.
//

// Need to change c++ standard lib to libstdc++

#include <iostream>
#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"

using namespace std;
using namespace cv;

int main(int argc, const char * argv[]) {

    
    Mat src = imread("/Users/sunyuyin/Desktop/test_img.png", 1);
    
    namedWindow("Unprocessed Image", WINDOW_AUTOSIZE);
    imshow("Unprocessed Image", src);
    
    Mat dst = src.clone();
    GaussianBlur(src, dst, Size(15, 15), 0, 0);
    
    namedWindow("Processed Image", WINDOW_AUTOSIZE);
    imshow("Processed Image", dst);
    
    waitKey();
    
    return 0;
}

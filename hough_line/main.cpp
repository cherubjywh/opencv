//
//  main.cpp
//  hough_line
//
//  Created by Yuyin Sun on 15-11-14.
//  Copyright (c) 2015年 Yuyin Sun. All rights reserved.
//

#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace std;
using namespace cv;

int main(int argc, const char * argv[]) {

    const string filename = "/Users/sunyuyin/Desktop/test_image_2.jpeg";
    Mat src = imread(filename);
    
    if (src.empty()) {
        CV_Error(CV_StsBadArg, "Cannot read image");
    }
    
    
    namedWindow("Image");
    imshow("Image", src);
    
    Mat mask;
    Canny(src, mask, 100, 200, 3);

    
    
    Mat dst;
    cvtColor(mask, dst, COLOR_GRAY2BGR);

    namedWindow("Canny");
    imshow("Canny", dst);
    
    
    vector<Vec4i> lines;
    HoughLinesP(mask, lines, 1, CV_PI/180, 50, 60, 5);
    
    for (size_t i = 0; i < lines.size(); ++i) {
        auto l = lines[i];
        line(dst, Point(l[0], l[1]), Point(l[2], l[3]), Scalar(0, 0, 255), 3);
    }
    
    namedWindow("Lines");
    imshow("Lines", dst);
    
    waitKey(0);
    destroyAllWindows();
    
    return 0;
}

//
//  main.cpp
//  read_display_image
//
//  Created by Yuyin Sun on 15-11-11.
//  Copyright (c) 2015年 Yuyin Sun. All rights reserved.
//

#include <iostream>
#include "opencv2/highgui/highgui.hpp"

using namespace std;
using namespace cv;


void read_show_image() {
    Mat src = imread("/Users/sunyuyin/Documents/Workspace/CPP/opencv_exercise/data/orl_faces/s1/1.pgm", CV_LOAD_IMAGE_UNCHANGED);
    
    if (src.empty()) {
        cout << "Error: Image cannot be loaded..!!" << endl;
        return;
    }
    
    namedWindow("MyWindow", CV_WINDOW_KEEPRATIO);
    imshow("MyWindow", src);
    
    waitKey(0);
    
    destroyWindow("MyWindow");
}

void generate_show_blank_image() {
    Mat src(100, 200, CV_8UC3, Scalar(0, 0, 100));
    
    namedWindow("MyWindow", CV_WINDOW_AUTOSIZE);
    imshow("MyWindow", src);
    
    waitKey(0);
    
    destroyWindow("MyWindow");
}

int main(int argc, const char * argv[]) {
    
    read_show_image();
//    generate_show_blank_image();
    
    
    return 0;
}

//
//  main.cpp
//  write_image_video
//
//  Created by Yuyin Sun on 15-11-12.
//  Copyright (c) 2015年 Yuyin Sun. All rights reserved.
//

#include <iostream>
#include <opencv2/highgui/highgui.hpp>

using namespace std;
using namespace cv;


void write_video_to_file () {
    VideoCapture cap(0);
    
    if (!cap.isOpened()) {
        cout << "Error: Cannot open the video file" << endl;
        return;
    }
    
    namedWindow("My Video", CV_WINDOW_FULLSCREEN);
    
    double dWidth = cap.get(CV_CAP_PROP_FRAME_WIDTH);
    double dHeight = cap.get(CV_CAP_PROP_FRAME_HEIGHT);
    
    cout << "Frame Size = " << dWidth << "x" << dHeight << endl;
    
    Size frameSize(static_cast<int>(dWidth), static_cast<int>(dHeight));
    
    
}

void write_image_to_file () {
    Mat src(100, 200, CV_8UC3, Scalar(0, 0, 100));
    
    vector<int> compression_params = {CV_IMWRITE_JPEG_QUALITY, 98};
    
    bool bSuccess = imwrite("/Users/sunyuyin/Desktop/test_output.jpg", src, compression_params);
    
    if (!bSuccess) {
        cout << "Error: Fail to save the image" << endl;
    }
    
    namedWindow("My Image", CV_WINDOW_FULLSCREEN);
    imshow("My Image", src);
    
    waitKey(0);
    
    destroyWindow("My Image");
}

int main(int argc, const char * argv[]) {
    
    
    return 0;
}

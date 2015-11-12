//
//  main.cpp
//  capture_video
//
//  Created by Yuyin Sun on 15-11-11.
//  Copyright (c) 2015年 Yuyin Sun. All rights reserved.
//

#include <iostream>

#include "opencv2/highgui/highgui.hpp"

using namespace std;
using namespace cv;


void use_camera() {
    
    // Due to opencv version problem, this code does not work on Mac.
    // I will try to fix it later, got fixed
    
    VideoCapture cap(-1);
    
    if (!cap.isOpened())
    {
        cout << "Cannot open the video cam" << endl;
        return;
    }
    
    
    double dWidth = cap.get(CV_CAP_PROP_FRAME_WIDTH);
    double dHeight = cap.get(CV_CAP_PROP_FRAME_HEIGHT);
    
    cout << "Frame size: " << dWidth << "x" << dHeight << endl;
    
    namedWindow("MyVideo", CV_WINDOW_FULLSCREEN);
    
    
    while (1) {
        Mat frame;
        
        
        bool bSuccess = cap.read(frame);
        
        if (!bSuccess)
        {
            cout << "Cannot read a frame from video stream" << endl;
            break;
        }
        
        imshow("MyVideo", frame);
        
        if (waitKey(30) == 27)
        {
            cout << "esc key is pressed by user" << endl;
            break;
        }
    }
    
    destroyWindow("MyVideo");
}

void read_mp4() {
    VideoCapture cap("/Users/sunyuyin/Desktop/ICRA16_1354_VI_i.mp4"); // open the video file for reading
    
    if ( !cap.isOpened() )  // if not success, exit program
    {
        cout << "Cannot open the video file" << endl;
        return;
    }
    
    //cap.set(CV_CAP_PROP_POS_MSEC, 300); //start the video at 300ms
    
    double fps = cap.get(CV_CAP_PROP_FPS); //get the frames per seconds of the video
    
    cout << "Frame per seconds : " << fps << endl;
    
    namedWindow("MyVideo",CV_WINDOW_AUTOSIZE); //create a window called "MyVideo"
    
    while(1)
    {
        Mat frame;
        
        bool bSuccess = cap.read(frame); // read a new frame from video
        
        if (!bSuccess) //if not success, break loop
        {
            cout << "Cannot read the frame from video file" << endl;
            break;
        }
        
        imshow("MyVideo", frame); //show the frame in "MyVideo" window
        
        if(waitKey(30) == 27) //wait for 'esc' key press for 30 ms. If 'esc' key is pressed, break loop
        {
            cout << "esc key is pressed by user" << endl;
            break;
        }
    }
}

int main(int argc, const char * argv[]) {
    

//    read_mp4();
    
    use_camera();
    
    return 0;
}

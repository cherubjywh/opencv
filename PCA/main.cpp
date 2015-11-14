//
//  main.cpp
//  PCA
//
//  Created by Yuyin Sun on 15-11-13.
//  Copyright (c) 2015年 Yuyin Sun. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace std;
using namespace cv;

static void read_imageList(string imgList, vector<Mat>& images) {
    ifstream file(imgList.c_str(), ifstream::in);
    if (!file) {
        string error_message = "No valid input file was given";
        CV_Error(CV_StsBadArg, error_message);
    }
    string line;
    while (getline(file, line)) {
        images.push_back(imread(line, CV_LOAD_IMAGE_UNCHANGED));
    }
}

static Mat formatImageForPCA(const vector<Mat>& data) {
    Mat dst(static_cast<int>(data.size()), data[0].rows*data[0].cols, CV_32F);
    
    for (unsigned int i = 0; i < data.size(); ++i) {
        Mat image_row = data[i].clone().reshape(1, 1);
        Mat row_i = dst.row(i);
        image_row.convertTo(row_i, CV_32F);
    }
    
    return dst;
}

static Mat toGrayscale(InputArray _src) {
    Mat src = _src.getMat();
    if (src.channels() != 1) {
        string error_message = "Only Matrices with one channel are supported";
        CV_Error(CV_StsBadArg, error_message);
    }
    
    Mat dst;
    cv::normalize(_src, dst, 0, 255, NORM_MINMAX, CV_8UC1);
    
    return dst;
}

int main(int argc, const char * argv[]) {
    
    string imgList = "/Users/sunyuyin/Documents/Workspace/CPP/opencv_exercise/data/img_list.txt";
    
    vector<Mat> images;
    
    try {
        read_imageList(imgList, images);
    } catch (cv::Exception& e) {
        cerr << "Error open file \"" << imgList << "\". Reason:" << e.msg << endl;
        exit(1);
    }
    
    if (imgList.size() <= 1) {
        string error_message = "This demo needs a least 2 images";
        CV_Error(CV_StsError, error_message);
    }
    
    Mat data = formatImageForPCA(images);
    PCA pca(data, cv::Mat(), 0, 0.95);
    
    unsigned int image_index = 10;
    
    Mat point = pca.project(data.row(image_index));
    Mat reconstruction = pca.backProject(point);
    
    reconstruction = reconstruction.reshape(images[image_index].channels(), images[image_index].rows);
    reconstruction = toGrayscale(reconstruction);
    
    string winName_ori = "Original Image";
    namedWindow(winName_ori, WINDOW_NORMAL);
    
    imshow(winName_ori, images[image_index]);
    
    string winName = "Reconstruction";
    namedWindow(winName, WINDOW_NORMAL);
    
    imshow(winName, reconstruction);
    
    int key = 0;
    while (key != 'q') {
        key = waitKey();
    }
    
    destroyAllWindows();
    
    return 0;
}

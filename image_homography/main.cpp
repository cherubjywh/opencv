//
//  main.cpp
//  image_homography
//
//  Created by Yuyin Sun on 15-11-23.
//  Copyright (c) 2015年 Yuyin Sun. All rights reserved.
//

#include <iostream>
#include <limits>
#include <opencv2/core/core.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/nonfree/nonfree.hpp>
#include <opencv2/imgproc/imgproc.hpp>

using namespace cv;
using namespace std;


void call_existing_code() {
    Mat image1= imread( "/Users/sunyuyin/Desktop/img_stiching/img1.JPG" );
    Mat image2= imread( "/Users/sunyuyin/Desktop/img_stiching/img2.JPG" );
    Mat gray_image1;
    Mat gray_image2;
    // Convert to Grayscale
    cvtColor( image1, gray_image1, CV_RGB2GRAY );
    cvtColor( image2, gray_image2, CV_RGB2GRAY );
    
    imshow("first image",image2);
    imshow("second image",image1);
    
    if( !gray_image1.data || !gray_image2.data )
    { std::cout<< " --(!) Error reading images " << std::endl; return; }
    
    //-- Step 1: Detect the keypoints using SURF Detector
    int minHessian = 400;
    
    SurfFeatureDetector detector( minHessian );
    
    std::vector< KeyPoint > keypoints_object, keypoints_scene;
    
    detector.detect( gray_image1, keypoints_object );
    detector.detect( gray_image2, keypoints_scene );
    
    //-- Step 2: Calculate descriptors (feature vectors)
    SurfDescriptorExtractor extractor;
    
    Mat descriptors_object, descriptors_scene;
    
    extractor.compute( gray_image1, keypoints_object, descriptors_object );
    extractor.compute( gray_image2, keypoints_scene, descriptors_scene );
    
    //-- Step 3: Matching descriptor vectors using FLANN matcher
    FlannBasedMatcher matcher;
    std::vector< DMatch > matches;
    matcher.match( descriptors_object, descriptors_scene, matches );
    
    double max_dist = 0; double min_dist = 100;
    
    //-- Quick calculation of max and min distances between keypoints
    for( int i = 0; i < descriptors_object.rows; i++ )
    { double dist = matches[i].distance;
        if( dist < min_dist ) min_dist = dist;
        if( dist > max_dist ) max_dist = dist;
    }
    
    printf("-- Max dist : %f \n", max_dist );
    printf("-- Min dist : %f \n", min_dist );
    
    //-- Use only "good" matches (i.e. whose distance is less than 3*min_dist )
    std::vector< DMatch > good_matches;
    
    for( int i = 0; i < descriptors_object.rows; i++ )
    { if( matches[i].distance < 3*min_dist )
    { good_matches.push_back( matches[i]); }
    }
    std::vector< Point2f > obj;
    std::vector< Point2f > scene;
    
    for( int i = 0; i < good_matches.size(); i++ )
    {
        //-- Get the keypoints from the good matches
        obj.push_back( keypoints_object[ good_matches[i].queryIdx ].pt );
        scene.push_back( keypoints_scene[ good_matches[i].trainIdx ].pt );
    }
    
    // Find the Homography Matrix
    Mat H = findHomography( obj, scene, CV_RANSAC );
    // Use the Homography Matrix to warp the images
    cv::Mat result;
    warpPerspective(image1,result,H,cv::Size(image1.cols+image2.cols,image1.rows));
    
    
    cv::Mat half(result,cv::Rect(0,0,image2.cols,image2.rows));
    image2.copyTo(half);
    imshow( "Result", result );
    
    waitKey(0);
}

void call_my_code() {
    
//    Mat img1_rgb = imread("/Users/sunyuyin/Desktop/img_stiching/img2_online.jpg");
//    Mat img2_rgb = imread("/Users/sunyuyin/Desktop/img_stiching/img1_online.jpg");
    
    
    Mat img1_rgb = imread("/Users/sunyuyin/Desktop/img_stiching/img1.JPG");
    Mat img2_rgb = imread("/Users/sunyuyin/Desktop/img_stiching/img2.JPG");

    
    if (img1_rgb.empty() || img2_rgb.empty()) {
        exit(-1);
    }
    
    Mat img1, img2;
    
    cvtColor(img1_rgb, img1, CV_RGB2GRAY);
    cvtColor(img2_rgb, img2, CV_RGB2GRAY);
    
    
    SiftFeatureDetector detector;
    vector<KeyPoint> keypoints_img1, keypoints_img2;
    detector.detect(img1, keypoints_img1);
    detector.detect(img2, keypoints_img2);
    
    SiftDescriptorExtractor extractor;
    
    Mat descriptors_img1, descriptors_img2;
    
    extractor.compute(img1, keypoints_img1, descriptors_img1);
    extractor.compute(img2, keypoints_img2, descriptors_img2);
    
    FlannBasedMatcher matcher;
    vector<DMatch> matches;
    matcher.match(descriptors_img1, descriptors_img2, matches);
    
    double max_dist = 0.0, min_dist = numeric_limits<double>::max();
    
    for (int i = 0; i < descriptors_img1.rows; ++i) {
        double dist = matches[i].distance;
        if (dist < min_dist) min_dist = dist;
        if (dist > max_dist) max_dist = dist;
    }
    
    cout << "-- Max dist: " << max_dist << endl;
    cout << "-- Min dist: " << min_dist << endl;
    
    
    vector<DMatch> good_matches;
    for (int i = 0; i < descriptors_img1.rows; ++i) {
        if (matches[i].distance < 3 * min_dist) {
            good_matches.push_back(matches[i]);
        }
    }
    
    vector<Point2f> img1_matches;
    vector<Point2f> img2_matches;
    
    for (int i = 0; i < good_matches.size(); ++i) {
        img1_matches.push_back(keypoints_img1[good_matches[i].queryIdx].pt);
        img2_matches.push_back(keypoints_img2[good_matches[i].trainIdx].pt);
    }
    
    Mat H = findHomography(img1_matches, img2_matches, CV_RANSAC);
    
    Mat result;
    
    cout << H << endl;
    
    
    warpPerspective(img1_rgb, result, H, Size(img1_rgb.cols+img2_rgb.cols,img1_rgb.rows));
    
    /*
     namedWindow("img1");
     imshow("img1", img1);
     namedWindow("result");
     imshow("result", result);
     */
    
    
    
     Mat half(result, Rect(0, 0, img2_rgb.cols, img2_rgb.rows));
     img2_rgb.copyTo(half);
     imshow("Result", result);
    
    /*
    
    Mat img_matches;
    drawMatches(img1, keypoints_img1, img2, keypoints_img2, good_matches, img_matches, Scalar::all(-1), Scalar::all(-1),
                vector<char>(), DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS);
    
    namedWindow("img_matches");
    imshow("img_matches", img_matches);
    
    */
    /*
     
     Mat img1_output, img2_output;
     
     drawKeypoints(img1, keypoints_img1, img1_output);
     namedWindow("Image 1 keypoints", WINDOW_AUTOSIZE);
     
     imshow("Image 1 keypoints", img1_output);
     
     
     drawKeypoints(img2, keypoints_img2, img2_output);
     namedWindow("Image 2 keypoints");
     
     imshow("Image 2 keypoints", img2_output);
     */

}

int main(int argc, const char * argv[]) {

    call_existing_code();
//    call_my_code();
    
    waitKey(0);
    destroyAllWindows();
    
    return 0;
}

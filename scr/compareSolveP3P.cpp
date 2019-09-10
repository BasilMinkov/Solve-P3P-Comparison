// Standard libraries.
#include <iostream>
#include <cmath>
#include <vector>

// Eigen main library.
#include <Eigen/Dense>

// OpenCV sublibraries.
#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/calib3d/calib3d.hpp"
#include "opencv2/highgui/highgui.hpp"

#include "P3P_MinkovSawada.h"
#include "TdS_Math.h"
#include "Groebner_P3P.h"

using namespace std;
using namespace Eigen;

#define pi 3.14159265358979323846264338327950288


/**
    Computes a random value in a range.

    @param min Lower border of a range.
    @param max Upper border of a range.
    @return A random integer in a given range.
*/
double RandNew(double vmin, double vmax){return (double)rand()*(vmax - vmin) / RAND_MAX + vmin;}


/**
    Computes the angle between two vectors.

    @param firstVector
    @param secondVector
    @return Angle between vectors in degrees.
*/
double angleBetweenVectors(Vector3d firstVector, Vector3d secondVector) 
{
	return acos((firstVector.transpose() * secondVector / (firstVector.norm() * secondVector.norm()))(0));
}


/**
    Computes the visual angles of a simulated triangle in a 3D scene. 
    The visual angle of an object is the angle formed by rays 
    projecting from the eye (0, 0, 0) to the sides of an object.

    @param vertexMatrix Matrix of triangle vertexes (coordinateXYZ X vertexPointN) in 3D scene.
    @return Vector of visual angles in degrees.
*/
Vector3d ComputeVisualAngles(Matrix3d vertexMatrix)
{
	Vector3d visualAngles;

	for (int vertex = 0; vertex < 3; vertex++)
		{
			int vertex1 = vertex; 
			int vertex2 = vertex + 1; 
			if (vertex2 > 2) {vertex2 = 0;}
			visualAngles(vertex) = angleBetweenVectors(vertexMatrix.col(vertex1), vertexMatrix.col(vertex2));
		}

	return visualAngles;
}


/**
    Computes the angles of a simulated triangle in a 3D scene.

    @param vertexMatrix Matrix of triangle vertexes (coordinateXYZ X vertexPointN) in 3D scene.
    @return Vector of angles of a simulated triangle.
*/
Vector2d ComputeSimulatedTriangleAngles(Matrix3d vertexMatrix)
{
	Vector2d simulatedTriangleAngles;

	for (int vertex = 0; vertex < 2; vertex++) 
	{
		int vertex1 = vertex; 
		int vertex2 = vertex + 1; 
		int vertex3 = vertex + 2; 
		if (vertex3 > 2) {vertex3 = 0;}
		simulatedTriangleAngles(vertex) = angleBetweenVectors(vertexMatrix.col(vertex1)-vertexMatrix.col(vertex2), vertexMatrix.col(vertex2)-vertexMatrix.col(vertex3));
	}

	return simulatedTriangleAngles;
}


/**
    Foo

    @param 
    @return 
*/
vector<cv::Point2f> GenerateImagePoints(Matrix3d vertexMatrix)
{
	vector<cv::Point2f> points;

	// Push the same points as in the vertex matrix.

	// for (int vertex = 0; vertex < 3; vertex++)
	// {
	// 	points.push_back(cv::Point2f(
	// 		vertexMatrix.col(vertex)(0), 
	// 		vertexMatrix.col(vertex)(1))
	// 	);
	// }

	// Push the points in discussed way.

	points.push_back(cv::Point2f(0, 0));
	points.push_back(cv::Point2f(0, 1));
	points.push_back(cv::Point2f(vertexMatrix.col(2)(0), vertexMatrix.col(2)(1)));

  return points;
}


/**
    Foo

    @param 
    @return 
*/
vector<cv::Point3f> GenerateObjectPoints(Matrix3d vertexMatrix)
{
	vector<cv::Point3f> points;

	// Push the same points as in the vertex matrix.
 
	// for (int vertex = 0; vertex < 3; vertex++)
	// {
	// 	points.push_back(cv::Point3f(
	// 		vertexMatrix.col(vertex)(0),
	// 		vertexMatrix.col(vertex)(1),
	// 		vertexMatrix.col(vertex)(2))
	// 	);
	// }

	// Push the points in discussed way.

	points.push_back(cv::Point3f(0, 0, 0));
	points.push_back(cv::Point3f(0, 1, 0));
	points.push_back(cv::Point3f(vertexMatrix.col(2)(0), vertexMatrix.col(2)(1), 0));
  
	return points;
}


void printImagePoints(vector<cv::Point2f> vector)
{
	cout << "Image Coordinates In 2D Space (coordinateXY X vertexPointN): << ";

	for (int i=0; i<vector.size(); i++)
	{
		cout << vector[i].x << ", " << vector[i].y << "; ";
	}
	cout << endl;
}


void printObjectPoints(vector<cv::Point3f> vector)
{
	cout << "Object Coordinates In 3D Space (coordinateXYZ X vertexPointN): << ";

	for (int i=0; i<vector.size(); i++)
	{
		cout << vector[i].x << ", " << vector[i].y << ", " << vector[i].z << "; ";
	}
	cout << endl;
}


void P3P_OpenCV(vector<cv::Point2f> imagePoints, vector<cv::Point3f> objectPoints)
{

	cv::Mat cameraMatrix(3, 3, cv::DataType<double>::type); // set camera matrix
	cv::setIdentity(cameraMatrix); 

	cv::Mat distCoeffs(4, 1, cv::DataType<double>::type); // set distortion coeffs
	distCoeffs.at<double>(0) = 0;
	distCoeffs.at<double>(1) = 0;
	distCoeffs.at<double>(2) = 0;
	distCoeffs.at<double>(3) = 0;

	cv::Mat rvec(3,1,cv::DataType<double>::type);
	cv::Mat tvec(3,1,cv::DataType<double>::type);

	cv::solvePnP(objectPoints, imagePoints, cameraMatrix, distCoeffs, rvec, tvec, "CV_P3P");

	std::vector<cv::Point2f> projectedPoints;
	cv::projectPoints(objectPoints, rvec, tvec, cameraMatrix, distCoeffs, projectedPoints);
 
	std::cout << "	rvec: " << rvec << std::endl;
	std::cout << "	tvec: " << tvec << std::endl;

	for(unsigned int i = 0; i < projectedPoints.size(); ++i)
	{
		std::cout << "	Image point: " << imagePoints[i] << " Projected to " << projectedPoints[i] << std::endl;
	}
}


int main()
{

	// Set random seed.
	int rseed = (unsigned int)time(0);
	srand(rseed);

	// Eigen output format constants. 
	IOFormat CommaInitFmt(StreamPrecision, DontAlignCols, ", ", "; ", "", "", " << ", ";"); // set the output format (vector [x1, x2, x3, x4] as "<< x1, x2, x3, x4;")
	string SEP = "----------------------------------------\n"; // separator to be able to tell one output from another

	// Simulation consts.
	int NUM_REPEAT = 10;
	int NUM_POINTS = 3;
	int MIN_Z = 1;
	int MAX_Z = 1000;

	// Variables.
	Matrix3d vertexMatrix; // matrix of triangle vertexes (coordinateXYZ X vertexPointN)
	Vector3d visualAngles; // vector of visual angles A&B, B&C, and C&A


	// Initialase the comparison simulation loop. The aim is to obtain the statistics 
	for (int trial = 0; trial < NUM_REPEAT; trial++)
	{

		cout << "Trial #" << trial << endl << endl;
		// int rseed = (unsigned int)time(0);
		// srand(rseed);

		// Generate vertex matrix in 3D space.
		for (int p = 0; p < NUM_POINTS; p++)
		{	
			// Generate Z coordinate. 
			vertexMatrix(2, p) = RandNew(MIN_Z, MAX_Z); // Z

			// Generate X and Y spherical coordinate based on Z. 
			double eccentricity = RandNew(0, 85) * pi / 180.0;
			double ori_xy = RandNew(0, 360) * pi / 180.0;
			vertexMatrix(0, p) = vertexMatrix(2, p) * tan(eccentricity) * cos(ori_xy); // X
			vertexMatrix(1, p) = vertexMatrix(2, p) * tan(eccentricity) * sin(ori_xy); // Y
		}

		cout << "Vertex Matrix (coordinateXYZ X vertexPointN):" << vertexMatrix.format(CommaInitFmt) << endl; 

		// Compute input for three simulations.

		cout << "\nInput:\n\n";

		//Input-A

		//A1) Compute visual angles between A&B, B&C, and C&A
		Vector3d visualAngles = ComputeVisualAngles(vertexMatrix);
		cout << "Visual Angles:" << visualAngles.format(CommaInitFmt) << endl; 

		//A2) Compute 2D image coordinates (camera matrix projecting a 3D scene to a 2D image plane)
  		vector<cv::Point2f> imagePoints = GenerateImagePoints(vertexMatrix); // matrix of triangle vertexes (coordinateXYZ X vertexPointN)
  		printImagePoints(imagePoints);


		//Input-B

		//B1) Shape of the simulated triangle in a 3D scene: 2 of 3 vertices of the triangles: [angle of vertex_A, angle of vertex_B]
		Vector2d simulatedTriangleAngles = ComputeSimulatedTriangleAngles(vertexMatrix);
		cout << "Simulated Triangle Angles:" << simulatedTriangleAngles.format(CommaInitFmt) << endl; 

		//B2) Shape of the simulated trianlge in a 3D scene: 3D coordinates: (0,0,0), (0,1,0), (x3,y3,0)
  		std::vector<cv::Point3f> objectPoints = GenerateObjectPoints(vertexMatrix); // matrix of triangle vertexes (coordinateXYZ X vertexPointN)
  		printObjectPoints(objectPoints);

		// Compute simulations. 

		cout << "\nOutput:\n\n";

		//P3P_MinkovSawada: A1 and B1

		double interpretations[10][3];
		cout << "P3P_MinkovSawada: " << P3P_MinkovSawadaB(
			simulatedTriangleAngles(0), 
			simulatedTriangleAngles(1), 
			visualAngles(0), 
			visualAngles(1), 
			visualAngles(2), 
			interpretations
			) << " sol." << endl; // Tada, why don't you define the last argument inside of the function? Is it for debuging 

		//P3P_Banno: ? and ?

		// CGroebner::CalculatePosition()


		//P3P_OpenCV: A2 and B2
		cout << "P3P_OpenCV:" << endl;
		P3P_OpenCV(imagePoints,objectPoints);

		cout << SEP;

	} // trial

	cout << "NOTE: For the case of matrix output, comma separates rows, while semicolon separates columns." << endl;

	return 0;
}
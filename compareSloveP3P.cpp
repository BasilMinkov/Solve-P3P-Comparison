#include <iostream>
#include <vector>
#include <cmath>

// OpenCV sublibraries
#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/calib3d/calib3d.hpp"
#include "opencv2/highgui/highgui.hpp"

using namespace std;

#define pi 3.14159265358979323846264338327950288

int RandNew(int min, int max){return rand() % max + min; }

int main()
{
	int num_repeat = 10000;
	int numPoints = 3;
	int minZ = 1;
	int maxZ = 1000;
	double px[3], py[3], pz[3];


	for (int trial = 0; trial < num_repeat; trial++)
	{
		int rseed = (unsigned int)time(0);
		srand(rseed);

		for (int p = 0; p < numPoints; p++)
		{
			pz[p] = RandNew(minZ, maxZ);

			double eccentricity = RandNew(0, 85) * pi / 180.0;
			double ori_xy = RandNew(0, 360) * pi / 180.0;
			px[p] = pz[p] * tan(eccentricity) * cos(ori_xy);
			py[p] = pz[p] * tan(eccentricity) * sin(ori_xy);
		}

		//Input-A
		//A1) Compute visual angles between A&B, B&C, and C&A
		//A2) Compute 2D image coordinates (camera matrix projecting a 3D scene to a 2D image plane)

		//Input-B
		//B1) Shape of the simulated triangle in a 3D scene: 2 of 3 vertices of the triangles: [angle of vertex_A, angle of vertex_B]
		//B2) Shape of the simulated trianlge in a 3D scene: 3D coordinates: (0,0,0), (0,1,0), (x3,y3,0)

		//P3P_MinkovSawada: A1 and B1
		//P3P_Banno: ? and ?
		//P3P_OpenCV: A2 and B2

	} // trial

	return 0;
}
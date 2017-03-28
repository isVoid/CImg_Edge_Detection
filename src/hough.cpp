/*
*	Hough Algorithm
*	Michael Wang SDCS SYSU 2017
*/


#include <iostream>
#include <vector>
#include <list>
#include <cmath>

#define cimg_use_jpeg
#include "CImg.h"
#include "canny.h"

#define HOUGH_SPACE_SIZE 600
#define hss HOUGH_SPACE_SIZE
#define NEARST_NEIGHBOUR 1

using namespace cimg_library;
using namespace std;

CImg<unsigned char> img;
CImg<unsigned char> resized;
CImg<unsigned char> cny;
CImg<int> hough_space;
CImg<unsigned char> result;

//Canny Parameter
int gfs = 3;
double g_sig = 1.0;
int thres_lo = 20;
int thres_hi = 40;

//Hough Parameter
int voting_thres = 64;
float resize_fac = 8;
float thres_fac = 0.4;
double PI = 3.1415;

struct param_space_point {
	int rho = 0;
	double theta = 0;
	int vote = 0;

	param_space_point(int r, double t, int v) {
		rho = r; theta = t; vote = v;
	}
};

struct point
{
	int x = 0;
	int y = 0;

	point() {}

	point(int _x, int _y) {
		x = _x;
		y = _y;
	}	
};

void plotPoint(int x, int y, CImg<unsigned char>& img) {

	int r = 30;

	printf("Plotting... x:%d, y:%d\n", x, y);
	printf("w: %d, h: %d\n", img._width, img._height);

	for (int i = x - r; i < x + r; i++) {
		for (int j = y - r; j < y + r; j++) {
			// if (x == 2248) {
			// 	printf("%d, %d\n", i, j);
			// }
			if (i >= 0 && i < img._width && j >= 0 && j < img._height) {
				if (sqrt(pow(i - x, 2) + pow(j - y, 2)) <= r) {
					if (x == 2248) {
						printf("Draw... i:%d, j:%d\n", i, j);
					}
					// printf("Draw... i:%d, j:%d\n", i, j);
					img(i, j, 0) = 255;
					img(i, j, 1) = 255;
					img(i, j, 2) = 0;
				}
			}

		}
	}

	// img.display();

}

void plotLine(vector<param_space_point> filtered, CImg<unsigned char>& result) {

	for (int i = 0; i < filtered.size(); i++) {

		param_space_point p = filtered[i];

			cimg_forXY(result, x, y) {
				if (abs(p.rho - (x/resize_fac) * cosf(p.theta) - (y/resize_fac) * sinf(p.theta)) < 1) {
					result(x, y, 0) = 255;
					result(x, y, 1) = 0;
					result(x, y, 2) = 0;
					// result(x, y) = 255;
						
				}
			}

			cout << "Finish Plotting Line" << endl;

	}
}

void plotLineFast(vector<param_space_point> filtered, CImg<unsigned char>& result) {

	for (int i = 0; i < filtered.size(); i++) {

		param_space_point p = filtered[i];
		int w = result._width;
		int h = result._height;

		if (abs(sin(p.theta)) > 0.5) {
			float k = cos(p.theta) / sin(p.theta);
			float b = p.rho / sin(p.theta);
			// printf("k: %f, b: %f\n", k, b);

			// cimg_forX(result, x) {
			for (float x = 0; x < w; x += 0.5){
				int y = k * x + b;
				// printf("x: %d, y: %d\n", x, y);

				if (x >= 0 && x <= w && y >= 0 && y <= h){
					// printf("x: %d, y: %d\n", x, y);
					// n_plot++;

					result(int(x), y, 0) = 255;
					result(int(x), y, 1) = 0;
					result(int(x), y, 2) = 0;	
				}
				
			}
		}
		else {
			float k = tan(p.theta);
			float b = -p.rho / cos(p.theta);
			// printf("k: %f, b: %f\n", k, b);

			// cimg_forY(result, y) {
			for (float y = 0; y < h; y += 0.5){
				int x = k * y + b;
				// printf("x: %d, y: %d\n", x, y);
				if (x >= 0 && x <= w && y >= 0 && y <= h){
					// printf("x: %d, y: %d\n", x, y);
					// n_plot++;

					result(x, int(y), 0) = 255;
					result(x, int(y), 1) = 0;
					result(x, int(y), 2) = 0;	
				}
			}
		}

	}

}

void localFiltering(vector<param_space_point> v, vector<param_space_point>& filtered) {

	for (int i = 0; i < v.size(); i++) {

		bool add = true;
		param_space_point p = v[i];

		for (int j = 0; j < filtered.size(); j++) {

			param_space_point p0 = filtered[j];
			float dist = sqrt(pow(p.theta - p0.theta, 2) + pow((p.rho - p0.rho), 2));

			if (dist < 200 ) {
				add = false;
			} else {
				printf("Param dist: %f\n", dist);
			}
		}

		if (add) {
			filtered.insert(filtered.end(), p);
		}

	}

	printf("Filtered Param points: %d\n", filtered.size());

}

void computeIntersects(vector<param_space_point> filtered, vector<point>& intersects) {

for (int i = 0; i < filtered.size(); i++) {

		for (int j = i + 1; j < filtered.size(); j++) {

			param_space_point p0 = filtered[i], p1 = filtered[j];
			float a = cos(p0.theta), b = sin(p0.theta), e = p0.rho,
				c = cos(p1.theta), d = sin(p1.theta), f = p1.rho;

			float det = (a*d-b*c);
			if (abs(det) > 0.01) {
			
				float x = (e*d-b*f)/det;
				float y = (a*f-e*c)/det;

				point intersect;
				intersect.x = int(x) * resize_fac;
				intersect.y = int(y) * resize_fac;
				plotPoint(intersect.x, intersect.y, result);

				intersects.insert(intersects.end(), intersect);
			}

		}

	}

	printf("Intersects: %d\n", intersects.size());
}

vector<param_space_point> v;
vector<param_space_point> v_debug;
vector<param_space_point> filtered;
vector<point> intersects;

void hough_line() {

	int w = cny._width, h = cny._height;
	int pmax = max(w, h);
	pmax = ceil(1.414 * pmax);

	printf("w:%d, h:%d, pmax:%d\n",w,h,pmax);

	hough_space.assign(360, 2 * pmax, 1, 1, 0);
 
	int max = 0;

	cimg_forXY(cny, x, y) {

		int v = cny(x, y);

		if (v > voting_thres) {
			// printf("x, y: %d %d\n", x, y);
			for (int th = 0; th < 360; th++) {

				float angle = 2*PI * th / 360;
				float xcos = x * cos(angle);
				float ysin = y * sin(angle);
				float rho = xcos + ysin;
				
				float rho_s = rho + pmax;
				// int rho_norm = float(rho_s) * hss / (2 * pmax);

				// if (rho_norm < 0 || rho_norm >= hss) {
				// 	cout << "x:" << x << ' ' << "y: " << y << endl;
				// 	cout << "rho_norm:" << rho_norm << ' ' << "theta:" << th << endl;	
				// 	printf("Angle: %f\n", angle);
				// 	printf("cos: %f, sin: %f\n", cos(angle), sin(angle));
				// 	printf("-xcos: %f, ysin: %f\n", -xcos, ysin);
				// 	printf("pmax: %d\n", pmax);
				// 	printf("raw rho: %d\n", rho);
				// 	printf("\n");
				// }

				// cout << "rho_norm:" << rho_norm << ' ' << "theta:" << th << endl;
				// printf("rho_norm: %d\n", rho_norm);
				hough_space(th, int(rho_s))++;

				if (hough_space(th, rho_s) > max) {
					max = hough_space(th, rho_s);
					// printf("max %u\n", max);
				}
			}

		}
	}

	hough_space.display();

	// exit(-1);
	cout << "max: " << max << endl;

	int thresh = max * thres_fac;
	cimg_forXY(hough_space, th, rho_norm) {

		if (hough_space(th, rho_norm) < thresh){
			hough_space(th, rho_norm) = 0;
		}
		else {
			int rho = rho_norm - pmax;
			double theta = 2 * PI * th / 360;
			v.insert(v.end(), param_space_point(rho, theta, hough_space(th, rho_norm)));
		}
	}

	printf("thres amount: %u\n", v.size());

	localFiltering(v, filtered);

	plotLine(filtered, result);

	computeIntersects(filtered, intersects);

}

int main(int argc, char** argv) {

	char* original_path = "./data/2.jpg";
	// char* canny_path = "./data/4c.jpg";

	img.assign(original_path);
	resized.assign(original_path);
	resized.resize(resized._width / resize_fac, resized._height / resize_fac);
	canny c(resized);

	int gfs = 5;
	double g_sig = 1.5;
	int thres_lo = 30;
	int thres_hi = 55;


	cny = c.process(gfs, g_sig, thres_lo, thres_hi);
	cny.display();
	// cny.assign(canny_path); 
	result.assign(img);

	hough_line();

	hough_space.display();
	result.display();
	// cny.display();

	return 0;
}
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

bool debug_disp = false;

double PI = 3.1415;

float resize_fac = 8;
//Canny Parameter
int gfs = 3;
double g_sig = 1.0;
int thres_lo = 20;
int thres_hi = 40;

//Hough Parameter
int voting_thres = 64;
float thres_fac = 0.4;
int filter_thres = 200;

/*
*	A Point in Hough Param Space
*	Properties:
*	Rho, Integer, from -rho_max to rho_max 
*	Theta, Float, from 0 to 2PI
*	Vote_count, Integer, Positive
*/
struct param_space_point {
	int rho = 0;
	double theta = 0;
	int vote = 0;

	param_space_point(int r, double t, int v) {
		rho = r; theta = t; vote = v;
	}
};

/*
*	A Point in Image Space
*	Properties:
*	x: Integer, from 0 to image.width
*	y: Integer, from 0 to image.height
*/

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

/*
*	@Override
*	getHist
*	Calculates a normalized gray image's histogram
*	
*	@param:
*	T hist[256]: An array to store the images histogram
*	CImg<float> img: Image to calculate
*/
template <typename T> void getHist(T hist[256], CImg<float> img) {
	cimg_forXY(img, x, y) {
		hist[int(img(x,y) * 255)] += 1;
	}
}

/*
*	histgramEq
*	Equalized a histogram
*	
*	@param:
*	T hist[256]: A(n) (image/channel)'s histogram
*	int size: total pixel number of the image/channel, usually w x h.
*/
template <typename T> void histgramEq(T hist[256], int size) {
	
	//Caculate cummulative density function
	for(int i = 1; i < 256; i++) {
		hist[i] = hist[i] + hist[i-1];
	}

	//Compute color transform function for equalization
	for(int i = 0; i < 256; i++) {
		hist[i] = float(hist[i]) / size;
		hist[i] = (unsigned int)(255 * hist[i]);
	}

}

/*
*	histgramEq_hsi
*	Histogram-equalize a color image using the intensity channel
*	on hsi space.
*	
*	@param:
*	CImg<unsigned char> in: Image to perform operation
*
*	@return:
*	CImg<unsigned char>: Histogram equalized image
*/
CImg<unsigned char> histgramEq_hsi(CImg<unsigned char> in) {

	CImg<float>in_hsi = in.get_RGBtoHSI();

	int size = in._width * in._height;

	float hist_i[256] = {};

	CImg<float> intensity = in_hsi.get_channel(2);
	getHist(hist_i, intensity);

	histgramEq(hist_i, size);

	//Applying color transform function to image
	cimg_forXY(intensity, x, y) {
		in_hsi.atXYZC(x,y,1,2) = hist_i[int(intensity(x,y) * 255)] / 255.0;
	}

	CImg<unsigned char> out = in_hsi.get_HSItoRGB();

	return out;
}

/*
*	plotPoint function: A function to plot points on the image
*	@param
*	x: x coordinate
*	y: y coordinate
*	img: image to plot on
*/

void plotPoint(int x, int y, CImg<unsigned char>& img) {

	int r = 30;

	printf("Plotting... x:%d, y:%d\n", x, y);
	printf("w: %d, h: %d\n", img._width, img._height);

	for (int i = x - r; i < x + r; i++) {
		for (int j = y - r; j < y + r; j++) {
			if (i >= 0 && i < img._width && j >= 0 && j < img._height) {
				if (sqrt(pow(i - x, 2) + pow(j - y, 2)) <= r) {
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

/*
*	plotLine function: A function to plot line on the image, using full image scan 
*	and hough param formula as line formula, O(wh)
*	@param
*	filtered: a vector of hough prameter space points, each represents a line in image space
*	result: the result image to plot on
*/

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

/*
*	plotLineFast function: A function to plot line faster (HAS NOT BEEN TESTED YET!!!)
*	Compute Slope-Intercept Form from Hough Param Formula, scan image horizontally once. O(w)
*	@param
*	filtered: a vector of hough prameter space points, each represents a line in image space
*	result: the result image to plot on
*/

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

/*
*	localFiltering: A naiive implementation of local point filtering, choose the first point
*	encountered in the local area. Going to be replaced by K-means classification
*	@param
*	v: list of hough space points
*	filtered: list of filtered points
*/

void localFiltering(vector<param_space_point> v, vector<param_space_point>& filtered) {

	for (int i = 0; i < v.size(); i++) {

		bool add = true;
		param_space_point p = v[i];

		for (int j = 0; j < filtered.size(); j++) {

			param_space_point p0 = filtered[j];
			float dist = sqrt(pow(p.theta - p0.theta, 2) + pow((p.rho - p0.rho), 2));

			if (dist < filter_thres ) {
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

/*
*	computeIntersects: Computing the intersects of different lines and the points in a vector
*	filtered: hough space points, each represents a line in image space
*	intersects: list of intersects
*/
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

/*
*	hough_line: hough line detection function
*	Step:
*	1. Traverse Canny image, vote on hough space
*	2. Traverse hough space, threshold usable points
*	3. Perform local filtering, reduce duplicates
*	4. Compute intersects of lines
*	5. Plot Lines and Intersect Points.
*/

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

	if (debug_disp)
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

/*
*	Main Function: Provides input of canny and hough parameters.
*	Parameter List, from left to right:
*
*	Path_Image_to_detect: string
*	Perform_Histogram_equalization: bool(0, 1)
*	Gaussian_filter_size: odd int
*	Gaussian_sigma: float
*	threshold_lower: int
*	threshold_higher: int
*	voting_threshold: int
*	threshold_factor: float
*	filtering_treshold: int
*	Show_debug_windows: bool(0, 1)
*	Path_Result_to_save: string
*
*	Sample usage: ./HOUGH_LINE "./data/5.jpg" 0 5 1.5 40 60 64 0.4 30 0 "./data/5.jpg"
*
*	For more information, refer to manual/report
*/

int main(int argc, char** argv) {

	char* original_path = argv[1];
	char* out_path = argv[11];

	bool do_hist_eq = bool(atoi(argv[2]));

	gfs = atoi(argv[3]);
	g_sig = atof(argv[4]);
	thres_lo = atoi(argv[5]);
	thres_hi = atoi(argv[6]);

	voting_thres = atoi(argv[7]);
	thres_fac = atof(argv[8]);
	filter_thres = atoi(argv[9]);

	debug_disp = bool(atoi(argv[10]));

	img.assign(original_path);
	resized.assign(original_path);
	resized.resize(resized._width / resize_fac, resized._height / resize_fac);

	if (do_hist_eq) {
		resized = histgramEq_hsi(resized);	
	}

	if (debug_disp)
		resized.display();
	canny c(resized);

	cny = c.process(gfs, g_sig, thres_lo, thres_hi);

	if (debug_disp)
		cny.display();

	result.assign(img);

	hough_line();

	if (debug_disp)
		hough_space.display();

	result.display();

	result.save_jpeg(out_path);
	printf("Result saved to ./output/\n");

	return 0;
}
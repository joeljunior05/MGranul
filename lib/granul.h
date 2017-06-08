#include "meikon.h"

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#if defined(ENABLE_OPENMP)
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_get_thread_num() { return 0;}
inline omp_int_t omp_get_max_threads() { return 1;}
inline omp_int_t omp_get_num_threads() { return 1;}
#endif

using namespace cv;
using namespace std;

// CONST VALUES

#define GRY_MAX 255
#define BOUND_PIXELS 2

// TYPES

typedef unsigned char BYTE;
typedef unsigned char GRY;
typedef double DBL;
typedef float FLT;
typedef Vec3b COR;
typedef complex<FLT> CPX;

// STRUCTS

struct MAXLOCAL {
  FLT r; // correlation (between -1 and +1)
  char forma; // c, q, r, ou e. Igual da FORMA
  int li,ci;
  double d1, d2, ang; // lado1, lado2, angulo em graus. lado1 e lado2 ja vem multiplicados pela escala.
  double area, area_int;
  bool    valid = true;

  MAXLOCAL() {}
  MAXLOCAL(FLT pr, char pforma, int pli, int pci, double pd1, double pd2, double pang)
    { r=pr; forma=pforma; li=pli; ci=pci; d1=pd1; d2=pd2; ang=pang; area=0; }
  MAXLOCAL(FLT pr, char pforma, int pli, int pci, double pd1, double pd2, double pang, double parea)
    { r=pr; forma=pforma; li=pli; ci=pci; d1=pd1; d2=pd2; ang=pang; area=parea; }
};

struct KERNEL
{
	Mat_<GRY> imgG;
	Mat_<FLT> imgF;
	char shape; // c(Circle), s(Square), r(Rectangle), e(Ellipse) or g(Generic)
	double d1, d2, ang; // dimension1(rows), dimension1(cols), angle in degrees
	KERNEL() {}
};

struct P_CORRELATION {
  int id;
  int li;
  int ci;

  P_CORRELATION() {}
  P_CORRELATION(int i, int l, int c) {
    id = i; li = l; ci = c;
  }
};

struct MM_CORRELATION {
  int id_g;
  int area_intersec;
  int area_complemento;
  int id_m;

  MM_CORRELATION() {}
  MM_CORRELATION(int g, int a_i, int a_c, int m) {
    id_g = g; area_intersec = a_i; area_complemento = a_c; id_m = m;
  }
};

typedef multimap<int,MAXLOCAL> MMAPA;
typedef MMAPA::iterator IMMAPA;

// AUX FUNCTIONS

int evenToOdd(double val); // An even number became an odd one.
void normalizeFLT(Mat_<FLT> &matF, FLT target, FLT overall);
void normalizeGRY(Mat_<GRY> &matG, GRY target, GRY overall);
void lcFromKernel(int l, int c, int& li, int& ci, KERNEL k);

void convertGRYToFLT(Mat in, Mat& out);
void convertFLTToGRY(Mat in, Mat& out);

void sortMAXLOCAL(vector<MAXLOCAL>& a, int l, int r);

bool isInsideKernel(int r, int c, MAXLOCAL m);
bool isMAXEqual(MAXLOCAL m, MAXLOCAL n);

void removeCloser(vector<MAXLOCAL>& v, vector<MAXLOCAL>& u, int maxDist);

void printMAXLOCAL(vector<MAXLOCAL>& v, Mat& mat, bool withTotal = false, bool withCorr = false);
void printTextMAXLOCAL(vector<MAXLOCAL>& v, char* fileName);
void printMAXLOCAL(MAXLOCAL m);

// KERNEL FUNCTIONS

void createKernelCircle(double diameter, vector<KERNEL>& output, int op=-1);
void createKernelRectangle(double lengthL, double lengthC, double angle, vector<KERNEL>& output, int op=-1);
void createKernelSquare(double length, double angle, vector<KERNEL>& output, int op=-1);
void createKernelEllipse(double lengthL, double lengthC, double angle, vector<KERNEL>& output, int op=-1);
void createKernelGeneric(Mat_<FLT> fltM, Mat_<GRY> gryM, vector<KERNEL>& output, double angle);
void createKernels(vector<KERNEL> &output, char type, double lengthL, double lengthC, double escOitava = 5, double escMax = 1.0, double escMin=0.5, double angStep = 10);

// GRANULOMETRY FUNCTIONS

void correlation(Mat in, Mat &out, KERNEL ker, int type=CV_TM_CCORR);
void correlationInBatch(Mat& img, vector<MAXLOCAL> &out, vector<KERNEL> kers, float minCorr=0.001, int maxDist=2, int type=CV_TM_CCORR);
void computeIntersection(Mat_<MAXLOCAL>& n_out, MAXLOCAL& loc, float maxInt=1);
void sift(vector<MAXLOCAL> v, vector<MAXLOCAL>& w, float minCorr, float maxInter, Size size);
void sift(vector<MAXLOCAL> in, vector<MAXLOCAL>& out, double minCorrCir, double maxInterCir, double minCorrRet, double maxInterRet,
              double minCorrQua, double maxInterQua, double minCorrEli, double maxInterEli, Size size);
void siftMSER(Mat img, vector<MAXLOCAL> in, vector<MAXLOCAL>& out, double minCorrMSER);
void chooseKmeans(Mat img, vector<MAXLOCAL> in, vector<MAXLOCAL>& out, double num_k);
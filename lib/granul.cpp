#include "granul.h"
#include <fstream>

// AUX FUNCTIONS

Mat_<GRY> global;

int evenToOdd(double val){
	int ret = ceil(val);
	if(ret%2 == 0)
	{
		return ret + 1;
	}
	return ret;

	//return 2*cvRound(ceil((val-1.0)/2.0))+1;
}

void sumAndShow(Mat_<FLT> &matF){
	float pos_counter = 0;
	float neg_counter = 0;

	for(int l = 0; l < matF.rows; l++){
		for(int c = 0; c < matF.cols; c++){
			if(matF.at<FLT>(l, c) > 0)
				pos_counter += matF.at<FLT>(l, c);
			else
				neg_counter += matF.at<FLT>(l, c);
		}
	}

	cout << "POS = " << pos_counter << endl;

	cout << "NEG = " << neg_counter << endl;
}

void NegPosUm(Mat_<FLT>& f)
// Faz a somatoria dos negativos dar -1 e dos positivos dar +1.
{ double pos=0.0;
  double neg=0.0;
  for (MatIterator_<FLT> fi=f.begin(); fi!=f.end(); fi++) {
    if (*fi>0.0) pos += *fi;
    else neg -= *fi;
  } 
  if (pos<=0.0 || neg<=0.0) erro("Erro: NegPosUm divisao zero");
  for (MatIterator_<FLT> fi=f.begin(); fi!=f.end(); fi++) {
    if (*fi>0.0) *fi /= pos;
    else *fi /= neg;
  }
}

void normalizeFLT(Mat_<FLT> &matF, FLT target, FLT overall){
	
	unsigned int counter = 0;

	for(int l = 0; l < matF.rows; l++){
		for(int c = 0; c < matF.cols; c++){
			if(matF.at<FLT>(l, c) == target)
				counter++;
		}
	}

	if(counter == 0)
		return;

	FLT piece = overall / ((float) counter);

	for(int l = 0; l < matF.rows; l++){
		for(int c = 0; c < matF.cols; c++){
			if(matF.at<FLT>(l, c) == target)
				matF.at<FLT>(l, c) = piece;
		}
	}
}


void normalizeGRY(Mat_<GRY> &matG, GRY target, GRY overall){
	unsigned int counter = 0;

	for(int l = 0; l < matG.rows; l++){
		for(int c = 0; c < matG.cols; c++){
			if(matG.at<GRY>(l, c) == target)
				counter++;
		}
	}

	if(counter == 0)
		return;

	FLT piece = overall / counter;

	for(int l = 0; l < matG.rows; l++){
		for(int c = 0; c < matG.cols; c++){
			if(matG.at<GRY>(l, c) == target)
				matG.at<GRY>(l, c) = piece;
		}
	}

}

void lcFromKernel(int l, int c, int& li, int& ci, KERNEL k){

	li = evenToOdd(k.d1 / 2) + l;
	ci = evenToOdd(k.d2 / 2) + c;
}

void convertGRYToFLT(Mat in, Mat& out){
	out.create(in.size(), CV_32F);

	for(int l = 0; l < in.rows; l++){
		for(int c = 0; c < in.cols; c++){
			out.at<FLT>(l, c) = (((FLT) in.at<GRY>(l, c)) / GRY_MAX);
		}
	}
}

void convertFLTToGRY(Mat in, Mat& out){
	out.create(in.size(), CV_8U);

	for(int l = 0; l < in.rows; l++){
		for(int c = 0; c < in.cols; c++){
			out.at<GRY>(l, c) = (in.at<FLT>(l, c) * GRY_MAX);
		}
	}
}

void sortMAXLOCAL(vector<MAXLOCAL>& a, int l, int r)
{
	assert(l<=r);
	int i=l; 
	int j=r;
	MAXLOCAL x=a[(l+r)/2];

	do {
		while (a[i].r>x.r) i++;
		while (x.r>a[j].r) j--;

		if (i<=j) {
			swap(a[i],a[j]); i++; j--;
		}
	} while (i<=j);

	if (l<j) sortMAXLOCAL(a,l,j);
	if (i<r) sortMAXLOCAL(a,i,r);
}


bool isInsideKernel(int r, int c, MAXLOCAL m){
	// c(Circle), s(Square), r(Rectangle), e(Ellipse) or g(Generic)
	char type = m.forma;

	if (type =='c'){ //case circle

		FLT distance = sqrt((m.li - r)*(m.li - r) + (m.ci - c)*(m.ci - c));

		if(distance <= ceil(m.d1/2)){
			return true;
		}

	} else if (type =='s' || type =='r'){ //case rectangle

		Mat rot_mat = getRotationMatrix2D(Point2f(m.ci, m.li), m.ang, 1.0);

		FLT rot_r, rot_c;

		rot_c = rot_mat.at<double>(0,0) * c + rot_mat.at<double>(0,1) * r + rot_mat.at<double>(0,2);
		rot_r = rot_mat.at<double>(1,0) * c + rot_mat.at<double>(1,1) * r + rot_mat.at<double>(1,2);

		if(rot_c >= ceil(m.ci - m.d2/2) && rot_c <= ceil(m.ci + m.d2/2) 
			&& rot_r >= ceil(m.li - m.d1/2) && rot_r <= ceil(m.li + m.d1/2)) return true;


	} else if (type =='e'){ //case ellipse

		Mat rot_mat = getRotationMatrix2D(Point2f(m.ci, m.li), m.ang, 1.0);

		FLT rot_r, rot_c;

		rot_c = rot_mat.at<double>(0,0) * c + rot_mat.at<double>(0,1) * r + rot_mat.at<double>(0,2);
		rot_r = rot_mat.at<double>(1,0) * c + rot_mat.at<double>(1,1) * r + rot_mat.at<double>(1,2);

		rot_c = abs(rot_c - m.ci);
		rot_r = abs(rot_r - m.li);

		return (rot_c * rot_c) / ceil(m.d2/2 * m.d2/2) + (rot_r * rot_r) / ceil(m.d1/2 * m.d1/2)
                <= 1.0;

	} else {
		cerr << type << " : Type Invalid" << endl;

	}

	return false;
}

bool isMAXEqual(MAXLOCAL m, MAXLOCAL n){
	return m.r == n.r 
	&& m.forma == n.forma
	&& m.li == n.li
	&& m.ci == n.ci
	&& m.d1 == n.d1
	&& m.d2 == n.d2;
}

void removeCloser(vector<MAXLOCAL>& v, vector<MAXLOCAL>& u, int maxDist){ 

	MMAPA m;

	sortMAXLOCAL(v,0,v.size()-1);

	for (unsigned i=0; i<v.size(); i++) {
		IMMAPA low=m.lower_bound(v[i].ci-maxDist);
		IMMAPA up=m.upper_bound(v[i].ci+maxDist);
		for (IMMAPA p=low; p!=up; p++) {
			if (v[i].forma==p->second.forma && abs(v[i].li - p->second.li)<=maxDist) {
				assert(abs(v[i].ci-p->second.ci)<=maxDist);
				goto saida;
			}
		}
		m.insert(make_pair(v[i].ci,v[i]));
		saida:;
	}
	u.clear();
	for (IMMAPA p=m.begin(); p!=m.end(); p++) 
		u.push_back(p->second);

	sortMAXLOCAL(u,0,u.size()-1);
}

void printMAXLOCAL(vector<MAXLOCAL>& v, Mat& mat, bool withTotal, bool withCorr){

	ostringstream strStream;

	for(unsigned int  i = 0; i < v.size(); i++){

		char type = v[i].forma;

		if (type =='c'){ //case circle

			circle(mat, Point_<int>(v[i].ci, v[i].li), (v[i].d1)/2, Scalar(0), 2);

		} else if (type =='s' || type =='r'){ //case square

			RotatedRect rRect = RotatedRect(Point2f(v[i].ci, v[i].li), Size2f(v[i].d2, v[i].d1), v[i].ang);
			Point2f vertices[4];
			rRect.points(vertices);
			for (int i = 0; i < 4; i++)
				line(mat, vertices[i], vertices[(i+1)%4], Scalar(0), 2);

		} else if (type =='e'){ //case ellipse

			ellipse(mat, Point(v[i].ci, v[i].li), Size((v[i].d2)/2, (v[i].d1)/2), v[i].ang, 0, 360, Scalar(0), 2);

		} else {
			cerr << type << " : Type Invalid" << endl;

		}

		if(withCorr){
			strStream.str("");
			strStream.clear();

			strStream << v[i].r;
			putText(mat, strStream.str(), Point_<int>(v[i].ci-4, v[i].li-4), FONT_HERSHEY_SIMPLEX, 0.7, Scalar(0));
		}
	}

	if(withTotal){
		strStream.str("");
		strStream.clear();
		strStream << v.size();

		putText(mat, strStream.str(), Point_<int>(15, 15), FONT_HERSHEY_SIMPLEX, 0.7, Scalar(0));
	}

}

void printMAXLOCAL(MAXLOCAL m){
	printf("%f %c %d %d %g %g %g\n",m.r,m.forma,m.li,m.ci,m.d1,m.d2,m.ang); 
}

void printTextMAXLOCAL(vector<MAXLOCAL>& v, char* fileName){

	FILE* arq=fopen(fileName,"wt");
	if (arq==NULL) printf("Can not open file %s\n", fileName);
	fprintf(arq,"%lu\n",v.size());
	fprintf(arq,"%8s %1s %4s %4s %8s %8s %3s\n","//correl","f","l","c","d1","d2","deg");
	for (unsigned i=0; i<v.size(); i++) {
		const MAXLOCAL& m=v[i];
		fprintf(arq,"%8.6f %c %4d %4d %8.5f %8.5f %3.0f\n",m.r,m.forma,m.li,m.ci,m.d1,m.d2,m.ang);
	}
	fclose(arq);

}


// KERNEL FUNCTIONS

void createKernelCircle(double diameter, vector<KERNEL>& output, int op){
	KERNEL tempK 	= KERNEL();
	int d			= evenToOdd(sqrt(2.0)*diameter) + BOUND_PIXELS; // TODO: Provide round function
	tempK.shape 	= 'c';
	tempK.d1 		= tempK.d2 = diameter;
	tempK.ang 	= 0;
	tempK.imgF	= Mat_<FLT>(d, d);
	tempK.imgF	= 0;
	tempK.imgG	= cv::Mat_<GRY>(d, d, 128);
	int cW = (tempK.imgF.cols)/2;
	int cH = (tempK.imgF.rows)/2;

	double radius = cvRound(diameter);

	if(op < 1){

		if(radius >= 5){
			circle(tempK.imgF, Point_<int>(cW, cH), ceil(sqrt(2.0)*diameter/2), Scalar(+1.0), CV_FILLED);
			circle(tempK.imgF, Point_<int>(cW, cH), ceil(diameter/2), Scalar(-1.0), CV_FILLED);
		} else if (radius == 4) {
			FLT vkernel[]={ 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 1, 0, 0, 0, 0,
							0, 0, 0, 1,-1, 1, 0, 0, 0,
							0, 0, 1,-1,-1,-1, 1, 0, 0,
							0, 1,-1,-1,-1,-1,-1, 1, 0,
							0, 0, 1,-1,-1,-1, 1, 0, 0,
							0, 0, 0, 1,-1, 1, 0, 0, 0,
							0, 0, 0, 0, 1, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0};
			Mat_<FLT> kernel(9,9,vkernel);
			tempK.imgF = kernel.clone();
		} else if (radius == 3) {
			FLT vkernel[]={ 0, 0, 0, 0, 0, 0, 0,
							0, 1, 1, 1, 1, 1, 0,
							0, 1,-1,-1,-1, 1, 0,
							0, 1,-1,-1,-1, 1, 0,
							0, 1,-1,-1,-1, 1, 0,
							0, 1, 1, 1, 1, 1, 0,
							0, 0, 0, 0, 0, 0, 0};
			Mat_<FLT> kernel(7,7,vkernel);
			tempK.imgF = kernel.clone();
		} else if (radius == 2) {
			FLT vkernel[]={ 0, 0, 1, 0, 0,
							0, 1,-1, 1, 0,
							1,-1,-1,-1, 1,
							0, 1,-1, 1, 0,
							0, 0, 1, 0, 0};
			Mat_<FLT> kernel(5,5,vkernel);
			tempK.imgF = kernel.clone();
		} else if (radius == 1) {
			FLT vkernel[]={ 1, 1, 1,
							1,-1, 1,
							1, 1, 1};
			Mat_<FLT> kernel(3,3,vkernel);
			tempK.imgF = kernel.clone();
		}

		normalizeFLT(tempK.imgF, +1.0, +1.0);
		normalizeFLT(tempK.imgF, -1.0, -1.0);

		//NegPosUm(tempK.imgF);

	}

	if(op == -1 || op == 1){

		if(radius >= 5){
			circle(tempK.imgG, Point_<int>(cW, cH), ceil(sqrt(2.0)*diameter/2), Scalar(255), CV_FILLED);
			circle(tempK.imgG, Point_<int>(cW, cH), ceil(diameter/2), Scalar(0), CV_FILLED);
		} else if (radius == 4) {
			GRY vkernel[]={ 128, 128, 128, 128, 128, 128, 128, 128, 128,
							128, 128, 128, 128, 255, 128, 128, 128, 128,
							128, 128, 128, 255,   0, 255, 128, 128, 128,
							128, 128, 255,   0,   0,   0, 255, 128, 128,
							128, 255,   0,   0,   0,   0,   0, 255, 128,
							128, 128, 255,   0,   0,   0, 255, 128, 128,
							128, 128, 128, 255,   0, 255, 128, 128, 128,
							128, 128, 128, 128, 255, 128, 128, 128, 128,
							128, 128, 128, 128, 128, 128, 128, 128, 128};
			Mat_<GRY> kernel(9,9,vkernel);
			tempK.imgG = kernel.clone();
		} else if (radius == 3) {
			GRY vkernel[]={ 128, 128, 128, 128, 128, 128, 128,
							128, 255, 255, 255, 255, 255, 128,
							128, 255,   0,   0,   0, 255, 128,
							128, 255,   0,   0,   0, 255, 128,
							128, 255,   0,   0,   0, 255, 128,
							128, 255, 255, 255, 255, 255, 128,
							128, 128, 128, 128, 128, 128, 128};
			Mat_<GRY> kernel(7,7,vkernel);
			tempK.imgG = kernel.clone();
		} else if (radius == 2) {
			GRY vkernel[]={ 128, 128, 255, 128, 128,
							128, 255,   0, 255, 128,
							255,   0,   0,   0, 255,
							128, 255,   0, 255, 128,
							128, 128, 255, 128, 128};
			Mat_<GRY> kernel(5,5,vkernel);
			tempK.imgG = kernel.clone();
		} else if(radius == 1){
			GRY vkernel[]={ 255, 255, 255,
							255,   0, 255,
							255, 255, 255};
			Mat_<GRY> kernel(3,3,vkernel);
			tempK.imgG = kernel.clone();
		}
	}

	output.push_back(tempK);
}

void createKernelRectangle(double lengthL, double lengthC, double angle, vector<KERNEL>& output, int op){ // angle in degrees
	KERNEL tempK 	= KERNEL();
	int d1		= evenToOdd(sqrt(2.0)*lengthL)+BOUND_PIXELS; // TODO: Provide round function
	int d2		= evenToOdd(sqrt(2.0)*lengthC)+BOUND_PIXELS; // TODO: Provide round function
	cout << d1 << " " << d2;
	tempK.shape 	= 'r';
	tempK.d1 		= lengthL;
	tempK.d2 		= lengthC;
	tempK.ang 	= angle;
	tempK.imgF	= Mat_<FLT>(d2, d1);
	tempK.imgF	= 0;
	tempK.imgG	= Mat_<GRY>(d2, d1, 128);

	int cW = (tempK.imgF.cols)/2;
	int cH = (tempK.imgF.rows)/2;

	int bL1 = ceil((sqrt(2.0)*lengthL/2)); // bound
	int bL2 = ceil((sqrt(2.0)*lengthC/2)); // bound

	int kL1 = ceil(lengthL/2);
	int kL2 = ceil(lengthC/2);

	if(op < 1){

		rectangle(tempK.imgF, Point_<int>(cW - bL1, cH - bL2), Point_<int>(cW + bL1, cH + bL2), Scalar(+1.0), CV_FILLED);
		rectangle(tempK.imgF, Point_<int>(cW - kL1, cH - kL2), Point_<int>(cW + kL1, cH + kL2), Scalar(-1.0), CV_FILLED);

		normalizeFLT(tempK.imgF, +1.0, +1.0);
		normalizeFLT(tempK.imgF, -1.0, -1.0);
	}

	if(op == -1 || op == 1){
		rectangle(tempK.imgG, Point_<int>(cW - bL1, cH - bL2), Point_<int>(cW + bL1, cH + bL2), Scalar(255), CV_FILLED);
		rectangle(tempK.imgG, Point_<int>(cW - kL1, cH - kL2), Point_<int>(cW + kL1, cH + kL2), Scalar(0), CV_FILLED);
	}

	Point2f src_center(tempK.imgF.cols/2.0F, tempK.imgF.rows/2.0F);
	Mat rot_mat = getRotationMatrix2D(src_center, angle, 1.0);
	Mat dst;

	Rect rect = RotatedRect (src_center, tempK.imgF.size(), angle).boundingRect();

	rot_mat.at<double>(0,2) += rect.width/2.0 - src_center.x;
	rot_mat.at<double>(1,2) += rect.height/2.0 - src_center.y;

	warpAffine(tempK.imgF, dst, rot_mat, rect.size(), INTER_LINEAR, BORDER_CONSTANT, Scalar(0));
	tempK.imgF = dst.clone();
	warpAffine(tempK.imgG, dst, rot_mat, rect.size(), INTER_LINEAR, BORDER_CONSTANT, Scalar(128));
	tempK.imgG = dst.clone();

	output.push_back(tempK);
}

void createKernelSquare(double length, double angle, vector<KERNEL>& output, int op){ // angle in degrees
	if((length/2) > 4)
		createKernelRectangle(length, length, angle, output, op);
	else
		createKernelCircle(length, output, op);

	output[output.size() - 1].shape 	= 's';
}

void createKernelEllipse(double lengthL, double lengthC, double angle, vector<KERNEL>& output, int op){
	KERNEL tempK 	= KERNEL();
	int d		= evenToOdd(sqrt(2.0)*max(lengthL,lengthC))+BOUND_PIXELS; // TODO: Provide round function
	tempK.shape 	= 'e';
	tempK.d1 		= lengthL;
	tempK.d2		= lengthC;
	tempK.ang 	= angle;
	tempK.imgF	= Mat_<FLT>(d, d);
	tempK.imgF	= 0;
	tempK.imgG	= Mat_<GRY>(d, d, 128);

	int cC = (tempK.imgF.cols)/2;
	int cL = (tempK.imgF.rows)/2;

	int bLL = ceil(sqrt(2.0)*lengthL/2); // bound
	int bLC = ceil(sqrt(2.0)*lengthC/2); // bound

	int kLL = ceil(lengthL/2);
	int kLC = ceil(lengthC/2);

	if(op < 1){

		ellipse(tempK.imgF, Point(cC,cL), Size(bLC,bLL), 0, 0, 360, Scalar(+1.0), CV_FILLED);
		ellipse(tempK.imgF, Point(cC,cL), Size(kLC,kLL), 0, 0, 360, Scalar(-1.0), CV_FILLED);

		normalizeFLT(tempK.imgF, +1.0, +1.0);
		normalizeFLT(tempK.imgF, -1.0, -1.0);

	}

	if(op == -1 || op == 1){
		ellipse(tempK.imgG, Point(cC,cL), Size(bLC,bLL), 0, 0, 360, Scalar(255), CV_FILLED);
		ellipse(tempK.imgG, Point(cC,cL), Size(kLC,kLL), 0, 0, 360, Scalar(0), CV_FILLED);
	}

	Point2f src_center(tempK.imgF.cols/2.0F, tempK.imgF.rows/2.0F);
	Mat rot_mat = getRotationMatrix2D(src_center, angle, 1.0);
	Mat dst;

	Rect rect = RotatedRect (src_center, tempK.imgF.size(), angle).boundingRect();

	rot_mat.at<double>(0,2) += rect.width/2.0 - src_center.x;
	rot_mat.at<double>(1,2) += rect.height/2.0 - src_center.y;

	warpAffine(tempK.imgF, dst, rot_mat, rect.size(), INTER_LINEAR, BORDER_CONSTANT, Scalar(0));
	tempK.imgF = dst.clone();
	warpAffine(tempK.imgG, dst, rot_mat, rect.size(), INTER_LINEAR, BORDER_CONSTANT, Scalar(128));
	tempK.imgG = dst.clone();

	output.push_back(tempK);
}

void createKernelGeneric(Mat_<FLT> fltM, Mat_<GRY> gryM, vector<KERNEL>& output, double angle=0){
	KERNEL tempK 	= KERNEL();
	tempK.shape 	= 'g';
	tempK.d1 		= fltM.rows;
	tempK.d2 		= fltM.cols;
	tempK.imgF	= fltM;
	tempK.imgG	= gryM;
	tempK.ang 	= angle;

	output.push_back(tempK);
}

void createKernels(vector<KERNEL> &out, char type, double lengthL, double lengthC, double escOitava, double escMax , double escMin, double angStep ){
	double nOitavas=gLog2(escMax/escMin);
	int nFormas= cvRound(nOitavas*escOitava+1);

	if (type =='c'){ //case circle
		for (int e=0; e<nFormas; e++) {
			double esc=escMax*pow(2, double(-e)/escOitava);

			createKernelCircle(lengthL*esc, out);
		}

	} else if (type =='s'){ //case square
	   	for (int e=0; e<nFormas; e++) {
			for (double a=0; a<90; a += angStep) {
				double esc=escMax*pow(2, double(-e)/escOitava);

				createKernelSquare(lengthL*esc, a, out);

			}
		}

	} else if (type =='r'){ //case rectangle

		for (int e=0; e<nFormas; e++) {
			for (double a=0; a<180; a += angStep) {
				double esc=escMax*pow(2, double(-e)/escOitava);

				createKernelRectangle(lengthL*esc, lengthC*esc, a, out);
			}
		}

	} else if (type =='e'){ //case ellipse

		for (int e=0; e<nFormas; e++) {
			for (double a=0; a<180; a += angStep) {
				double esc=escMax*pow(2, double(-e)/escOitava);

				createKernelEllipse(lengthL*esc, lengthC*esc, a, out);
			}
		}

	} else {
		cerr << type << " : Type Invalid" << endl;

	}
}

// GRANULOMETRY FUNCTIONS

//Correlation

static const cv::Scalar ONE(1);

static const cv::Scalar ZERO(0);

static const cv::Scalar WHITE(255, 255, 255);

// Fourier transform performance is not a monotonic function of a vector
// size - matrices whose dimensions are powers of two are the fastest to
// process, and multiples of 2, 3 and 5 (for example, 300 = 5*5*3*2*2) are
// also processed quite efficiently. Therefore it makes sense to pad input
// data with zeros to get a bit larger matrix that can be transformed much
// faster than the original one.
cv::Size fit(const cv::Size &size)
{
  return cv::Size(cv::getOptimalDFTSize(size.width),
                  cv::getOptimalDFTSize(size.height));
}

cv::Mat F_fwd(const cv::Mat &I, const cv::Size &size)
{
  // Pad input matrix to given size.
  cv::Mat P;
  int m = size.height - I.rows;
  int n = size.width - I.cols;
  cv::copyMakeBorder(I, P, 0, m, 0, n, cv::BORDER_CONSTANT, ZERO);

  // Compute Fourier transform for input data. The last argument
  // informs the dft() function of how many non-zero rows are there,
  // so it can handle the rest of the rows more efficiently and save
  // some time.
  cv::Mat F;
  cv::dft(P, F, 0, I.rows);

  return F;
}

cv::Mat F_inv(const cv::Mat &F, const cv::Size &size)
{
  // Compute inverse Fourier transform for input data. The last
  // argument informs the dft() function of how many non-zero
  // rows are expected in the output, so it can handle the rest
  // of the rows more efficiently and save some time.
  cv::Mat I;
  cv::dft(F, I, cv::DFT_INVERSE + cv::DFT_SCALE, size.height);

  return I(cv::Rect(0, 0, size.width, size.height));
}

Mat fF_I;

cv::Mat C(const cv::Mat &T, const cv::Mat &I, const cv::Size &size)
{
  // Compute the Fourier transforms of template and image.
  cv::Mat F_T = F_fwd(T, size);
  cv::Mat F_I = F_fwd(I, size);

  // Compute the cross correlation in the frequency domain.
  cv::Mat F_TI;
  cv::mulSpectrums(F_I, F_T, F_TI, 0, true);

  // Compute the inverse Fourier transform of the cross-correlation,
  // dismissing those rows and columns of the cross-correlation
  // matrix that would require the template to "roll over" the image.
  cv::Size clipped;
  clipped.width = I.cols - T.cols;
  clipped.height = I.rows - T.rows;
  return F_inv(F_TI, clipped);
}

cv::Mat FC(const cv::Mat &I, const cv::Mat &T, const cv::Size &size)
{
  // Compute the Fourier transforms of template and image.
  cv::Mat F_T = F_fwd(T, size);

  // Compute the cross correlation in the frequency domain.
  cv::Mat F_TI;
  cv::mulSpectrums(fF_I, F_T, F_TI, 0, true);

  // Compute the inverse Fourier transform of the cross-correlation,
  // dismissing those rows and columns of the cross-correlation
  // matrix that would require the template to "roll over" the image.
  cv::Size clipped;
  clipped.width = I.cols - T.cols;
  clipped.height = I.rows - T.rows;
  return F_inv(F_TI, clipped);
}


cv::Mat W(const cv::Mat &T)
{
  return cv::Mat(T.size(), CV_64F, ONE);
}

cv::Point3f matchTemplate(const cv::Mat &T, const cv::Mat &I, cv::Mat &out)
{
  // Compute the optimal size for DFT computing.
  cv::Size size = fit(I.size());

  //Compute the cross-correlation and normalizing matrix.
  cv::Mat C_TI = C(T, I, size);
  cv::Mat M_I = C(W(T), I.mul(I), size);

  int i_s, j_s;
  float r = 0;
  int rows = C_TI.rows;
  int cols = C_TI.cols;
  /*for (int i = 0; i < rows; i++)
  {
    for (int j = 0; j < cols; j++)
    {
      float v = C_TI.at<FLT>(i, j) / sqrt(M_I.at<FLT>(i, j));
      if (r < v)
      {
        r = v;
        i_s = i;
        j_s = j;
      }
    }
  }*/

  out = C_TI.clone();

  return cv::Point3f(j_s, i_s, r);
}

cv::Mat L(const cv::Mat &I)
{
  cv::Mat L_I, L_F;
  cv::cvtColor(I, L_I, CV_BGR2GRAY);
  L_I.convertTo(L_F, CV_64F);
  return L_F;
}

//Correlation

void correlation(Mat in, Mat &out, KERNEL ker, int type){
	Mat inF;

	in.convertTo(inF, CV_32F, 1.0/255.0, 0.0);

	if(in.rows < ker.imgF.rows || in.cols < ker.imgF.cols){
		out = Scalar(0);

		return;
	}

	matchTemplate(inF, ker.imgF, out, type);

	//filter2D(inF, out, CV_32F, ker.imgF, Point(-1,-1), 0, BORDER_REPLICATE);
	//matchTemplate( inF, ker.imgF, out, CV_TM_SQDIFF );
	//matchTemplate( inF, ker.imgF, out, CV_TM_SQDIFF_NORMED );
	//matchTemplate( inF, ker.imgF, out, CV_TM_CCORR );
	//matchTemplate( inF, ker.imgF, out, CV_TM_CCORR_NORMED );
	//matchTemplate( inF, ker.imgF, out, CV_TM_CCOEFF );
	//matchTemplate( inF, ker.imgF, out, CV_TM_CCOEFF_NORMED );
}

void correlationFFTInBatch(Mat& img, vector<MAXLOCAL> &out, vector<KERNEL> kers, float minCorr, int maxDist){

	bool isMaxExtremal = false;
	KERNEL tmp;
	Mat_<FLT> corrMat;

	Mat tF;
	Mat inF; 
	convertGRYToFLT(img, inF);
	inF.convertTo(inF, CV_64F);

	cv::Size size = fit(inF.size());

	fF_I = F_fwd(inF, size);

	std::vector<MAXLOCAL> v;

	for(unsigned int i = 0; i < kers.size(); i++){
		tmp = kers[i];

		printf("Applying filter #%u-#%lu\r", i+1, kers.size());
		fflush(stdout); 
		
		tmp.imgF.convertTo(tF, CV_64F);
		corrMat = FC(inF, tF, size);

		for (int l=1; l<corrMat.rows-1; l++){
			for (int c=1; c<corrMat.cols-1; c++) {

				FLT valor = corrMat.at<FLT>(l,c);

				if (valor>=minCorr) {
					isMaxExtremal = true;

					for (int l2=-1; l2<=1 && isMaxExtremal; l2++)
						for (int c2=-1; c2<=1 && isMaxExtremal; c2++)
							if (valor<corrMat(l+l2,c+c2)) {
								isMaxExtremal = false;
							}

							if(isMaxExtremal){

								int li,ci;
								li = l + (img.rows - corrMat.rows)/2;
								ci = c + (img.cols - corrMat.cols)/2;
								out.push_back(MAXLOCAL(valor, tmp.shape, li, ci, tmp.d1, tmp.d2, tmp.ang));
							}

				}
			}
		}
	}

	removeCloser(out, out, maxDist);
}

void correlationInBatch(Mat& img, vector<MAXLOCAL> &out, vector<KERNEL> kers, float minCorr, int maxDist, int type){

#pragma omp parallel shared(out)
{
	bool isMaxExtremal = false;
	KERNEL tmp;
	Mat_<FLT> corrMat;

	std::vector<MAXLOCAL> v;

	int tN = omp_get_thread_num();
	int T  = omp_get_num_threads();

	uint work = (kers.size() / T);

	uint end;

	if(tN == T-1){
		end = kers.size();
	}else{
		end = work + work*tN;
	}

	for(unsigned int i = work*tN; i < end; i++){
		tmp = kers[i];

		printf("Applying filter #%u-#%lu\r", i+1, kers.size());
		fflush(stdout); 

		correlation(img, corrMat, tmp, type);

		for (int l=1; l<corrMat.rows-1; l++){
			for (int c=1; c<corrMat.cols-1; c++) {

				FLT valor = corrMat.at<FLT>(l,c);

				if (valor>=minCorr) {
					isMaxExtremal = true;

					for (int l2=-1; l2<=1 && isMaxExtremal; l2++)
						for (int c2=-1; c2<=1 && isMaxExtremal; c2++)
							if (valor<corrMat(l+l2,c+c2)) {
								isMaxExtremal = false;
							}

							if(isMaxExtremal){

								int li,ci;
								li = l + (img.rows - corrMat.rows)/2;
								ci = c + (img.cols - corrMat.cols)/2;
								v.push_back(MAXLOCAL(valor, tmp.shape, li, ci, tmp.d1, tmp.d2, tmp.ang));
							}

				}
			}
		}
	}

	for(uint i = 0; i < v.size(); i++){
		out.push_back(v[i]);
	}
}
	removeCloser(out, out, maxDist);
}

void computeIntersection(MATRIZ<MAXLOCAL>& n_out, MAXLOCAL &loc, float maxInt){
	std::vector<MAXLOCAL> others;

	RotatedRect rRect = RotatedRect(Point2f(loc.ci, loc.li), Size2f(loc.d2, loc.d1), loc.ang);
	Rect brect = rRect.boundingRect();

	uint row_orig = brect.y;
	uint col_orig = brect.x;

	uint row_dest = row_orig + brect.height;
	uint col_dest = col_orig + brect.width;

	if (row_orig < 0) row_orig = 0;
	if (col_orig < 0) col_orig = 0;

	if (row_dest > n_out.rows) row_dest = n_out.rows;
	if (col_dest > n_out.cols) col_dest = n_out.cols;

	double intersecArea = 0;
	int area 		 = 0;

	for(uint r = row_orig; r < row_dest; r++){
		for(uint c = col_orig; c < col_dest; c++){

			if(!isInsideKernel(r, c, loc)){
				continue;
			}

			area++;

			if(n_out(r, c).r >= -1 && !isMAXEqual(n_out(r, c), loc)){
				uint s;

				for (s = 0; s < others.size(); s++){
					if(isMAXEqual(n_out(r, c), others[s]))
						break;
				}

				if(s >= others.size()){
					others.push_back(n_out(r, c));
					others[s].area_int = 1;
				}
				else{
					others[s].area_int++;
				}

				if(others[s].area_int/others[s].area > maxInt){
					r = row_dest;
					c = col_dest;
					intersecArea = area;
					loc.valid = false;
				}

				intersecArea++;
			}
		}
	}

	loc.area = area;
	loc.area_int = intersecArea;

	intersecArea = intersecArea/area;


	if(intersecArea <= maxInt){
		for(uint r = row_orig; r < row_dest; r++){
			for(uint c = col_orig; c < col_dest; c++){
				n_out(r, c) = loc;
			}
		}
	}
	else{
		loc.valid = false;
	}

}

void sift(vector<MAXLOCAL> in, vector<MAXLOCAL>& out, float minCorr, float maxInter, Size size){
	if(in.size() <= 1)
		return;
	
	MATRIZ<MAXLOCAL> mAux(size.height, size.width);
	global = Mat_<GRY>(size); global = 0;
	std::vector<MAXLOCAL> aux;

	MAXLOCAL m = MAXLOCAL(-10, 'c', 0, 0, 0, 0, 0);

	for(uint l = 0; l < mAux.rows; l++){
		for(uint c = 0; c < mAux.cols; c++){
			mAux(l, c) = m;
		}
	}

	sortMAXLOCAL(in, 0, in.size()-1);

	for(uint i = 0; i < in.size(); i++){

		in[i].area = 0;
		in[i].area_int = 0;

		if(in[i].r >= minCorr){
			global = 0;

			computeIntersection(mAux, in[i], maxInter);

			if(in[i].valid){

				aux.push_back(in[i]);
			}
		}
		else
			i = in.size();
	}

	out = aux;
}

void sift(vector<MAXLOCAL> in, vector<MAXLOCAL>& out, double minCorrCir, double maxInterCir, double minCorrRet, double maxInterRet,
              double minCorrQua, double maxInterQua, double minCorrEli, double maxInterEli, Size size){
	if(in.size() <= 1)
		return;
	
	MATRIZ<MAXLOCAL> mAux(size.height, size.width);
	global = Mat_<GRY>(size); global = 0;
	std::vector<MAXLOCAL> aux;

	MAXLOCAL m = MAXLOCAL(-10, 'c', 0, 0, 0, 0, 0);

	for(uint l = 0; l < mAux.rows; l++){
		for(uint c = 0; c < mAux.cols; c++){
			mAux(l, c) = m;
		}
	}

	sortMAXLOCAL(in, 0, in.size()-1);

	for(uint i = 0; i < in.size(); i++){

		char type = in[i].forma;

		in[i].area = 0;
		in[i].area_int = 0;

		if((type == 'c' && in[i].r >= minCorrCir) ||
			(type == 's' && in[i].r >= minCorrQua) ||
			(type == 'r' && in[i].r >= minCorrRet) ||
			(type == 'e' && in[i].r >= minCorrEli)){

			global = 0;

			double maxInter;

			if (type =='c'){ //case circle
				maxInter = maxInterCir;

			} else if (type =='s'){ //case square
				maxInter = maxInterQua;

			} else if (type =='r'){ //case rectangle
				maxInter = maxInterRet;

			} else if (type =='e'){ //case ellipse
				maxInter = maxInterEli;
				
			} 

			computeIntersection(mAux, in[i], maxInter);

			if(in[i].valid){

				aux.push_back(in[i]);
			}
		}
		else
			i = in.size();
	}

	out = aux;
}

Mat_<COR> mostrar_MSER;

void putMAXLOCALOnMatrix(MATRIZ<vector<DBL>>& n_out, MAXLOCAL &loc, DBL value){

	RotatedRect rRect = RotatedRect(Point2f(loc.ci, loc.li), Size2f(loc.d2, loc.d1), loc.ang);
	Rect brect = rRect.boundingRect();

	uint row_orig = brect.y;
	uint col_orig = brect.x;

	uint row_dest = row_orig + brect.height;
	uint col_dest = col_orig + brect.width;

	if (row_orig < 0) row_orig = 0;
	if (col_orig < 0) col_orig = 0;

	if (row_dest > n_out.rows) row_dest = n_out.rows;
	if (col_dest > n_out.cols) col_dest = n_out.cols;

	loc.area = 0;

	for(uint r = row_orig; r < row_dest; r++){
		for(uint c = col_orig; c < col_dest; c++){

			if(!isInsideKernel(r, c, loc)){
				continue;
			}

			loc.area++;

			n_out(r, c).push_back(value);

			mostrar_MSER.at<COR>(r ,c) = COR(255,0,0);
		}
	}

}

void siftMSER(Mat img, vector<MAXLOCAL> in, vector<MAXLOCAL>& out, double minCorrMSER){
	if(in.size() <= 1)
		return;

	Mat_<GRY> grayImg;

	cv::cvtColor(img, grayImg, CV_BGR2GRAY);

	MATRIZ<vector<DBL>> mAux_G(img.rows, img.cols); //Aux matrix to Granulometry

	DBL minLocArea, maxLocArea;

	mostrar_MSER = img.clone();

	minLocArea = std::numeric_limits<double>::infinity();
	maxLocArea = 0;

	for (uint i = 0; i < in.size(); i++) {

		putMAXLOCALOnMatrix(mAux_G, in[i], i);
		in[i].valid = true;

		if(in[i].area > maxLocArea)
			maxLocArea = in[i].area;
		if(in[i].area < minLocArea)
			minLocArea = in[i].area;
	}

	vector<vector<Point>> regions;
	cv::MSER descriptor = cv::MSER(1, (int) minLocArea/2, (int) maxLocArea*2);
	//cv::MSER descriptor = cv::MSER();
	descriptor(grayImg, regions);

	MATRIZ<vector<int>> mAux_M(img.rows, img.cols); //Aux matrix to MSER
	std::vector<MAXLOCAL> aux;

	sortMAXLOCAL(in, 0, in.size()-1);
	
	for (uint i = 0; i < regions.size(); i++) {
		for (uint j = 0; j < regions[i].size(); j++) {

			mAux_M(regions[i][j].y, regions[i][j].x).push_back(i);
		}

	}


	MATRIZ<int> relation(regions.size(), in.size());

	for (int r = 0; r < img.rows; r++) {
		for (int c = 0; c < img.cols; c++) {

			for(uint i=0; i < mAux_M(r, c).size(); i++){
				int region = mAux_M(r, c)[i];

				for(uint j=0; j < mAux_G(r, c).size(); j++){
					int grain = mAux_G(r, c)[j];

					relation(region, grain)++;
				}
				
			}
		}
	}

	for (uint r = 0; r < relation.rows; r++) {
		double totalAreaGrain = 0;
		double totalAreaInter = 0; 
		double totalAreaMSER = regions[r].size(); 
		for (uint c = 0; c < relation.cols; c++) {
			if(relation(r,c) > 0){
				totalAreaGrain += in[c].area;
				totalAreaInter += relation(r,c);
			}
		}

		double ratio = totalAreaInter / (totalAreaMSER + totalAreaGrain - totalAreaInter);
		
		if(ratio >= minCorrMSER)
		{
			for (uint c = 0; c < relation.cols; c++) {
				if(relation(r, c) > 0 && in[c].valid){
					in[c].valid = false;
					aux.push_back(in[c]);
				}
			}
		}
	}

	out = aux;
}

void chooseKmeans(Mat img, vector<MAXLOCAL> in, vector<MAXLOCAL>& out, double num_k){

	if(in.size() <= 1)
		return;

	namedWindow("Amostra S(Considerar)", WINDOW_NORMAL);

	vector<vector<FLT>> msamp;

	Mat_<COR> ent_cor = img.clone();

	vector<MAXLOCAL> aux;

	for (int i = 0; i < in.size(); i++) {

		queue<int> fila;
		fila.push(in[i].li); fila.push(in[i].ci);
		int max_distancia = (int)in[i].d1 / 2;
		max_distancia++;

		for (int l = in[i].li - in[i].d1; l < in[i].li + in[i].d1; l++) {
			for (int c = in[i].ci - in[i].d2; c < in[i].ci + in[i].d2; c++) {

				double distancia = sqrt(pow(in[i].li - l, 2) + pow(in[i].ci - c, 2));
				abs(distancia);

				if (distancia < max_distancia) {
					msamp.push_back(vector<FLT>(5));

					msamp[msamp.size() - 1][0] = ent_cor(l, c)[0];
					msamp[msamp.size() - 1][1] = ent_cor(l, c)[1];
					msamp[msamp.size() - 1][2] = ent_cor(l, c)[2];
					msamp[msamp.size() - 1][3] = in[0].r / in[i].r;
					msamp[msamp.size() - 1][4] = i;
				}
			}
		}


	}

	Mat_<FLT> kmeansSample(msamp.size(), 3);

	for (int i = 0;i < kmeansSample.rows; i++) {
		for (int j = 0; j < kmeansSample.cols; j++) {

			kmeansSample(i, j) = msamp[i][j];

		}
	}

	int clusterCount = num_k;
	Mat labels;
	int attempts = 1;
	Mat centers;

	kmeans(kmeansSample, clusterCount, labels, TermCriteria(CV_TERMCRIT_ITER | CV_TERMCRIT_EPS, 1000000, 0.0001), attempts, KMEANS_PP_CENTERS, centers);

	vector<int> v_cluster(clusterCount, 0);
	int a_id = 0;

	for (int i = 0; i < in.size(); i++) {

		vector<FLT> sample_vet(clusterCount, 0);

		for (int l = 0; l < kmeansSample.rows; l++) {
			for (int c = 0; c < kmeansSample.cols; c++) {
				if (msamp[l][4] == i) {
					int cluster_idx = labels.at<int>(l, 0);
					sample_vet[cluster_idx]++;
				}
			}
		}

		a_id = 0;

		for (int c = 1; c < sample_vet.size(); c++) {
			if (sample_vet[a_id] <= sample_vet[c]) {
				a_id = c;
			}
		}

		if (v_cluster[a_id] == 0) {

			Mat_<COR> tcor;
			tcor = ent_cor(Rect(in[i].ci - 5 - in[i].d2 / 2, in[i].li - 5  - in[i].d1 / 2, in[i].d2 + 5 , in[i].d1 + 5 )).clone();

			cv::resize(tcor, tcor, Size(80, 80));
			
			imshow("Amostra S(Considerar)", tcor);
			char press = waitKey(0);

			if (press == 'S' || press == 's') {
				v_cluster[a_id] = 1;
			}
			else {
				v_cluster[a_id] = -1;
			}
		}

		if (v_cluster[a_id] == 1) aux.push_back(in[i]);

	}

	out.clear();

	for(uint i=0; i < aux.size(); i++){
		out.push_back(aux[i]);
	}

}
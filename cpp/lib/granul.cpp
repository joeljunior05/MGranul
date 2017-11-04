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

//CrossCorrelation FFT

void crossCorr( const Mat& img, const Mat& _templ, Mat& corr,
                Size corrsize, int ctype,
                Point anchor, double delta, int borderType )
{
    const double blockScale = 4.5;
    const int minBlockSize = 256;
    std::vector<uchar> buf;

    Mat templ = _templ;
    int depth = img.depth(), cn = img.channels();
    int tdepth = templ.depth(), tcn = templ.channels();
    int cdepth = CV_MAT_DEPTH(ctype), ccn = CV_MAT_CN(ctype);

    CV_Assert( img.dims <= 2 && templ.dims <= 2 && corr.dims <= 2 );

    if( depth != tdepth && tdepth != std::max(CV_32F, depth) )
    {
        _templ.convertTo(templ, std::max(CV_32F, depth));
        tdepth = templ.depth();
    }

    CV_Assert( depth == tdepth || tdepth == CV_32F);
    CV_Assert( corrsize.height <= img.rows + templ.rows - 1 &&
               corrsize.width <= img.cols + templ.cols - 1 );

    CV_Assert( ccn == 1 || delta == 0 );

    corr.create(corrsize, ctype);

    int maxDepth = depth > CV_8S ? CV_64F : std::max(std::max(CV_32F, tdepth), cdepth);
    Size blocksize, dftsize;

    blocksize.width = cvRound(templ.cols*blockScale);
    blocksize.width = std::max( blocksize.width, minBlockSize - templ.cols + 1 );
    blocksize.width = std::min( blocksize.width, corr.cols );
    blocksize.height = cvRound(templ.rows*blockScale);
    blocksize.height = std::max( blocksize.height, minBlockSize - templ.rows + 1 );
    blocksize.height = std::min( blocksize.height, corr.rows );

    dftsize.width = std::max(getOptimalDFTSize(blocksize.width + templ.cols - 1), 2);
    dftsize.height = getOptimalDFTSize(blocksize.height + templ.rows - 1);
    if( dftsize.width <= 0 || dftsize.height <= 0 )
        CV_Error( CV_StsOutOfRange, "the input arrays are too big" );

    // recompute block size
    blocksize.width = dftsize.width - templ.cols + 1;
    blocksize.width = MIN( blocksize.width, corr.cols );
    blocksize.height = dftsize.height - templ.rows + 1;
    blocksize.height = MIN( blocksize.height, corr.rows );

    Mat dftTempl( dftsize.height*tcn, dftsize.width, maxDepth );
    Mat dftImg( dftsize, maxDepth );

    int i, k, bufSize = 0;
    if( tcn > 1 && tdepth != maxDepth )
        bufSize = templ.cols*templ.rows*CV_ELEM_SIZE(tdepth);

    if( cn > 1 && depth != maxDepth )
        bufSize = std::max( bufSize, (blocksize.width + templ.cols - 1)*
            (blocksize.height + templ.rows - 1)*CV_ELEM_SIZE(depth));

    if( (ccn > 1 || cn > 1) && cdepth != maxDepth )
        bufSize = std::max( bufSize, blocksize.width*blocksize.height*CV_ELEM_SIZE(cdepth));

    buf.resize(bufSize);

    // compute DFT of each template plane
    for( k = 0; k < tcn; k++ )
    {
        int yofs = k*dftsize.height;
        Mat src = templ;
        Mat dst(dftTempl, Rect(0, yofs, dftsize.width, dftsize.height));
        Mat dst1(dftTempl, Rect(0, yofs, templ.cols, templ.rows));

        if( tcn > 1 )
        {
            src = tdepth == maxDepth ? dst1 : Mat(templ.size(), tdepth, &buf[0]);
            int pairs[] = {k, 0};
            mixChannels(&templ, 1, &src, 1, pairs, 1);
        }

        if( dst1.data != src.data )
            src.convertTo(dst1, dst1.depth());

        if( dst.cols > templ.cols )
        {
            Mat part(dst, Range(0, templ.rows), Range(templ.cols, dst.cols));
            part = Scalar::all(0);
        }
        dft(dst, dst, 0, templ.rows);
    }

    int tileCountX = (corr.cols + blocksize.width - 1)/blocksize.width;
    int tileCountY = (corr.rows + blocksize.height - 1)/blocksize.height;
    int tileCount = tileCountX * tileCountY;

    Size wholeSize = img.size();
    Point roiofs(0,0);
    Mat img0 = img;

    if( !(borderType & BORDER_ISOLATED) )
    {
        img.locateROI(wholeSize, roiofs);
        img0.adjustROI(roiofs.y, wholeSize.height-img.rows-roiofs.y,
                       roiofs.x, wholeSize.width-img.cols-roiofs.x);
    }
    borderType |= BORDER_ISOLATED;

    // calculate correlation by blocks
    for( i = 0; i < tileCount; i++ )
    {
        int x = (i%tileCountX)*blocksize.width;
        int y = (i/tileCountX)*blocksize.height;

        Size bsz(std::min(blocksize.width, corr.cols - x),
                 std::min(blocksize.height, corr.rows - y));
        Size dsz(bsz.width + templ.cols - 1, bsz.height + templ.rows - 1);
        int x0 = x - anchor.x + roiofs.x, y0 = y - anchor.y + roiofs.y;
        int x1 = std::max(0, x0), y1 = std::max(0, y0);
        int x2 = std::min(img0.cols, x0 + dsz.width);
        int y2 = std::min(img0.rows, y0 + dsz.height);
        Mat src0(img0, Range(y1, y2), Range(x1, x2));
        Mat dst(dftImg, Rect(0, 0, dsz.width, dsz.height));
        Mat dst1(dftImg, Rect(x1-x0, y1-y0, x2-x1, y2-y1));
        Mat cdst(corr, Rect(x, y, bsz.width, bsz.height));

        for( k = 0; k < cn; k++ )
        {
            Mat src = src0;
            dftImg = Scalar::all(0);

            if( cn > 1 )
            {
                src = depth == maxDepth ? dst1 : Mat(y2-y1, x2-x1, depth, &buf[0]);
                int pairs[] = {k, 0};
                mixChannels(&src0, 1, &src, 1, pairs, 1);
            }

            if( dst1.data != src.data )
                src.convertTo(dst1, dst1.depth());

            if( x2 - x1 < dsz.width || y2 - y1 < dsz.height )
                copyMakeBorder(dst1, dst, y1-y0, dst.rows-dst1.rows-(y1-y0),
                               x1-x0, dst.cols-dst1.cols-(x1-x0), borderType);

            dft( dftImg, dftImg, 0, dsz.height );
            Mat dftTempl1(dftTempl, Rect(0, tcn > 1 ? k*dftsize.height : 0,
                                         dftsize.width, dftsize.height));
            mulSpectrums(dftImg, dftTempl1, dftImg, 0, true);
            dft( dftImg, dftImg, DFT_INVERSE + DFT_SCALE, bsz.height );

            src = dftImg(Rect(0, 0, bsz.width, bsz.height));

            if( ccn > 1 )
            {
                if( cdepth != maxDepth )
                {
                    Mat plane(bsz, cdepth, &buf[0]);
                    src.convertTo(plane, cdepth, 1, delta);
                    src = plane;
                }
                int pairs[] = {0, k};
                mixChannels(&src, 1, &cdst, 1, pairs, 1);
            }
            else
            {
                if( k == 0 )
                    src.convertTo(cdst, cdepth, 1, delta);
                else
                {
                    if( maxDepth != cdepth )
                    {
                        Mat plane(bsz, cdepth, &buf[0]);
                        src.convertTo(plane, cdepth);
                        src = plane;
                    }
                    add(src, cdst, cdst);
                }
            }
        }
    }
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

	std::vector<MAXLOCAL> v;

	for(unsigned int i = 0; i < kers.size(); i++){
		tmp = kers[i];

		printf("Applying filter #%u-#%lu\r", i+1, kers.size());
		fflush(stdout); 

		Size corrSize(img.cols - tmp.imgF.cols + 1, img.rows - tmp.imgF.rows + 1);
		corrMat.create(corrSize);

		crossCorr( img, tmp.imgF, corrMat, corrMat.size(), corrMat.type(), Point(0,0), 0, 0);

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
		int T = omp_get_num_threads();

		uint work = (kers.size() / T);

		uint end;

		if (tN == T - 1)
		{
			end = kers.size();
		}
		else
		{
			end = work + work * tN;
		}

		for (unsigned int i = work * tN; i < end; i++)
		{
			tmp = kers[i];

			printf("Applying filter #%u-#%lu\r", i + 1, kers.size());
			fflush(stdout);

			correlation(img, corrMat, tmp, type);

			for (int l = 1; l < corrMat.rows - 1; l++)
			{
				for (int c = 1; c < corrMat.cols - 1; c++)
				{

					FLT valor = corrMat.at<FLT>(l, c);

					if (valor >= minCorr)
					{
						isMaxExtremal = true;

						for (int l2 = -1; l2 <= 1 && isMaxExtremal; l2++)
							for (int c2 = -1; c2 <= 1 && isMaxExtremal; c2++)
								if (valor < corrMat(l + l2, c + c2))
								{
									isMaxExtremal = false;
								}

						if (isMaxExtremal)
						{

							int li, ci;
							li = l + (img.rows - corrMat.rows) / 2;
							ci = c + (img.cols - corrMat.cols) / 2;
							v.push_back(MAXLOCAL(valor, tmp.shape, li, ci, tmp.d1, tmp.d2, tmp.ang));
						}
					}
				}
			}
		}

		for (uint i = 0; i < v.size(); i++)
		{
			out.push_back(v[i]);
		}
}

removeCloser(out, out, maxDist);
}

void computeIntersection(MATRIZ<MAXLOCAL>& n_out, MAXLOCAL &loc, float maxInt){
	std::vector<MAXLOCAL> others;

	RotatedRect rRect = RotatedRect(Point2f(loc.ci, loc.li), Size2f(loc.d2, loc.d1), loc.ang);
	Rect brect = rRect.boundingRect();

	int row_orig = brect.y;
	int col_orig = brect.x;

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
				if(isInsideKernel(r, c, loc))
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

void putMAXLOCALOnMatrix(MATRIZ<vector<DBL>>& n_out, MAXLOCAL &loc, DBL value){

	RotatedRect rRect = RotatedRect(Point2f(loc.ci, loc.li), Size2f(loc.d2, loc.d1), loc.ang);
	Rect brect = rRect.boundingRect();

	int row_orig = brect.y;
	int col_orig = brect.x;

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

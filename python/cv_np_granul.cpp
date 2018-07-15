#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include "granultypes.h"
#include <numpy/ndarraytypes.h>
#include <numpy/ndarrayobject.h>
#include <granul.h>

#define NUMPY_IMPORT_ARRAY_RETVAL

/***** Functions copied from OpenCV source in order to provide cv <-> numpy capabilities ******/

static PyObject* opencv_error = 0;

void init_numpy()
{
    import_array();
}

static size_t REFCOUNT_OFFSET = (size_t)&(((PyObject*)0)->ob_refcnt) +
    (0x12345678 != *(const size_t*)"\x78\x56\x34\x12\0\0\0\0\0")*sizeof(int);

static inline PyObject* pyObjectFromRefcount(const int* refcount)
{
    return (PyObject*)((size_t)refcount - REFCOUNT_OFFSET);
}

static inline int* refcountFromPyObject(const PyObject* obj)
{
    return (int*)((size_t)obj + REFCOUNT_OFFSET);
}

static inline int* refcountFromPyObject(const PyArrayObject* obj)
{
    return refcountFromPyObject((const PyObject*) obj);
}

static int failmsg(const char *fmt, ...)
{
    char str[1000];

    va_list ap;
    va_start(ap, fmt);
    vsnprintf(str, sizeof(str), fmt, ap);
    va_end(ap);

    PyErr_SetString(PyExc_TypeError, str);
    return 0;
}

class PyAllowThreads
{
public:
    PyAllowThreads() : _state(PyEval_SaveThread()) {}
    ~PyAllowThreads()
    {
        PyEval_RestoreThread(_state);
    }
private:
    PyThreadState* _state;
};


#define ERRWRAP2(expr) \
try \
{ \
    PyAllowThreads allowThreads; \
    expr; \
} \
catch (const cv::Exception &e) \
{ \
    PyErr_SetString(opencv_error, e.what()); \
    return 0; \
}

class PyEnsureGIL
{
public:
    PyEnsureGIL() : _state(PyGILState_Ensure()) {}
    ~PyEnsureGIL()
    {
        PyGILState_Release(_state);
    }
private:
    PyGILState_STATE _state;
};

class NumpyAllocator : public MatAllocator
{
public:
    NumpyAllocator() {}
    ~NumpyAllocator() {}

    void allocate(int dims, const int* sizes, int type, int*& refcount,
                  uchar*& datastart, uchar*& data, size_t* step)
    {
        PyEnsureGIL gil;

        int depth = CV_MAT_DEPTH(type);
        int cn = CV_MAT_CN(type);
        const int f = (int)(sizeof(size_t)/8);
        int typenum = depth == CV_8U ? NPY_UBYTE : depth == CV_8S ? NPY_BYTE :
                      depth == CV_16U ? NPY_USHORT : depth == CV_16S ? NPY_SHORT :
                      depth == CV_32S ? NPY_INT : depth == CV_32F ? NPY_FLOAT :
                      depth == CV_64F ? NPY_DOUBLE : f*NPY_ULONGLONG + (f^1)*NPY_UINT;
        int i;
        npy_intp _sizes[CV_MAX_DIM+1];
        
        for( i = 0; i < dims; i++ ){
            _sizes[i] = sizes[i];
        }
        if( cn > 1 )
        {
            /*if( _sizes[dims-1] == 1 )
                _sizes[dims-1] = cn;
            else*/
                _sizes[dims++] = cn;
        }

        PyObject* o = PyArray_SimpleNew(dims, _sizes, typenum);

        if(!o)
            CV_Error_(CV_StsError, ("The numpy array of typenum=%d, ndims=%d can not be created", typenum, dims));
        refcount = refcountFromPyObject(o);
        npy_intp* _strides = PyArray_STRIDES((PyArrayObject*) o);
        for( i = 0; i < dims - (cn > 1); i++ )
            step[i] = (size_t)_strides[i];
        datastart = data = (uchar*)PyArray_DATA((PyArrayObject*) o);
    }

    void deallocate(int* refcount, uchar*, uchar*)
    {
        PyEnsureGIL gil;
        if( !refcount )
            return;
        PyObject* o = pyObjectFromRefcount(refcount);
        Py_INCREF(o);
        Py_DECREF(o);
    }
};

NumpyAllocator g_numpyAllocator;

PyObject* cvMatToPyObject(Mat& m){
    if( !m.data )
        Py_RETURN_NONE;
    cv::Mat temp = cv::Mat();
    cv::Mat *p = (Mat*)&m;
    if(!p->refcount || p->allocator != &g_numpyAllocator)
    {
        temp.allocator = &g_numpyAllocator;
        ERRWRAP2(m.copyTo(temp));
        p = &temp;
    }
    p->addref();
    return pyObjectFromRefcount(p->refcount);
}

cv::Mat PyArrayObjTocvMat(PyArrayObject* o){

    int typenum = PyArray_TYPE(o);
    int type =  typenum == NPY_UBYTE    ? CV_8U : 
                typenum == NPY_BYTE     ? CV_8S :
                typenum == NPY_USHORT   ? CV_16U : 
                typenum == NPY_SHORT    ? CV_16S :
                typenum == NPY_INT      ? CV_32S :
                typenum == NPY_LONG     ? CV_32S :
                typenum == NPY_FLOAT    ? CV_32F :
                typenum == NPY_DOUBLE   ? CV_64F : -1;

    if( type < 0 )
    {
        failmsg("toMat: Data type = %d is not supported", typenum);
        return cv::Mat();
    }

    int ndims = PyArray_NDIM(o);

    if(ndims >= CV_MAX_DIM)
    {
        failmsg("toMat: Dimensionality (=%d) is too high", ndims);
        return cv::Mat();
    }

    int size[CV_MAX_DIM+1];
    size_t step[CV_MAX_DIM+1], elemsize = CV_ELEM_SIZE1(type);
    const npy_intp* _sizes = PyArray_DIMS(o);
    const npy_intp* _strides = PyArray_STRIDES(o);

    for(int i = 0; i < ndims; i++)
    {
        size[i] = (int)_sizes[i];
        step[i] = (size_t)_strides[i];
    }

    bool transposed = false;
    if( ndims >= 2 && step[0] < step[1] )
    {
        std::swap(size[0], size[1]);
        std::swap(step[0], step[1]);
        transposed = true;
    }

    if( ndims == 3 && size[2] <= CV_CN_MAX && step[1] == elemsize*size[2] )
    {
        ndims--;
        type |= CV_MAKETYPE(0, size[2]);
    }

    if( ndims > 2)
    {
        failmsg("toMat: Object has more than 2 dimensions");
    }

    cv::Mat ret = Mat(ndims, size, type, PyArray_DATA(o), step);

    if( ret.data )
    {
        ret.refcount = refcountFromPyObject(o);
        ret.addref(); // protect the original numpy array from deallocation
                    // (since Mat destructor will decrement the reference counter)
    };
    ret.allocator = &g_numpyAllocator;

    if( transposed )
    {
        Mat tmp;
        tmp.allocator = &g_numpyAllocator;
        transpose(ret, tmp);
        ret = tmp;
    }
    return ret;
}

cv::Mat PyObjectTocvMat(PyObject* o){

    if( !PyArray_Check(o) )
    {
        failmsg("toMat: Object is not a numpy array");
        return cv::Mat();
    }

    return PyArrayObjTocvMat((PyArrayObject*) o); 
}


/************ Auxiliary MSGranul Python functions definition ************/

void pyMaxLocalObjectFromMAXLOCAL(PyMaxLocalObject *maxlocalobj, MAXLOCAL local){
    maxlocalobj->correlation = local.r;
    maxlocalobj->shape      = local.forma;
    maxlocalobj->x          = local.ci;
    maxlocalobj->y          = local.li;
    maxlocalobj->w          = local.d1;
    maxlocalobj->h          = local.d2;
    maxlocalobj->ang        = local.ang;
}

void MAXLOCALFromPyMaxLocalObject(MAXLOCAL *local, PyMaxLocalObject maxlocalobj){
    local->r     = maxlocalobj.correlation;
    local->forma = maxlocalobj.shape      ;
    local->ci    = maxlocalobj.x          ;
    local->li    = maxlocalobj.y          ;
    local->d1    = maxlocalobj.w          ;
    local->d2    = maxlocalobj.h          ;
    local->ang   = maxlocalobj.ang        ;
} 

void pyKernelObjectFromKernel(PyKernelObject *pykernelobj, KERNEL kernel){

    pykernelobj->shape  = kernel.shape;
    pykernelobj->w      = kernel.d1;
    pykernelobj->h      = kernel.d2;
    pykernelobj->ang    = kernel.ang;

    pykernelobj->imgry   =   (PyObject *) cvMatToPyObject(kernel.imgG);
    pykernelobj->imfloat =   (PyObject *) cvMatToPyObject(kernel.imgF);
}

void kernelfromPyKernelObject(KERNEL* kernel, PyKernelObject pykernelobj){

    kernel->shape    =   pykernelobj.shape  ;
    kernel->d1       =   pykernelobj.w      ;
    kernel->d2       =   pykernelobj.h      ;
    kernel->ang      =   pykernelobj.ang    ;
    kernel->imgG     =   PyObjectTocvMat(pykernelobj.imgry);
    kernel->imgF     =   PyObjectTocvMat(pykernelobj.imfloat); 
}

PyObject* convertKernelsToPyList(vector<KERNEL> kernels){
    PyObject* retList = PyList_New(kernels.size());
    int list_index = 0;

    for (std::vector<KERNEL>::iterator it = kernels.begin(); it != kernels.end(); ++it){
        PyKernelObject *item = (PyKernelObject *) (PyKernelType).tp_alloc(&PyKernelType, 0);

        pyKernelObjectFromKernel(item, *it);

        if(PyList_SetItem(retList, list_index, (PyObject*) item) < 0)
            return NULL;
        ++list_index;
    }

    return retList;
}

vector<KERNEL> convertPyListToKernels(PyObject* list){
    vector<KERNEL> retvector = vector<KERNEL>(PyList_Size(list));

    for (int idx = 0; idx < retvector.size(); ++idx){
        PyKernelObject *item = (PyKernelObject *) PyList_GetItem(list, idx);

        retvector[idx] = KERNEL();
        kernelfromPyKernelObject(&retvector[idx], *item);
    }

    return retvector;
}

#define PyMaxLocalObject(X) ((PyMaxLocalObject*) X)

PyObject* convertMaxLocalsToPyList(vector<MAXLOCAL> maxlocals){
    PyObject* retList = PyList_New(maxlocals.size());
    int list_index = 0;

    for (std::vector<MAXLOCAL>::iterator it = maxlocals.begin(); it != maxlocals.end(); ++it){
        PyMaxLocalObject *item = PyMaxLocalObject((PyMaxLocalType).tp_alloc(&PyMaxLocalType, 0));

        pyMaxLocalObjectFromMAXLOCAL(item, *it);

        if(PyList_SetItem(retList, list_index, (PyObject*) item) < 0)
            return NULL;
        ++list_index;
    }

    return retList;
}

vector<MAXLOCAL> convertPyListToMaxLocals(PyObject* list){

    vector<MAXLOCAL> retvector = vector<MAXLOCAL>(PyList_Size(list));
    
    for (int idx = 0; idx < retvector.size(); ++idx){
        PyMaxLocalObject* item = (PyMaxLocalObject *) PyList_GetItem(list, idx);

        retvector[idx] = MAXLOCAL();
        MAXLOCALFromPyMaxLocalObject(&retvector[idx], *item);
    }

    return retvector;
}

void sortPyMAXLOCALList(PyListObject* list, int l, int r){
	int i=l; 
	int j=r;
	PyMaxLocalObject* x= PyMaxLocalObject(list->ob_item[(l+r)/2]);

	do {
		while (PyMaxLocalObject(list->ob_item[i])->correlation > x->correlation) i++;
		while (x->correlation  >  PyMaxLocalObject(list->ob_item[j])->correlation) j--;

		if (i<=j) {
			swap(list->ob_item[i],list->ob_item[j]); i++; j--;
		}
	} while (i<=j);

	if (l<j) sortPyMAXLOCALList(list,l,j);
	if (i<r) sortPyMAXLOCALList(list,i,r);
}

PyObject* sortPyMAXLOCAL(PyObject* list)
{

    if( !PyList_Check(list) )
    {
        failmsg("list: Object is not a Python list");
        return NULL;
    }

    int low = 0, high = PyList_Size(list) - 1;

    PyObject* retList = PyList_GetSlice(list, low, high);

    sortPyMAXLOCALList((PyListObject*) retList, low, high);

    return retList;
}

PyObject* correlate(PyObject* img, PyObject* kers, float minCorr, int maxDist, int type){

    if( !PyList_Check(kers) )
    {
        failmsg("list: Object is not a Python list");
        return NULL;
    }

    vector<KERNEL> kernels = convertPyListToKernels(kers);
    vector<MAXLOCAL> output, ret_locals;
    cv::Mat in = PyObjectTocvMat(img);

    if (in.channels() == 3) {
        cv::cvtColor(in, in, CV_RGB2GRAY);  
    } else if (in.channels() == 4) {
        cv::cvtColor(in, in, CV_RGBA2GRAY);  
    }

    correlationInBatch(in, output, kernels, minCorr, type);
    
    removeCloser(output, ret_locals, maxDist);

    return convertMaxLocalsToPyList(ret_locals);
}

PyObject* remove_maxlocal_closer(PyObject* maxLocalList, int maxDist){

    if( !PyList_Check(maxLocalList) )
    {
        failmsg("list: Object is not a Python list");
        return NULL;
    }

    vector<MAXLOCAL> ret = vector<MAXLOCAL>();

    vector<MAXLOCAL> maxlocals = convertPyListToMaxLocals(maxLocalList);

    removeCloser(maxlocals, ret, maxDist);

    return convertMaxLocalsToPyList(ret);
}

PyObject* apply_granul(PyObject* img, PyObject* maxLocalList, double minCorrDef, double maxInterDef){

    if( !PyList_Check(maxLocalList) )
    {
        failmsg("list: Object is not a Python list");
        return NULL;
    }

    vector<MAXLOCAL> ret = vector<MAXLOCAL>();

    vector<MAXLOCAL> maxlocals = convertPyListToMaxLocals(maxLocalList);

    cv::Mat input_img = PyObjectTocvMat(img);

    double minCorrCir, minCorrRet, minCorrQua, minCorrEli;
    minCorrCir = minCorrRet = minCorrQua = minCorrEli = minCorrDef;

    double maxInterCir, maxInterRet, maxInterQua, maxInterEli;
    maxInterCir = maxInterRet = maxInterQua = maxInterEli = maxInterDef;

    sift(maxlocals, ret, minCorrCir, maxInterCir, minCorrRet, maxInterRet, minCorrQua, maxInterQua, minCorrEli, maxInterEli, input_img.size());

    Mat output_img = input_img.clone();

    printMAXLOCAL(ret, output_img);

    PyObject* retList = PyList_New(2);

    if(PyList_SetItem(retList, 0, convertMaxLocalsToPyList(ret)) < 0)
            return NULL;

    if(PyList_SetItem(retList, 1, cvMatToPyObject( output_img)) < 0)
            return NULL;

    return retList;
}

PyObject* apply_msgranul(PyObject* img, PyObject* maxLocalList, double minCorrMSER){

    if( !PyList_Check(maxLocalList) )
    {
        failmsg("list: Object is not a Python list");
        return NULL;
    }

    vector<MAXLOCAL> ret = vector<MAXLOCAL>();

    vector<MAXLOCAL> maxlocals = convertPyListToMaxLocals(maxLocalList);

    cv::Mat input_img = PyObjectTocvMat(img);

    siftMSER(input_img, maxlocals, ret, minCorrMSER);

    Mat output_img = input_img.clone();

    printMAXLOCAL(ret, output_img);

    PyObject* retList = PyList_New(2);

    if(PyList_SetItem(retList, 0, convertMaxLocalsToPyList(ret)) < 0)
            return NULL;

    if(PyList_SetItem(retList, 1, cvMatToPyObject( output_img)) < 0)
            return NULL;

    return retList;
}

PyObject* extract_locals(PyObject* img, PyObject* maxLocalList){

    if( !PyList_Check(maxLocalList) )
    {
        failmsg("list: Object is not a Python list");
        return NULL;
    }

    vector<Mat> ret = vector<Mat>();

    vector<MAXLOCAL> maxlocals = convertPyListToMaxLocals(maxLocalList);

    cv::Mat input_img = PyObjectTocvMat(img);

    extractLocals(input_img, maxlocals, ret);

    PyObject* retList = PyList_New(ret.size());

    
    for(size_t idx = 0; idx < ret.size(); idx++)
    {
        if(PyList_SetItem(retList, idx, cvMatToPyObject(ret[idx])) < 0)
            return NULL;
    }

    return retList;
} 
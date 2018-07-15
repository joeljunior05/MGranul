#include "cv_np_granul.cpp"

/*
    This function will be called from Python to create Kernels [msgranul.createkernels]
*/
static PyObject *msgranul_createkernels(PyObject *self, PyObject *args){

    char shape;
    int len_rows, len_cols; 
    int octave = 5, steps = 10;
    float min = 0.5, max = 1;
    

    if (!PyArg_ParseTuple(args, "Cii|iffi", &shape, &len_rows, &len_cols, //Mandatory parameters (shape, length_of_rows, length_of_cols)
                                            &octave, &min, &max, &steps)) //Optional parameters (octave, min, max, angle_steps)
        return NULL;

    vector<KERNEL> kers;

    createKernels(kers, shape, len_rows, len_cols, octave, max, min, steps);

    return convertKernelsToPyList(kers);
}

/*
    This function will be called from Python to sort a list of MaxLocals [msgranul.sortlocal]
*/
static PyObject *msgranul_sortlocal(PyObject *self, PyObject *args){
    PyObject* maxLocalList;

    if (!PyArg_ParseTuple(args, "O", &maxLocalList)) //If no arg is provided it returns with error
        return NULL;

    return sortPyMAXLOCAL(maxLocalList);
}

/*
    This function will be called from Python to correlate an Image and a list of kernels [msgranul.correlate]
*/
static PyObject *msgranul_correlate(PyObject *self, PyObject *args){
    PyObject* img;
    PyObject* kers;
    float minCorr=0.0001;
    int maxDist=2, type = CV_TM_CCORR;

    if (!PyArg_ParseTuple(args, "OO|fii", &img, &kers, &minCorr, &maxDist, &type)) //If no arg is provided it returns with error
        return NULL;

    return correlate(img, kers, minCorr, maxDist, type);
}

/*
    This function will be called from Python to remove MaxLocals which are closer to each other [msgranul.removecloser]
*/
static PyObject *msgranul_removecloser(PyObject *self, PyObject *args){
    PyObject* maxLocalList;
    int maxDist;

    if (!PyArg_ParseTuple(args, "Oi", &maxLocalList, &maxDist)) //If no arg is provided it returns with error
        return NULL;

    return remove_maxlocal_closer(maxLocalList, maxDist);
}

/*
    This function will be called from Python to apply granulometry based on correlation [msgranul.apply]
*/
static PyObject *msgranul_apply(PyObject *self, PyObject *args){

    PyObject* img;
    PyObject* maxLocalList;

    double minCorrDef=0.1;

    double maxInterDef=0.3;

    if (!PyArg_ParseTuple(args, "OO|dd", &img, &maxLocalList, &minCorrDef, &maxInterDef)) //If no arg is provided it returns with error
        return NULL;


    return apply_granul(img, maxLocalList, minCorrDef, maxInterDef);
}

/*
    This function will be called from Python to apply MSGranul [msgranul.applyMSER]
*/
static PyObject *msgranul_applyMSER(PyObject *self, PyObject *args){
    PyObject* img;
    PyObject* maxLocalList;

    double minCorrMSER = 0.9;

    if (!PyArg_ParseTuple(args, "OO|d", &img, &maxLocalList, &minCorrMSER)) //If no arg is provided it returns with error
        return NULL;

    return apply_msgranul(img, maxLocalList, minCorrMSER);
}

static PyObject *msgranul_extract(PyObject *self, PyObject *args){
    PyObject* img;
    PyObject* maxLocalList;

    if (!PyArg_ParseTuple(args, "OO", &img, &maxLocalList)) //If no arg is provided it returns with error
        return NULL;

    return extract_locals(img, maxLocalList);
}

static PyObject *msgranul_printmaxlocals(PyObject *self, PyObject *args){
    PyObject* img;
    PyObject* maxLocalList;

    if (!PyArg_ParseTuple(args, "OO", &img, &maxLocalList)) //If no arg is provided it returns with error
        return NULL;

    return print_locals(img, maxLocalList);
}

static PyMethodDef MSGranulMethods[] = {
    {"createkernels",   msgranul_createkernels, METH_VARARGS, "Create a list of kernels."},
    {"sortlocal",       msgranul_sortlocal, METH_VARARGS, "Sort a list of locals by correlation value."},
    {"correlate",       msgranul_correlate, METH_VARARGS, "Correlate a list of kernels and an image."},
    {"removecloser",    msgranul_removecloser, METH_VARARGS, "Remove from list the maxlocals closer."},
    {"apply",           msgranul_apply, METH_VARARGS, "Apply Granulometry based on correlation using kernels."},
    {"applyMSER",       msgranul_applyMSER, METH_VARARGS, "Apply MSGranul using kernels."},
    {"extract",         msgranul_extract, METH_VARARGS, "Extract a piece of image where MaxLocal is located."},
    {"printmaxlocals",  msgranul_printmaxlocals, METH_VARARGS, "Prints MaxLocals on the image."},

    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef msgranulmodule = {
    PyModuleDef_HEAD_INIT,
    "msgranul",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    MSGranulMethods
};

int checkIfAllTypeAreReady(){
    return PyType_Ready(&PyMaxLocalType) + 
           PyType_Ready(&PyKernelType);
}

PyMODINIT_FUNC
PyInit_msgranul(void)
{
    PyObject *m;

    if(checkIfAllTypeAreReady() < 0)
        return NULL;

    m = PyModule_Create(&msgranulmodule);

    if (m == NULL)
        return NULL;

    Py_INCREF(&PyMaxLocalType);
    PyModule_AddObject(m, "MaxLocal", (PyObject *) &PyMaxLocalType);
    init_numpy();
    return m;
}
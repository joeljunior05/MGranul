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

    /*KERNEL ker = kers[0];

    printf("%f\n", ker.d1);

    PyKernelObject *local = (PyKernelObject *) (PyKernelType).tp_alloc(&PyKernelType, 0);

    pyKernelObjectFromKernel(local, ker);

    Py_INCREF(local);
    return (PyObject *) local;*/
    //Py_RETURN_NONE;
}

/*
    This function will be called from Python to sort a list of MaxLocals [msgranul.sortlocal]
*/
static PyObject *msgranul_sortlocal(PyObject *self, PyObject *args){
    const char *command;

    if (!PyArg_ParseTuple(args, "s", &command)) //If no arg is provided it return with error
        return NULL;

    Py_RETURN_NONE;
}

/*
    This function will be called from Python to correlate an Image and a list of kernels [msgranul.correlate]
*/
static PyObject *msgranul_correlate(PyObject *self, PyObject *args){
    const char *command;

    if (!PyArg_ParseTuple(args, "s", &command)) //If no arg is provided it return with error
        return NULL;

    Py_RETURN_NONE;
}

/*
    This function will be called from Python to remove MaxLocals which are closer to each other [msgranul.removecloser]
*/
static PyObject *msgranul_removecloser(PyObject *self, PyObject *args){
    const char *command;

    if (!PyArg_ParseTuple(args, "s", &command)) //If no arg is provided it return with error
        return NULL;

    Py_RETURN_NONE;
}

/*
    This function will be called from Python to apply granulometry based on correlation [msgranul.apply]
*/
static PyObject *msgranul_apply(PyObject *self, PyObject *args){
    const char *command;

    if (!PyArg_ParseTuple(args, "s", &command)) //If no arg is provided it return with error
        return NULL;

    Py_RETURN_NONE;
}

/*
    This function will be called from Python to apply MSGranul [msgranul.applyMSER]
*/
static PyObject *msgranul_applyMSER(PyObject *self, PyObject *args){
    const char *command;

    if (!PyArg_ParseTuple(args, "s", &command)) //If no arg is provided it return with error
        return NULL;

    Py_RETURN_NONE;
}

static PyMethodDef MSGranulMethods[] = {
    {"createkernels",   msgranul_createkernels, METH_VARARGS, "Create a list of kernels."},
    {"sortlocal",       msgranul_createkernels, METH_VARARGS, "Sort a list of locals by correlation value."},
    {"correlate",       msgranul_createkernels, METH_VARARGS, "Correlate a list of kernels and an image."},
    {"removecloser",    msgranul_createkernels, METH_VARARGS, "Remove from list the maxlocals closer."},
    {"apply",           msgranul_createkernels, METH_VARARGS, "Apply Granulometry based on correlation using kernels."},
    {"applyMSER",       msgranul_createkernels, METH_VARARGS, "Apply MSGranul using kernels."},
    
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
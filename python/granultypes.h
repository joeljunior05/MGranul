#include <Python.h>
#include <structmember.h>

/************ Types definition ************/

/* Begin - MaxLocal definition */

typedef struct PyMaxLocalObject{ // Object represent the object itself.
    PyObject_HEAD
    float correlation;
    char shape;
    int x;
    int y;
    int w;
    int h;
    float ang;
} PyMaxLocalObject;

static PyMemberDef maxlocal_members[] = {
    {"correlation", T_FLOAT, offsetof(PyMaxLocalObject, correlation), 0,
     "correlation value"},
    {"shape", T_CHAR, offsetof(PyMaxLocalObject, shape), 0,
     "kind of shape"},
    {"x", T_INT, offsetof(PyMaxLocalObject, x), 0,
     "center -> x"},
    {"y", T_INT, offsetof(PyMaxLocalObject, y), 0,
     "center -> y"},
    {"w", T_INT, offsetof(PyMaxLocalObject, w), 0,
     "width"},
    {"h", T_INT, offsetof(PyMaxLocalObject, h), 0,
     "heigh"},
    {"ang", T_INT, offsetof(PyMaxLocalObject, ang), 0,
     "angle in degree"},
    {NULL}  /* Sentinel */
};

static PyTypeObject PyMaxLocalType = { // Specify how an object will behave.
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name    = "msgranul.maxlocal",
    .tp_doc     = "MaxLocal objects",
    .tp_basicsize = sizeof(PyMaxLocalObject),
    .tp_itemsize = 0,
    .tp_flags   = Py_TPFLAGS_DEFAULT,
    .tp_new     = PyType_GenericNew,
    .tp_members = maxlocal_members
};

/* End - MaxLocal definition */

/* Begin - Kernel definition */

typedef struct PyKernelObject{
    PyObject_HEAD
    PyObject *imgry;
    PyObject *imfloat;
    char shape;
    int w;
    int h;
    float ang;
} PyKernelObject;

static PyMemberDef kernel_members[] = {
    {"imgry", T_OBJECT_EX, offsetof(PyKernelObject, imgry), 0,
     "grayscale image"},
    {"imfloat", T_OBJECT_EX, offsetof(PyKernelObject, imfloat), 0,
     "float image"},
    {"shape", T_CHAR, offsetof(PyKernelObject, shape), 0,
     "kind of shape"},
    {"w", T_INT, offsetof(PyKernelObject, w), 0,
     "width"},
    {"h", T_INT, offsetof(PyKernelObject, h), 0,
     "heigh"},
    {"ang", T_INT, offsetof(PyKernelObject, ang), 0,
     "angle in degree"},
    {NULL}  /* Sentinel */
};

static PyTypeObject PyKernelType = { // Specify how an object will behave.
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name    = "msgranul.kernel",
    .tp_doc     = "Kernel objects",
    .tp_basicsize = sizeof(PyKernelObject),
    .tp_itemsize = 0,
    .tp_flags   = Py_TPFLAGS_DEFAULT,
    .tp_new     = PyType_GenericNew,
    .tp_members = kernel_members
};

/* End - Kernel definition */
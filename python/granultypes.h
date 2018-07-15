#include <Python.h>
#include <structmember.h>

/************ Types definition ************/

#define safe_cvt(string_constant) ((char*) string_constant)

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
    {safe_cvt("correlation"), T_FLOAT, offsetof(PyMaxLocalObject, correlation), 0,
     safe_cvt("correlation value")},
    {safe_cvt("shape"), T_CHAR, offsetof(PyMaxLocalObject, shape), 0,
     safe_cvt("kind of shape")},
    {safe_cvt("x"), T_INT, offsetof(PyMaxLocalObject, x), 0,
     safe_cvt("center -> x")},
    {safe_cvt("y"), T_INT, offsetof(PyMaxLocalObject, y), 0,
     safe_cvt("center -> y")},
    {safe_cvt("w"), T_INT, offsetof(PyMaxLocalObject, w), 0,
     safe_cvt("width")},
    {safe_cvt("h"), T_INT, offsetof(PyMaxLocalObject, h), 0,
     safe_cvt("heigh")},
    {safe_cvt("ang"), T_INT, offsetof(PyMaxLocalObject, ang), 0,
     safe_cvt("angle in degree")},
    {NULL}  /* Sentinel */
};

static PyTypeObject PyMaxLocalType = { // Specify how an object will behave.
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name    = "msgranul.maxlocal",
    .tp_basicsize = sizeof(PyMaxLocalObject),
    .tp_itemsize = 0,
    .tp_flags   = Py_TPFLAGS_DEFAULT,
    .tp_doc     = "MaxLocal objects",
    .tp_members = maxlocal_members,
    .tp_new     = PyType_GenericNew
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
    {safe_cvt("imgry"), T_OBJECT_EX, offsetof(PyKernelObject, imgry), 0,
     safe_cvt("grayscale image")},
    {safe_cvt("imfloat"), T_OBJECT_EX, offsetof(PyKernelObject, imfloat), 0,
     safe_cvt("float image")},
    {safe_cvt("shape"), T_CHAR, offsetof(PyKernelObject, shape), 0,
     safe_cvt("kind of shape")},
    {safe_cvt("w"), T_INT, offsetof(PyKernelObject, w), 0,
     safe_cvt("width")},
    {safe_cvt("h"), T_INT, offsetof(PyKernelObject, h), 0,
     safe_cvt("heigh")},
    {safe_cvt("ang"), T_INT, offsetof(PyKernelObject, ang), 0,
     safe_cvt("angle in degree")},
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
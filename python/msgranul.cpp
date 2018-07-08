#include <Python.h>

typedef struct PyMaxLocalObject{
    PyObject_HEAD
    float correlation;
    char shape;
    int x;
    int y;
    int w;
    int h;
    float ang;
} PyMaxLocalObject;

typedef struct PyKernelObject{
    PyObject_HEAD
    char shape;
    int w;
    int h;
    float ang;
} PyKernelObject;
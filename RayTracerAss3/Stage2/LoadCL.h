#ifndef _CLLOADSOURCE_H
#define _CLLOADSOURCE_H

#include <CL/cl.h>

cl_program clLoadSource(cl_context context, char* filename, cl_int* err);

#endif
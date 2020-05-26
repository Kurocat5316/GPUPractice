#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdio.h>
#include <stdlib.h>
#include "LoadCL.h"

cl_program clLoadSource(cl_context context, char* filename, cl_int* err)
{
	cl_program program;
	FILE *program_handle;
	char *program_buffer;
	size_t program_size;

	program_handle = fopen(filename, "rb");
	if (program_handle == NULL) {
		printf("Couldn't find the program file\n");
		exit(1);
	}
	fseek(program_handle, 0, SEEK_END);
	program_size = ftell(program_handle);
	rewind(program_handle);
	program_buffer = (char*)malloc(program_size + 1);
	program_buffer[program_size] = '\0';
	fread(program_buffer, sizeof(char), program_size, program_handle);
	fclose(program_handle);

	program = clCreateProgramWithSource(context, 1, (const char**)&program_buffer, &program_size, err);
	free(program_buffer);

	return program;
}

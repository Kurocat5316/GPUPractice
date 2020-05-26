/*  The following code is a VERY heavily modified from code originally sourced from:
Ray tracing tutorial of http://www.codermind.com/articles/Raytracer-in-C++-Introduction-What-is-ray-tracing.html
It is free to use for educational purpose and cannot be redistributed outside of the tutorial pages. */

#define TARGET_WINDOWS

#pragma warning(disable: 4996)
#include "Timer.h"
#include "Primitives.h"
#include "Scene.h"
#include "Lighting.h"
#include "Intersection.h"
#include "ImageIO.h"
#include <CL/cl.h>

#include "LoadCL.h"

unsigned int buffer[MAX_WIDTH * MAX_HEIGHT];

// reflect the ray from an object
Ray calculateReflection(const Ray* viewRay, const Intersection* intersect)
{
	// reflect the viewRay around the object's normal
	Ray newRay = { intersect->pos, viewRay->dir - (intersect->normal * intersect->viewProjection * 2.0f) };

	return newRay;
}


// refract the ray through an object
Ray calculateRefraction(const Ray* viewRay, const Intersection* intersect, float* currentRefractiveIndex)
{
	// change refractive index depending on whether we are in an object or not
	float oldRefractiveIndex = *currentRefractiveIndex;
	*currentRefractiveIndex = intersect->insideObject ? DEFAULT_REFRACTIVE_INDEX : intersect->material->density;

	// calculate refractive ratio from old index and current index
	float refractiveRatio = oldRefractiveIndex / *currentRefractiveIndex;

	// Here we take into account that the light movement is symmetrical from the observer to the source or from the source to the oberver.
	// We then do the computation of the coefficient by taking into account the ray coming from the viewing point.
	float fCosThetaT;
	float fCosThetaI = fabsf(intersect->viewProjection);

	// glass-like material, we're computing the fresnel coefficient.
	if (fCosThetaI >= 1.0f)
	{
		// In this case the ray is coming parallel to the normal to the surface
		fCosThetaT = 1.0f;
	}
	else
	{
		float fSinThetaT = refractiveRatio * sqrtf(1 - fCosThetaI * fCosThetaI);

		// Beyond the angle (1.0f) all surfaces are purely reflective
		fCosThetaT = (fSinThetaT * fSinThetaT >= 1.0f) ? 0.0f : sqrtf(1 - fSinThetaT * fSinThetaT);
	}

	// Here we compute the transmitted ray with the formula of Snell-Descartes
	Ray newRay = { intersect->pos, (viewRay->dir + intersect->normal * fCosThetaI) * refractiveRatio - (intersect->normal * fCosThetaT) };

	return newRay;
}


// follow a single ray until it's final destination (or maximum number of steps reached)
Colour traceRay(const Scene* scene, Ray viewRay)
{
	Colour output(0.0f, 0.0f, 0.0f); 								// colour value to be output
	float currentRefractiveIndex = DEFAULT_REFRACTIVE_INDEX;		// current refractive index
	float coef = 1.0f;												// amount of ray left to transmit
	Intersection intersect;											// properties of current intersection

																	// loop until reached maximum ray cast limit (unless loop is broken out of)
	for (int level = 0; level < MAX_RAYS_CAST; ++level)
	{
		// check for intersections between the view ray and any of the objects in the scene
		// exit the loop if no intersection found
		if (!objectIntersection(scene, &viewRay, &intersect)) break;

		// calculate response to collision: ie. get normal at point of collision and material of object
		calculateIntersectionResponse(scene, &viewRay, &intersect);

		// apply the diffuse and specular lighting 
		if (!intersect.insideObject) output += coef * applyLighting(scene, &viewRay, &intersect);

		// if object has reflection or refraction component, adjust the view ray and coefficent of calculation and continue looping
		//if (intersect.material->reflection)
		//{
		//	viewRay = calculateReflection(&viewRay, &intersect);
		//	coef *= intersect.material->reflection;
		//}
		//else if (intersect.material->refraction)
		//{
		//	viewRay = calculateRefraction(&viewRay, &intersect, &currentRefractiveIndex);
		//	coef *= intersect.material->refraction;
		//}
		//else
		//{
			// if no reflection or refraction, then finish looping (cast no more rays)
			return output;
		//}
	}

	// if the calculation coefficient is non-zero, read from the environment map
	if (coef > 0.0f)
	{
		Material& currentMaterial = scene->materialContainer[scene->skyboxMaterialId];

		output += coef * currentMaterial.diffuse;
	}

	return output;
}


// render scene at given width and height and anti-aliasing level
void render(Scene* scene, const int width, const int height, const int aaLevel)
{
	// angle between each successive ray cast (per pixel, anti-aliasing uses a fraction of this)
	const float dirStepSize = 1.0f / (0.5f * width / tanf(PIOVER180 * 0.5f * scene->cameraFieldOfView));

	// pointer to output buffer
	unsigned int* out = buffer;

	// loop through all the pixels
	for (int y = -height / 2; y < height / 2; ++y)
	{
		for (int x = -width / 2; x < width / 2; ++x)
		{
			Colour output(0.0f, 0.0f, 0.0f);

			// calculate multiple samples for each pixel
			const float sampleStep = 1.0f / aaLevel, sampleRatio = 1.0f / (aaLevel * aaLevel);

			// loop through all sub-locations within the pixel
			for (float fragmentx = float(x); fragmentx < x + 1.0f; fragmentx += sampleStep)
			{
				for (float fragmenty = float(y); fragmenty < y + 1.0f; fragmenty += sampleStep)
				{
					// direction of default forward facing ray
					Vector dir = { fragmentx * dirStepSize, fragmenty * dirStepSize, 1.0f };

					// rotated direction of ray
					Vector rotatedDir = {
						dir.x * cosf(scene->cameraRotation) - dir.z * sinf(scene->cameraRotation),
						dir.y,
						dir.x * sinf(scene->cameraRotation) + dir.z * cosf(scene->cameraRotation) };

					// view ray starting from camera position and heading in rotated (normalised) direction
					Ray viewRay = { scene->cameraPosition, normalise(rotatedDir) };

					// follow ray and add proportional of the result to the final pixel colour
					output += sampleRatio * traceRay(scene, viewRay);
				}
			}

			// store saturated final colour value in image buffer
			*out++ = output.convertToPixel(scene->exposure);
		}
	}
}

void renderGPU(Scene scene, const int width, const int height, const int aaLevel)
{
	cl_int err;
	cl_platform_id platform;
	cl_device_id device;
	cl_context context;
	cl_command_queue queue;
	cl_program program;
	cl_kernel kernel;
	cl_mem clBufferScene;
	cl_mem clBufferCameraPosition;
	cl_mem clBufferMaterialContainer;
	cl_mem clBufferSphereContainer;
	cl_mem clBufferTriangleContainer;
	cl_mem clBufferLightContainer;
	cl_mem clBufferAaLevel;
	cl_mem clBufferOut;
	size_t workOffset[] = { 0, 0 };
	size_t workSize[] = { width, height };

	err = clGetPlatformIDs(1, &platform, NULL);
	if (err != CL_SUCCESS)
	{
		printf("\nError calling clGetPlatformIDs. Error code: %d\n", err);
		exit(1);
	}

	err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, NULL);
	if (err != CL_SUCCESS) {
		printf("Couldn't find any devices\n");
		exit(1);
	}

	context = clCreateContext(NULL, 1, &device, NULL, NULL, &err);
	if (err != CL_SUCCESS) {
		printf("Couldn't create a context\n");
		exit(1);
	}

	queue = clCreateCommandQueue(context, device, 0, &err);
	if (err != CL_SUCCESS) {
		printf("Couldn't create the command queue\n");
		exit(1);
	}

	char fileName[10];
	strcpy(fileName, "stage3.cl");

	program = clLoadSource(context, fileName, &err);
	if (err != CL_SUCCESS) {
		printf("Couldn't load/create the program\n");
		exit(1);
	}

	err = clBuildProgram(program, 0, NULL, NULL /*"-save-temps=C:/TEMP/OCL"*/, NULL, NULL);
	if (err != CL_SUCCESS) {
		char *program_log;
		size_t log_size;

		clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
		program_log = (char*)malloc(log_size + 1);
		program_log[log_size] = '\0';
		clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, log_size + 1, program_log, NULL);
		printf("%s\n", program_log);
		free(program_log);
		exit(1);
	}

	kernel = clCreateKernel(program, "GetDateFromScene", &err);
	if (err != CL_SUCCESS) {
		printf("Couldn't create the kernel\n");
		exit(1);
	}

	clBufferScene = clCreateBuffer(context, CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(Scene), &scene, &err);
	if (err != CL_SUCCESS) {
		printf("Couldn't create a clBufferScene object\n");
		exit(1);
	}

	clBufferCameraPosition = clCreateBuffer(context, CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(Point), &scene.cameraPosition, &err);
	if (err != CL_SUCCESS) {
		printf("Couldn't create a clBufferScene object\n");
		exit(1);
	}

	clBufferMaterialContainer = clCreateBuffer(context, CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(Material) * scene.numMaterials, scene.materialContainer, &err);
	if (err != CL_SUCCESS) {
		printf("Couldn't create a clBufferScene object\n");
		exit(1);
	}

	clBufferSphereContainer = clCreateBuffer(context, CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(Sphere) * scene.numSpheres, scene.sphereContainer, &err);
	if (err != CL_SUCCESS) {
		printf("Couldn't create a clBufferScene object\n");
		exit(1);
	}

	clBufferTriangleContainer = clCreateBuffer(context, CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(Triangle) * scene.numTriangles, scene.triangleContainer, &err);
	if (err != CL_SUCCESS) {
		printf("Couldn't create a clBufferScene object\n");
		exit(1);
	}

	clBufferLightContainer = clCreateBuffer(context, CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(Light) * scene.numLights, scene.lightContainer, &err);
	if (err != CL_SUCCESS) {
		printf("Couldn't create a clBufferScene object\n");
		exit(1);
	}

	clBufferAaLevel = clCreateBuffer(context, CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(int), (void*)&aaLevel, &err);
	if (err != CL_SUCCESS) {
		printf("Couldn't create a clBufferAaLevel object\n");
		exit(1);
	}

	clBufferOut = clCreateBuffer(context, CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(unsigned int) * width * height, buffer, &err);
	if (err != CL_SUCCESS) {
		printf("Couldn't create a clBufferScene object\n");
		exit(1);
	}


	// pass in all the kernel arguments 
	// (one argument per basic scalar -- alternatively these could be wrapped in a struct and passed in as one argument)
	err = clSetKernelArg(kernel, 0, sizeof(cl_mem), &clBufferScene);
	if (err != CL_SUCCESS) {
		printf("Couldn't set the kernel 0 argument\n");
		exit(1);
	}

	err = clSetKernelArg(kernel, 1, sizeof(cl_mem), &clBufferCameraPosition);
	if (err != CL_SUCCESS) {
		printf("Couldn't set the kernel 0 argument\n");
		exit(1);
	}

	err = clSetKernelArg(kernel, 2, sizeof(cl_mem), &clBufferMaterialContainer);
	if (err != CL_SUCCESS) {
		printf("Couldn't set the kernel 1 argument\n");
		exit(1);
	}

	err = clSetKernelArg(kernel, 3, sizeof(cl_mem), &clBufferSphereContainer);
	if (err != CL_SUCCESS) {
		printf("Couldn't set the kernel 2 argument\n");
		exit(1);
	}

	err = clSetKernelArg(kernel, 4, sizeof(cl_mem), &clBufferTriangleContainer);
	if (err != CL_SUCCESS) {
		printf("Couldn't set the kernel 3 argument\n");
		exit(1);
	}

	err = clSetKernelArg(kernel, 5, sizeof(cl_mem), &clBufferLightContainer);
	if (err != CL_SUCCESS) {
		printf("Couldn't set the kernel 4 argument\n");
		exit(1);
	}

	err = clSetKernelArg(kernel, 6, sizeof(cl_mem), &clBufferAaLevel);
	if (err != CL_SUCCESS) {
		printf("Couldn't set the kernel 5 argument\n");
		exit(1);
	}

	err = clSetKernelArg(kernel, 7, sizeof(cl_mem), &clBufferOut);
	if (err != CL_SUCCESS) {
		printf("Couldn't set the kernel 6 argument\n");
		exit(1);
	}


	err = clEnqueueNDRangeKernel(queue, kernel, 2, workOffset, workSize, NULL, 0, NULL, NULL);
	if (err != CL_SUCCESS) {
		printf("Couldn't enqueue the kernel execution command\n");
		exit(1);
	}

	err = clEnqueueReadBuffer(queue, clBufferOut, CL_TRUE, 0, sizeof(int) * width * height, buffer, 0, NULL, NULL);
	if (err != CL_SUCCESS) {
		printf("Couldn't enqueue the read buffer command\n");
		exit(1);
	}



	// clean up
	clReleaseMemObject(clBufferOut);
	clReleaseMemObject(clBufferAaLevel);
	clReleaseMemObject(clBufferScene);
	clReleaseMemObject(clBufferCameraPosition);
	clReleaseMemObject(clBufferMaterialContainer);
	clReleaseMemObject(clBufferSphereContainer);
	clReleaseMemObject(clBufferTriangleContainer);
	clReleaseMemObject(clBufferLightContainer);
	clReleaseCommandQueue(queue);
	clReleaseProgram(program);
	clReleaseKernel(kernel);
	clReleaseContext(context);

}


// read command line arguments, render, and write out BMP file
int main(int argc, char* argv[])
{
	int width = 500;
	int height = 300;
	int samples = 1;

	// rendering options
	int times = 1;
	unsigned int threads = 1;			// currently unused
	bool colourise = false;				// currently unused
	unsigned int blockSize = -1;		// currently unused

	// default input / output filenames
	const char* inputFilename = "../Scenes/cornell.txt";
	char outputFilenameBuffer[1000];
	char* outputFilename = outputFilenameBuffer;

	// do stuff with command line args
	for (int i = 1; i < argc; i++)
	{
		if (strcmp(argv[i], "-size") == 0)
		{
			width = atoi(argv[++i]);
			height = atoi(argv[++i]);
		}
		else if (strcmp(argv[i], "-samples") == 0)
		{
			samples = atoi(argv[++i]);
		}
		else if (strcmp(argv[i], "-input") == 0)
		{
			inputFilename = argv[++i];
		}
		else if (strcmp(argv[i], "-output") == 0)
		{
			outputFilename = argv[++i];
		}
		else if (strcmp(argv[i], "-runs") == 0)
		{
			times = atoi(argv[++i]);
		}
		else if (strcmp(argv[i], "-threads") == 0)
		{
			threads = atoi(argv[++i]);
		}
		else if (strcmp(argv[i], "-colourise") == 0)
		{
			colourise = true;
		}
		else if (strcmp(argv[i], "-blockSize") == 0)
		{
			blockSize = atoi(argv[++i]);
		}
		else
		{
			fprintf(stderr, "unknown argument: %s\n", argv[i]);
		}
	}

	// nasty (and fragile) kludge to make an ok-ish default output filename (can be overriden with "-output" command line option)
	sprintf(outputFilenameBuffer, "../Outputs/%s_%dx%dx%d_%s.bmp", (strrchr(inputFilename, '/') + 1), width, height, samples, (strrchr(argv[0], '\\') + 1));

	// read scene file
	Scene scene;
	if (!init(inputFilename, scene))
	{
		fprintf(stderr, "Failure when reading the Scene file.\n");
		return -1;
	}

	// total time taken to render all runs (used to calculate average)
	int totalTime = 0;
	for (int i = 0; i < times; i++)
	{
		Timer timer;									// create timer
		//render(&scene, width, height, samples);			// raytrace scene
		renderGPU(scene, width, height, samples);			// raytrace scene
		timer.end();									// record end time
		totalTime += timer.getMilliseconds();			// record total time taken
	}

	// output timing information (times run and average)
	printf("average time taken (%d run(s)): %ums\n", times, totalTime / times);

	// output BMP file
	write_bmp(outputFilename, buffer, width, height, width);
}

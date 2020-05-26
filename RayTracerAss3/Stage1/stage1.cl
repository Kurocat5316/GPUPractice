typedef struct Point
{
	float x, y, z;
}Point;


typedef struct Colour
{
	float red, green, blue;
}Colour;

typedef struct Vector
{
	float x, y, z;
} Vector;

typedef struct Ray
{
	__kernel Point start;
	__kernel Vector dir;
} Ray;

// material
typedef struct Material
{
	// type of colouring/texturing
	enum { GOURAUD, CHECKERBOARD, CIRCLES, WOOD } type;

	__kernel Colour diffuse;				// diffuse colour
	__kernel Colour diffuse2;			// second diffuse colour, only for checkerboard types

	__kernel Vector offset;				// offset of generated texture
	float size;					// size of generated texture

	__kernel Colour specular;			// colour of specular lighting
	float power;				// power of specular reflection

	float reflection;			// reflection amount
	float refraction;			// refraction amount
	float density;				// density of material (affects amount of defraction)
} Material;


// sphere object
typedef struct Sphere
{
	__kernel Point pos;					// a point on the plane
	float size;					// radius of sphere
	unsigned int materialId;	// material id
} Sphere;


// light object
typedef struct Light
{
	__kernel Point pos;					// location
	__kernel Colour intensity;			// brightness and colour
} Light;


// triangle object
typedef struct Triangle
{
	__kernel Point p1, p2, p3;			// the three points of the triangle
	__kernel Vector normal;				// normal of the triangle
	unsigned int materialId;	// material id
} Triangle;

typedef struct Scene{
	__kernel Point cameraPosition;					// camera location
	float cameraRotation;					// direction camera points
	float cameraFieldOfView;				// field of view for the camera

	float exposure;							// image exposure

	unsigned int skyboxMaterialId;

	// scene object counts
	unsigned int numMaterials;
	unsigned int numSpheres;
	unsigned int numTriangles;
	unsigned int numLights;

	// scene objects
	__global Material* materialContainer;
	__global Sphere* sphereContainer;
	__global Triangle* triangleContainer;
	__global Light* lightContainer;
}Scene;

// convert colour to pixel (in 0x00BBGGRR format) with respect to an exposure level 
inline unsigned int convertToPixel(Colour colour, float exposure){

	unsigned int result = ((unsigned char) (255 * (min(1.0f - exp(colour.blue * exposure), 1.0f))) << 16) + 
		((unsigned char) (255 * (min(1.0f - exp(colour.green * exposure), 1.0f))) << 8) +
		((unsigned char) (255 * (min(1.0f - exp(colour.red * exposure), 1.0f))) << 0);
	return result;
}

__kernel void GetDateFromScene(__global struct Scene* scene,
	__global struct Point* cameraPosition,
	__global struct Material* materialContainer,
	__global struct Sphere* sphereContainer,
	__global struct Triangle* triangleContainer,
	__global struct Light* lightContainer,
	__global unsigned int* out)
{
	unsigned int width = get_global_size(0);
	unsigned int ix = get_global_id(0);
	unsigned int iy = get_global_id(1);

	scene->materialContainer = materialContainer;
	scene->sphereContainer = sphereContainer;
	scene->triangleContainer = triangleContainer;
	scene->lightContainer = lightContainer;
	
	if (ix == 0 && iy == 0) {
		printf("CameraPosition - x: %f \n", scene->cameraPosition.x);
		printf("CameraPosition - y: %f \n", scene->cameraPosition.y);
		printf("CameraPosition - z: %f \n", scene->cameraPosition.z);


		printf("CameraRotation: %f \n", scene->cameraRotation);
		printf("CameraFieldOfView: %f \n", scene->cameraFieldOfView);
		printf("Exposure: %f \n", scene->exposure);
		printf("SkyboxMaterialId: %d \n", scene->skyboxMaterialId);
		printf("NumMaterials: %d \n", scene->numMaterials);
		printf("NumSpheres: %d \n", scene->numSpheres);
		printf("NumTriangles: %d \n", scene->numTriangles);
		printf("NumLights: %d \n", scene->numLights);
		printf("\n");




		for (int i = 0; i < scene->numMaterials; i++) {
			printf("NumMaterialsMaterialContainer[%d]: %d\n", i, scene->materialContainer[i].size);
		}

		printf("\n");


		for (unsigned int i = 0; i < scene->numSpheres; i++) {
			printf("sphereContainerMaterialId[%d]: %d \n", i, scene->sphereContainer[i].materialId);
			printf("sphereContainerPos[%d]: (%f, %f, %f) \n", i, scene->sphereContainer[i].pos.x, scene->sphereContainer[i].pos.y, scene->sphereContainer[i].pos.z);
			printf("sphereContainerSize[%d]: %f\n", i, scene->sphereContainer[i].size);


		}



		printf("\n");





		for (int i = 0; i < scene->numTriangles; i++) {
			printf("triangleContainerP1[%d]: (%f, %f, %f) \n", i, scene->triangleContainer[i].p1.x, scene->triangleContainer[i].p1.y, scene->triangleContainer[i].p1.z);
			printf("triangleContainerP2[%d]: (%f, %f, %f) \n", i, scene->triangleContainer[i].p2.x, scene->triangleContainer[i].p2.y, scene->triangleContainer[i].p2.z);
			printf("triangleContainerP3[%d]: (%f, %f, %f) \n", i, scene->triangleContainer[i].p3.x, scene->triangleContainer[i].p3.y, scene->triangleContainer[i].p3.z);

			printf("triangleContainerNormal[%d]: (%f, %f, %f) \n", i, scene->triangleContainer[i].normal.x, scene->triangleContainer[i].normal.y, scene->triangleContainer[i].normal.z);
			printf("triangleContainerMaterialId[%d]: %f \n", i, scene->triangleContainer[i].materialId);
		}

		printf("\n");


		for (int i = 0; i < scene->numLights; i++) {
			printf("lightContainerPos[%d]: (%f, %f, %f) \n", i, scene->lightContainer[i].pos.x, scene->lightContainer[i].pos.y, scene->lightContainer[i].pos.z);
			printf("lightContainerIntensity[%d]: (%f, %f, %f) \n", i, scene->lightContainer[i].intensity.red, scene->lightContainer[i].intensity.green, scene->lightContainer[i].intensity.blue);
		}
	}

	//link output
	//out->cameraPosition = scene->cameraPosition;
	//out->cameraRotation = scene->cameraRotation;
	//out->cameraFieldOfView = scene->cameraFieldOfView;
	//out->exposure = scene->exposure;
	//out->skyboxMaterialId = scene->skyboxMaterialId;
	//out->numMaterials = scene->numMaterials;
	//out->numSpheres = scene->numSpheres;
	//out->numTriangles = scene->numTriangles;
	//out->numLights = scene->numLights;
	//
	//out->materialContainer = scene->materialContainer;
	//out->sphereContainer = scene->sphereContainer;
	//out->triangleContainer = scene->triangleContainer;
	//out->lightContainer = scene->lightContainer;
	Colour output;
	output.red = 255.0f;
	output.blue = 255.0f;
	output.green = 255.0f;

	out[width * iy + ix] = convertToPixel(output, scene->exposure);
}
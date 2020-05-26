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

typedef struct Scene {
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
inline unsigned int convertToPixel(Colour colour, float exposure) {
	unsigned int result = ((unsigned char) (255 * (min(1.0f - exp(colour.blue * exposure), 1.0f))) << 16) +
		((unsigned char) (255 * (min(1.0f - exp(colour.green * exposure), 1.0f))) << 8) +
		((unsigned char) (255 * (min(1.0f - exp(colour.red * exposure), 1.0f))) << 0);
	return result;
}

inline float invsqrtf(const float x)
{
	return 1.0f / sqrt(x);
}

inline Vector normalise(const Vector v)
{
	float dot = v.x * v.x + v.y * v.y + v.z * v.z;
	Vector tmp;
	tmp.x = v.x * invsqrtf(dot);
	tmp.y = v.y * invsqrtf(dot);
	tmp.z = v.z * invsqrtf(dot);
	//printf("tmp: %f, %f, %f \n", tmp.x, tmp.y, tmp.z);
	return tmp;
}

typedef struct Intersection
{
	enum { NONE, SPHERE, TRIANGLE } objectType;			// type of object intersected with

	Point pos;											// point of intersection
	Vector normal;										// normal at point of intersection
	float viewProjection;								// view projection 
	bool insideObject;									// whether or not inside an object

	__global Material* material;									// material of object

														// object collided with

	__global struct Sphere* sphere;
	__global struct Triangle* triangle;
} Intersection;

inline bool isSphereIntersected(__global Sphere* s, const Ray* r, float* t)
{
	const float EPSILON = 0.01f;

	// Intersection of a ray and a sphere, check the articles for the rationale

	Vector dist;
	dist.x = s->pos.x - r->start.x;
	dist.y = s->pos.y - r->start.y;
	dist.z = s->pos.z - r->start.z;


	float B = r->dir.x * dist.x + r->dir.y * dist.y + r->dir.z * dist.z;

	float D = B * B - (dist.x * dist.x + dist.y * dist.y + dist.z * dist.z) + s->size * s->size;

	// if D < 0, no intersection, so don't try and calculate the point of intersection
	if (D < 0.0f) return false;

	// calculate both intersection times(/distances)
	float t0 = B - sqrt(D);
	float t1 = B + sqrt(D);

	// check to see if either of the two sphere collision points are closer than time parameter
	if ((t0 > EPSILON) && (t0 < *t))
	{
		*t = t0;
		return true;
	}
	else if ((t1 > EPSILON) && (t1 < *t))
	{
		*t = t1;
		return true;
	}

	return false;
}

inline bool isTriangleIntersected(__global Triangle* tri, const Ray* r, float* t)
{
	const float EPSILON = 0.01f;

	Vector e1;
	e1.x = tri->p2.x - tri->p1.x;
	e1.y = tri->p2.y - tri->p1.y;
	e1.z = tri->p2.z - tri->p1.z;

	Vector e2;
	e2.x = tri->p3.x - tri->p1.x;
	e2.y = tri->p3.y - tri->p1.y;
	e2.z = tri->p3.z - tri->p1.z;

	Vector h;
	h.x = r->dir.y * e2.z - r->dir.z * e2.y;
	h.y = r->dir.z * e2.x - r->dir.x * e2.z;
	h.z = r->dir.x * e2.y - r->dir.y * e2.x;

	float det = e1.x * h.x + e1.y * h.y + e1.z * h.z;

	if (det > -EPSILON && det < EPSILON) return false;

	float invDet = 1.0f / det;

	Vector s;
	s.x= r->start.x - tri->p1.x;
	s.y = r->start.y - tri->p1.y;
	s.z = r->start.z - tri->p1.z;

	float u = invDet * (s.x * h.x + s.y * h.y + s.z * h.z);

	if (u < 0.0f || u > 1.0f) return false;

	Vector q;
	q.x = s.y * e1.z - s.z * e1.y;
	q.y = s.z * e1.x - s.x * e1.z;
	q.z = s.x * e1.y - s.y * e1.x;

	float v = invDet * (q.x * r->dir.x + q.y * r->dir.y + q.z * r->dir.z);

	

	if (v < 0.0f || u + v > 1.0f) return false;

	float t0 = invDet * (e2.x * q.x + e2.y * q.y + e2.z * q.z);


	if (t0 > EPSILON && t0 < *t)
	{
		*t = t0;
		return true;
	}

	return false;
}

inline bool objectIntersection(__global Scene* scene, const Ray* viewRay, Intersection* intersect)
{
	const float MAX_RAY_DISTANCE = 2000000.0f;

	// set default distance to be a long long way away
	float t = MAX_RAY_DISTANCE;

	// no intersection found by default
	intersect->objectType = 0;

	// search for sphere collisions, storing closest one found
	for (unsigned int i = 0; i < scene->numSpheres; ++i)
	{
		if (isSphereIntersected(&scene->sphereContainer[i], viewRay, &t))
		{
			intersect->objectType = 1;
			intersect->sphere = &scene->sphereContainer[i];
		}
	}

	//search for triangle collisions, storing closest one found
	for (unsigned int i = 0; i < scene->numTriangles; ++i)
	{
		if (isTriangleIntersected(&scene->triangleContainer[i], viewRay, &t))
		{
			intersect->objectType = TRIANGLE;
			intersect->triangle = &scene->triangleContainer[i];
		}
	}
	
	//// nothing detected, return false
	if (intersect->objectType == 0)
	{
		return false;
	}
	////
	////// calculate the point of the intersection
	intersect->pos.x = viewRay->start.x + viewRay->dir.x * t;
	intersect->pos.y = viewRay->start.y + viewRay->dir.y * t;
	intersect->pos.z = viewRay->start.z + viewRay->dir.z * t;

	

	return true;
}

// calculate collision normal, viewProjection, object's material, and test to see if inside collision object
inline void calculateIntersectionResponse(__global Scene* scene, const Ray* viewRay, Intersection* intersect)
{
	Vector tmp;
	switch (intersect->objectType)
	{
		case 1:
			tmp.x = intersect->pos.x - intersect->sphere->pos.x;
			tmp.y = intersect->pos.y - intersect->sphere->pos.y;
			tmp.z = intersect->pos.z - intersect->sphere->pos.z;
			
			intersect->normal = normalise(tmp);
			intersect->material = &scene->materialContainer[intersect->sphere->materialId];
			break;
		case 2:
			intersect->normal = intersect->triangle->normal;
			intersect->material = &scene->materialContainer[intersect->triangle->materialId];
	}

	

	// calculate view projection
	intersect->viewProjection = viewRay->dir.x * intersect->normal.x + viewRay->dir.y * intersect->normal.y + viewRay->dir.z * intersect->normal.z;

	// detect if we are inside an object (needed for refraction)
	float cal = intersect->normal.x * viewRay->dir.x + intersect->normal.y * viewRay->dir.y + intersect->normal.z * viewRay->dir.z;
	intersect->insideObject = (cal > 0.0f);

	// if inside an object, reverse the normal
	if (intersect->insideObject)
	{
		intersect->normal.x = intersect->normal.x * -1.0f;
		intersect->normal.y = intersect->normal.y * -1.0f;
		intersect->normal.z = intersect->normal.z * -1.0f;
	}
}

inline Colour applySpecular(const Ray* lightRay, __global Light* currentLight, const float fLightProjection, const Ray* viewRay, const Intersection* intersect)
{
	Vector blinnDir;
	blinnDir.x = lightRay->dir.x - viewRay->dir.x;
	blinnDir.y = lightRay->dir.y - viewRay->dir.y;
	blinnDir.z = lightRay->dir.z - viewRay->dir.z;

	float tmp = blinnDir.x * blinnDir.x + blinnDir.y * blinnDir.y + blinnDir.z * blinnDir.z;
	float blinn = invsqrtf(tmp) * max(fLightProjection - intersect->viewProjection, 0.0f);
	blinn = pow(blinn, intersect->material->power);

	Colour output;
	output.red = blinn * intersect->material->specular.red * currentLight->intensity.red;
	output.blue = blinn * intersect->material->specular.blue * currentLight->intensity.blue;
	output.green = blinn * intersect->material->specular.green * currentLight->intensity.green;

	return output;
}

inline bool isInShadow(__global Scene* scene, const Ray* lightRay, const float lightDist)
{
	float t = lightDist;

	

	// search for sphere collision
	for (unsigned int i = 0; i < scene->numSpheres; ++i)
	{
		if (isSphereIntersected(&scene->sphereContainer[i], lightRay, &t))
		{
			
			return true;
		}
	}

	// search for triangle collision
	//for (unsigned int i = 0; i < scene->numTriangles; ++i)
	//{
	//	if (isTriangleIntersected(&scene->triangleContainer[i], lightRay, &t))
	//	{
	//		
	//		return true;
	//	}
	//}

	return false;
}

inline Colour applyCheckerboard(const Intersection* intersect)
{
	Point p;
	p.x = (intersect->pos.x - intersect->material->offset.x) / intersect->material->size;
	p.y = (intersect->pos.y - intersect->material->offset.y) / intersect->material->size;
	p.z = (intersect->pos.z - intersect->material->offset.z) / intersect->material->size;

	int tmpX = convert_int(floor(p.x));
	int tmpY = convert_int(floor(p.y));
	int tmpZ = convert_int(floor(p.z));

	int which = (tmpX + tmpY + tmpZ) & 1;

	return (which ? intersect->material->diffuse : intersect->material->diffuse2);
}

inline Colour applyCircles(const Intersection* intersect)
{
	Point p;
	p.x = (intersect->pos.x - intersect->material->offset.x) / intersect->material->size;
	p.y = (intersect->pos.x - intersect->material->offset.y) / intersect->material->size;
	p.z = (intersect->pos.x - intersect->material->offset.z) / intersect->material->size;

	float tmp = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
	int tmp2 = convert_int(floor(tmp));
	int which = tmp2 & 1;
	
	return (which ? intersect->material->diffuse : intersect->material->diffuse2);
}

inline Colour applyWood(const Intersection* intersect)
{
	Point p;
	p.x = (intersect->pos.x - intersect->material->offset.x) / intersect->material->size;
	p.y = (intersect->pos.y - intersect->material->offset.y) / intersect->material->size;
	p.z = (intersect->pos.z - intersect->material->offset.z) / intersect->material->size;

	// squiggle up where the point is
	//p = { p.x * cosf(p.y * 0.996f) * sinf(p.z * 1.023f), cosf(p.x) * p.y * sinf(p.z * 1.211f), cosf(p.x * 1.473f) * cosf(p.y * 0.795f) * p.z };
	p.x = p.x * cos(p.y * 0.996f) * sin(p.z * 1.023f);
	p.y = cos(p.x) * p.y * sin(p.z * 1.211f);
	p.z = cos(p.x * 1.473f) * cos(p.y * 0.795f) * p.z;


	float tmp = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
	int tmp2 = convert_int(floor(tmp));
	int which = tmp2 & 1;

	return (which ? intersect->material->diffuse : intersect->material->diffuse2);
}

inline Colour applyDiffuse(const Ray* lightRay, __global Light* currentLight, const Intersection* intersect)
{
	Colour tmp;


	switch (intersect->material->type)
	{
	case GOURAUD:
		tmp = intersect->material->diffuse;
		break;
	
	}

	float lambert = lightRay->dir.x * intersect->normal.x + lightRay->dir.y * intersect->normal.y + lightRay->dir.z * intersect->normal.z;

	//printf("%f %f %f\n", intersect->normal.x, intersect->normal.y, intersect->normal.z);

	Colour output;
	output.red = lambert * currentLight->intensity.red * tmp.red;
	output.green = lambert * currentLight->intensity.green * tmp.green;
	output.blue = lambert * currentLight->intensity.blue * tmp.blue;

	return output;
}


inline Colour applyLighting(__global Scene* scene, const Ray* viewRay, const Intersection* intersect)
{
	// colour to return (starts as black)
	Colour output;
	output.red = 0.0f;
	output.blue = 0.0f;
	output.green = 0.0f;

	// same starting point for each light ray
	Ray lightRay = { intersect->pos };

	for (unsigned int j = 0; j < scene->numLights; ++j)
	{
		// get reference to current light
		__global Light* currentLight = &scene->lightContainer[j];

		// light ray direction need to equal the normalised vector in the direction of the current light
		// as we need to reuse all the intermediate components for other calculations, 
		// we calculate the normalised vector by hand instead of using the normalise function
		lightRay.dir.x = currentLight->pos.x - intersect->pos.x;
		lightRay.dir.y = currentLight->pos.y - intersect->pos.y;
		lightRay.dir.z = currentLight->pos.z - intersect->pos.z;


		float angleBetweenLightAndNormal = lightRay.dir.x * intersect->normal.x + lightRay.dir.y * intersect->normal.y + lightRay.dir.z * intersect->normal.z;

		// skip this light if it's behind the object (ie. both light and normal pointing in the same direction)
		if (angleBetweenLightAndNormal <= 0.0f)
		{
			continue;
		}

		// distance to light from intersection point (and it's inverse)
		float lightDist = lightRay.dir.x * lightRay.dir.x + lightRay.dir.y * lightRay.dir.y + lightRay.dir.z * lightRay.dir.z;
		float invLightDist = 1.0f / sqrt(lightDist);

		// light ray projection
		float lightProjection = invLightDist * angleBetweenLightAndNormal;

		lightRay.dir.x = lightRay.dir.x * invLightDist;
		lightRay.dir.y = lightRay.dir.y * invLightDist;
		lightRay.dir.z = lightRay.dir.z * invLightDist;

		if (!isInShadow(scene, &lightRay, lightDist))
		{
			// add diffuse lighting from colour / texture
			output.red += applyDiffuse(&lightRay, currentLight, intersect).red;
			output.blue += applyDiffuse(&lightRay, currentLight, intersect).blue;
			output.green += applyDiffuse(&lightRay, currentLight, intersect).green;

			// add specular lighting
			output.red += applySpecular(&lightRay, currentLight, lightProjection, viewRay, intersect).red;
			output.blue += applySpecular(&lightRay, currentLight, lightProjection, viewRay, intersect).blue;
			output.green += applySpecular(&lightRay, currentLight, lightProjection, viewRay, intersect).green;
		}
	}

	return output;
}

inline Ray calculateReflection(const Ray* viewRay, const Intersection* intersect)
{
	// reflect the viewRay around the object's normal
	Vector tmp;
	tmp.x = viewRay->dir.x - (intersect->normal.x * intersect->viewProjection * 2.0f);
	tmp.y = viewRay->dir.y - (intersect->normal.y * intersect->viewProjection * 2.0f);
	tmp.z = viewRay->dir.z - (intersect->normal.z * intersect->viewProjection * 2.0f);
	Ray newRay = { intersect->pos, tmp};

	return newRay;
}

inline Ray calculateRefraction(const Ray* viewRay, const Intersection* intersect, float* currentRefractiveIndex)
{
	const float DEFAULT_REFRACTIVE_INDEX = 1.0f;

	// change refractive index depending on whether we are in an object or not
	float oldRefractiveIndex = *currentRefractiveIndex;
	*currentRefractiveIndex = intersect->insideObject ? DEFAULT_REFRACTIVE_INDEX : intersect->material->density;

	// calculate refractive ratio from old index and current index
	float refractiveRatio = oldRefractiveIndex / *currentRefractiveIndex;

	// Here we take into account that the light movement is symmetrical from the observer to the source or from the source to the oberver.
	// We then do the computation of the coefficient by taking into account the ray coming from the viewing point.
	float fCosThetaT;
	float fCosThetaI = fabs(intersect->viewProjection);


	// glass-like material, we're computing the fresnel coefficient.
	if (fCosThetaI >= 1.0f)
	{
		// In this case the ray is coming parallel to the normal to the surface
		fCosThetaT = 1.0f;
	}
	else
	{
		float fSinThetaT = refractiveRatio * sqrt(1 - fCosThetaI * fCosThetaI);

		// Beyond the angle (1.0f) all surfaces are purely reflective
		fCosThetaT = (fSinThetaT * fSinThetaT >= 1.0f) ? 0.0f : sqrt(1 - fSinThetaT * fSinThetaT);
	}

	Vector tmp;
	tmp.x = (viewRay->dir.x + intersect->normal.x * fCosThetaI) * refractiveRatio - (intersect->normal.x * fCosThetaT);
	tmp.y = (viewRay->dir.y + intersect->normal.y * fCosThetaI) * refractiveRatio - (intersect->normal.y * fCosThetaT);
	tmp.z = (viewRay->dir.z + intersect->normal.z * fCosThetaI) * refractiveRatio - (intersect->normal.z * fCosThetaT);

	// Here we compute the transmitted ray with the formula of Snell-Descartes
	Ray newRay = { intersect->pos, tmp };

	return newRay;
}

inline Colour traceRay(__global Scene* scene, Ray viewRay)
{
	Colour output;
	output.red = 0.0f;
	output.blue = 0.0f;
	output.green = 0.0f;
	const float DEFAULT_REFRACTIVE_INDEX = 1.0f;

	Colour white;
	white.red = 1.0f;
	white.green = 1.0f;
	white.blue = 1.0f;

	float currentRefractiveIndex = DEFAULT_REFRACTIVE_INDEX;		// current refractive index
	float coef = 1.0f;												// amount of ray left to transmit
	
	const int MAX_RAYS_CAST = 10;

	Intersection intersect;

	for (int level = 0; level < MAX_RAYS_CAST; ++level)
	{
		//// check for intersections between the view ray and any of the objects in the scene
		//// exit the loop if no intersection found
		if (!objectIntersection(scene, &viewRay, &intersect)) break;
		//
		//// calculate response to collision: ie. get normal at point of collision and material of object
		calculateIntersectionResponse(scene, &viewRay, &intersect);
		//
		//// apply the diffuse and specular lighting 
		
		if (!intersect.insideObject) {
			output.red += coef * applyLighting(scene, &viewRay, &intersect).red;
			output.green += coef * applyLighting(scene, &viewRay, &intersect).green;
			output.blue += coef * applyLighting(scene, &viewRay, &intersect).blue;
		}
		//
		//// if object has reflection or refraction component, adjust the view ray and coefficent of calculation and continue looping
		if (intersect.material->reflection)
		{
			viewRay = calculateReflection(&viewRay, &intersect);
			coef *= intersect.material->reflection;
		}
		else if (intersect.material->refraction)
		{
			viewRay = calculateRefraction(&viewRay, &intersect, &currentRefractiveIndex);
			coef *= intersect.material->refraction;
		}
		else
		{
			// if no reflection or refraction, then finish looping (cast no more rays)
			return output;
		}
	}

	if (coef > 0.0f)
	{
		output.red += coef * scene->materialContainer[scene->skyboxMaterialId].diffuse.red;
		output.green += coef * scene->materialContainer[scene->skyboxMaterialId].diffuse.green;
		output.blue += coef * scene->materialContainer[scene->skyboxMaterialId].diffuse.blue;
	}

	return output;
}



__kernel void GetDateFromScene(__global struct Scene* scene,
	__global struct Point* cameraPosition,
	__global struct Material* materialContainer,
	__global struct Sphere* sphereContainer,
	__global struct Triangle* triangleContainer,
	__global struct Light* lightContainer,
	__global int* aaLevel,
	__global unsigned int* out)
{
	unsigned int width = get_global_size(0);
	unsigned int height = get_global_size(1);
	unsigned int ix = get_global_id(0);
	unsigned int iy = get_global_id(1);

	scene->materialContainer = materialContainer;
	scene->sphereContainer = sphereContainer;
	scene->triangleContainer = triangleContainer;
	scene->lightContainer = lightContainer;

	Colour output;
	output.red = 0.0f;
	output.blue = 0.0f;
	output.green = 0.0f;

	const float PIOVER180 = 0.017453292519943295769236907684886f;

	const float dirStepSize = 1.0f / (0.5f * width / tan(PIOVER180 * 0.5f * scene->cameraFieldOfView));

	const float sampleStep = 1.0f / *aaLevel, sampleRatio = 1.0f / (*aaLevel * *aaLevel);


	// loop through all sub-locations within the pixel
	//if (ix == 0 && iy == 0) {
	//	printf("%d   ", ix - (width / 2));
	//	int i = ix - (width / 2);
	//	printf("%f \n", convert_float( i ));
	//}

	int iWidth = ix - (width / 2);
	int iHeigh = iy - (height / 2);

	for (float fragmentx = convert_float(iWidth); fragmentx < iWidth + 1.0f; fragmentx += sampleStep)
	{
		
		for (float fragmenty = convert_float(iHeigh); fragmenty < iHeigh + 1.0f; fragmenty += sampleStep)
		{
			// direction of default forward facing ray
			Vector dir;
			dir.x = fragmentx * dirStepSize;
			dir.y = fragmenty * dirStepSize;
			dir.z = 1.0f;

			

			// rotated direction of ray
			Vector rotatedDir;

			rotatedDir.x = dir.x * cos(scene->cameraRotation) - dir.z * sin(scene->cameraRotation);
			rotatedDir.y = dir.y;
			rotatedDir.z = dir.x * sin(scene->cameraRotation) + dir.z * cos(scene->cameraRotation);

			// view ray starting from camera position and heading in rotated (normalised) direction
			Ray viewRay;
			viewRay.start = scene->cameraPosition;
			viewRay.dir = normalise(rotatedDir);

			

			//= { scene->cameraPosition, normalise(rotatedDir) };

			// follow ray and add proportional of the result to the final pixel colour
			
			output.red += sampleRatio * traceRay(scene, viewRay).red;
			
			output.green += sampleRatio * traceRay(scene, viewRay).green;
			
			output.blue += sampleRatio * traceRay(scene, viewRay).blue;

			//printf("%f \n", traceRay(scene, viewRay).blue);

		}
	}

	out[width * iy + ix] = convertToPixel(output, scene->exposure);
}
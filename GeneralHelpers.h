
#ifndef  GENERAL_HELPERS
#define  GENERAL_HELPERS

#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"

#define PI (float) 3.14159265358979323846

//----------------------------------------------------------------------------------------------------------------------------------------------------------------//
// SMALL HELPERS	||
//------------------//

#define KA				GZ_AMBIENT_COEFFICIENT	
#define KD				GZ_DIFFUSE_COEFFICIENT	
#define KS				GZ_SPECULAR_COEFFICIENT	
#define SPECULAR_POWER  GZ_DISTRIBUTION_COEFFICIENT 

#define max(a, b)		((a > b) ? a : b)
#define min(a, b)		((a < b) ? a : b)

#define COPY_COORD(A, B)	A[0] = B[0];\
							A[1] = B[1];\
							A[2] = B[2];

const float degToRad = PI / 180;

const GzMatrix I_4X4 = { 1, 0, 0, 0,
0, 1, 0, 0,
0, 0, 1, 0,
0, 0, 0, 1 };

template< typename type >
inline type clampvalue(type val, type min, type max)
{
	return (val > max) ? max : (val < min ? min : val);
}

void CrossProduct(float a[3], float b[3], float productVector[3]);

void NormalizeVector(float vec[3]);

void MatrixMultiply(GzMatrix matA, GzMatrix matB, GzMatrix matRes);

void MatrixMulScalar(GzMatrix mat, float scalar);

float Vec4Dot(float a[4], float b[4]);

float Vec3Dot(float a[3], float b[3]);

void CopyMatrix(const GzMatrix from, GzMatrix to);

void TransformCoord(GzCoord &positionVector, GzMatrix transformationMatrix);

void CoordScalarMul(GzCoord coord, float scalar);

void CoordDiff(GzCoord coord1, GzCoord coord2);

void ConvertToPureRotationMatrix(GzMatrix mat);

#endif // !GENERALHELPERS
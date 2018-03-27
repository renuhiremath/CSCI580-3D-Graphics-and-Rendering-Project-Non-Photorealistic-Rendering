#include	"GeneralHelpers.h"
#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"

const GzMatrix I_4X4 = { 1, 0, 0, 0,
0, 1, 0, 0,
0, 0, 1, 0,
0, 0, 0, 1 };

template< typename type >
inline type clampvalue(type val, type min, type max)
{
	return (val > max) ? max : (val < min ? min : val);
}

void CrossProduct(float a[3], float b[3], float productVector[3])
{
	// productVector will have the components of aXb
	productVector[0] = a[1] * b[2] - a[2] * b[1];
	productVector[1] = a[2] * b[0] - a[0] * b[2];
	productVector[2] = a[0] * b[1] - a[1] * b[0];
}

void NormalizeVector(float vec[3])
{
	float sumOfSquares = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
	sumOfSquares = sqrt(sumOfSquares);
	vec[0] = vec[0] / sumOfSquares;
	vec[1] = vec[1] / sumOfSquares;
	vec[2] = vec[2] / sumOfSquares;
}

void MatrixMultiply(GzMatrix matA, GzMatrix matB, GzMatrix matRes)
{
	for (int row = 0; row < 4; ++row)
	{
		for (int col = 0; col < 4; ++col)
		{
			matRes[row][col] = 0;
			for (int it = 0; it < 4; ++it)
			{
				matRes[row][col] += (matA[row][it] * matB[it][col]);
			}
		}
	}
}

void MatrixMulScalar(GzMatrix mat, float scalar)
{
	for (int row = 0; row < 4; ++row)
		for (int col = 0; col < 4; ++col)
			mat[row][col] = mat[row][col] * scalar;
}

float Vec4Dot(float a[4], float b[4])
{
	float dot = 0;
	for (int it = 0; it < 4; ++it)
	{
		dot += a[it] * b[it];
	}
	return dot;
}

float Vec3Dot(float a[3], float b[3])
{
	float dot = 0;
	for (int it = 0; it < 3; ++it)
	{
		dot += a[it] * b[it];
	}
	return dot;
}

void CopyMatrix(const GzMatrix from, GzMatrix to)
{
	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j)
			to[i][j] = from[i][j];
}

void TransformCoord(GzCoord &positionVector, GzMatrix transformationMatrix)
{
	float vecHomogenous[4] = { positionVector[0], positionVector[1], positionVector[2], 1 };
	float w = Vec4Dot(vecHomogenous, transformationMatrix[3]);
	positionVector[0] = Vec4Dot(vecHomogenous, transformationMatrix[0]) / w;
	positionVector[1] = Vec4Dot(vecHomogenous, transformationMatrix[1]) / w;
	positionVector[2] = Vec4Dot(vecHomogenous, transformationMatrix[2]) / w;
}

void CoordScalarMul(GzCoord coord, float scalar)
{
	for (int i = 0; i < 3; ++i)
		coord[i] *= scalar;
}

void CoordDiff(GzCoord coord1, GzCoord coord2)
{
	for (int i = 0; i < 3; ++i)
		coord1[i] = coord1[i] - coord2[i];
}

void ConvertToPureRotationMatrix(GzMatrix mat)
{
	// assuming uniform scale
	// remove scale
	float scale = sqrt(mat[0][0] * mat[0][0] + mat[0][1] * mat[0][1] + mat[0][2] * mat[0][2]);
	if (scale != 1)
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				mat[i][j] = mat[i][j] / scale;

	// remove translation
	mat[0][3] = 0;
	mat[1][3] = 0;
	mat[2][3] = 0;
}

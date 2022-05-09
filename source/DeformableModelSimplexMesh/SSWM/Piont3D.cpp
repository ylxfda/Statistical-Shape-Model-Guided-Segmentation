
#include "Piont3D.h"

Point3D::Point3D()
{
	this->x = 0;
	this->y = 0;
	this->z = 0;
	this->X = 0;
	this->Y = 0;
	this->Z = 0;
	this->theta = 0.0;
	this->phi = 0.0;
	this->MeanX = 0.0; 
	this->StdX = 0.0;
	this->MeanY = 0.0; 
	this->StdY = 0.0;
	this->MeanZ = 0.0; 
	this->StdZ = 0.0;
	this->OpiScores[0]=0.0;
	this->OpiScores[1]=0.0;
	this->OpiScores[2]=0.0;
	this->Index = 0;
}

Point3D::Point3D(float x, float y, float z)
{
	this->x = x;
	this->y = y;
	this->z = z;
}

Point3D::~Point3D()
{
}

// operator overloading

// add two Points

Point3D Point3D::operator+(const Point3D & b) const

{

	return Point3D(x + b.x, y + b.y, z+b.z);

}



// subtract point b from a

Point3D Point3D::operator-(const Point3D & b) const

{

	return Point3D(x - b.x, y - b.y, z-b.z);

}

// reverse sign of Vector

Point3D Point3D::operator-() const

{

	return Point3D(-x, -y, -z);

}



// multiple vector by n

Point3D Point3D::operator*(float n) const

{

	return Point3D(n * x, n * y, n * z);

}

// divide vector by n

Point3D Point3D::operator/(float n) const

{

	return Point3D(x/n, y/n, z/n);

}

void Point3D::SetLocation(double x, double y, double z)
{
	this->x = x;
	this->y = y;
	this->z = z;
}

Point3D& Point3D::operator =(Point3D &temp)
{
	this->SetLocation(temp.x, temp.y, temp.z);
	return *this;
}
#ifndef Point3D_h  
#define Point3D_h       1

#define NUM_Sample 20
class Point3D
{
public:
	float x;   //the values for computation
	float y;
	float z;
	float X;   //the values for backup
	float Y;
	float Z;
	float theta; //the polar coordinates of this point
	float phi;
	float MeanX; //mean position and standard deviation
	float StdX;
	float MeanY; //mean position and standard deviation
	float StdY;
	float MeanZ; //mean position and standard deviation
	float StdZ;
	float PriorPC[3][3]; //the principle component of prior model, every row is a pc.
	float PriorMEAN[3]; //the mean in the prior.
	float PriorMINScore[3]; //the min score in the prior.
	float PriorMAXScore[3]; //the max in the prior.
	float PriorSTD[3]; //the std in the prior.
	float OpiScores[3]; //the scores obtained through optimization.
	float Scores[NUM_Sample][3]; //scores of in the form of components from the sample set; 
	

	// operator overloading

	Point3D operator+(const Point3D & b) const;

	Point3D operator-(const Point3D & b) const;

	Point3D operator-() const;

	Point3D operator*(float n) const;

	Point3D operator/(float n) const;

	Point3D& operator =(Point3D &temp);

	long Index;
	
	Point3D();
	Point3D(float x, float y, float z);
	~Point3D();

	void SetLocation(double a, double b, double c);
};

#endif
/**
 * Copyright (C) 2024 Lehigh University.
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see http://www.gnu.org/licenses/.
 *
 * Author: David Rutkowski (dmr518@lehigh.edu)
 */
 
#include <math.h>
#include <vector>
#include <queue>
#include <deque>
#include <fstream>
#include <iostream>
#include <random>
#include <time.h>
#include <sstream>


const double pi = 3.14159265359;

struct Coordinate2D
{
	double x,z;
};

struct Coordinate
{
	double x, y, z;
	
	struct Coordinate getUnitCoord()
	{
		double mag = sqrt(x*x + y*y + z*z);
		
		double invMag = 1.0/mag;
		
		struct Coordinate uCoord = {x*invMag, y*invMag, z*invMag};
		
		return uCoord;
	}
	
	double getMagnitude()
	{
		return sqrt(x*x + y*y + z*z);
	}
	
	Coordinate crossProduct(Coordinate b)
	{
		Coordinate c;
		
		c.x = this->y * b.z - this->z * b.y;
		c.y = this->z * b.x - this->x * b.z;
		c.z = this->x * b.y - this->y * b.x;
		
		return c;
	}
	
	
	Coordinate operator+(const Coordinate& b)
	{
		Coordinate c;
		c.x = this->x + b.x;
		c.y = this->y + b.y;
		c.z = this->z + b.z;
		return c;
	}
	
	Coordinate operator-(const Coordinate& b)
	{
		Coordinate c;
		c.x = this->x - b.x;
		c.y = this->y - b.y;
		c.z = this->z - b.z;
		return c;
	}
	
	Coordinate operator-() const
	{
		Coordinate c;
		
		c.x = -this->x;
		c.y = -this->y;
		c.z = -this->z;
		return c;
	}
	
	double operator*(Coordinate b) const
	{
		return this->x * b.x + this->y * b.y + this->z * b.z;
	}
	
	Coordinate operator*(double a) const
	{
		Coordinate c;
		c.x = this->x * a;
		c.y = this->y * a;
		c.z = this->z * a;
		return c;
	}
	
	Coordinate operator/(double a) const
	{
		Coordinate c;
		c.x = this->x / a;
		c.y = this->y / a;
		c.z = this->z / a;
		return c;
	}
};

inline Coordinate operator*(const double a, const Coordinate& b)
{
	Coordinate c = b * a;
	return c;
}

struct Particle
{
    Coordinate pos;
    int name;
    //int type;
    double timeLastGlobalUpdate;
    double posTime;
    
    int secondname;
    
    Particle(Coordinate c, int n, double t, double pt)
        : pos(c)
        , name(n)
        , timeLastGlobalUpdate(t)
        , posTime(pt)
        {}
};

// periodic wrap for positions varying between 0 and length
double periodicWrap(double pos, double length)
{
	pos = pos - length * 0.5;
	return pos - length * round(pos / length) + length*0.5;
}

double calcDistance(Coordinate posi, Coordinate posj)
{
    Coordinate rij = posj - posi;
    
    return rij.getMagnitude();
}

double calcDistanceSphere(Coordinate posi, Coordinate posj, double sphereRadius)
{
    double currDistance = calcDistance(posi, posj);
    
    double asin_argument = 0.5* currDistance / sphereRadius;
    
    // this needs to be here for rounding error issues, if it is not, then some points are included as exo/endo that shouldnt be
    if(asin_argument > 1.0 && asin_argument < 1.01)
    {
        asin_argument = 1.0;
    }
    
    double arclength = 2.0*sphereRadius * asin(asin_argument);
    
    return arclength;
}

double calcDistanceCylinder(Coordinate posi, Coordinate posj, double radius)
{
	double posi_z = posi.z;
	double posj_z = posj.z;
	
	posi.z = 0.0;
	posj.z = 0.0;
	
	double currDistance_xy = calcDistance(posi, posj);
    
	double asin_argument = 0.5* currDistance_xy / radius;
    
    // this needs to be here for rounding error issues, if it is not, then some points are included as exo/endo that shouldnt be
    if(asin_argument > 1.0 && asin_argument < 1.01)
    {
        asin_argument = 1.0;
    }
	
	double arclength_xy = 2.0 * radius * asin(asin_argument);
	
	double length_z = abs(posi_z - posj_z);
	
	return sqrt(arclength_xy*arclength_xy + length_z*length_z);	
}


double calcArcLengthDistance(Coordinate posi, Coordinate posj, double sphereRadius)
{
	if(posi.z > 0.0 && posj.z > 0.0)
	{
		// both points are on the sphere
		return calcDistanceSphere(posi, posj, sphereRadius);
	}
	else if(posi.z <= 0.0 && posj.z <= 0.0)
	{
		// both points are on the cylinder
		return calcDistanceCylinder(posi, posj, sphereRadius);
	}
	else
	{
		// one of the points is on the sphere and one is on the cylinder
		if(posi.z > 0.0)
		{
			// project this point onto the cylinder
			double initialPhi = atan2(posi.y, posi.x);
			
			posi = Coordinate{sphereRadius*cos(initialPhi), sphereRadius*sin(initialPhi), posi.z};
		}
		else if(posj.z > 0.0)
		{
			// project this point onto the cylinder
			double initialPhi = atan2(posj.y, posj.x);
			
			posj = Coordinate{sphereRadius*cos(initialPhi), sphereRadius*sin(initialPhi), posj.z};
		}
		
		return calcDistanceCylinder(posi, posj, sphereRadius);
	}	
}
	

Coordinate diffuseParticleSphere(Coordinate pos, double delPosX, double delPosY, double radius, double cylinderLength);
Coordinate diffuseParticleCylinder(Coordinate pos, double delPosX, double delPosY, double radius, double cylinderLength);
Coordinate diffuseParticleSphere_fromCyl(Coordinate pos, Coordinate tangent_3d, double length, double radius);

Coordinate diffuseParticle(Coordinate pos, double delPosX, double delPosY, double radius, double cylinderLength)
{
	if(pos.z > 0.0)
	{
		return diffuseParticleSphere(pos, delPosX, delPosY, radius, cylinderLength);
	}
	else
	{		
		return diffuseParticleCylinder(pos, delPosX, delPosY, radius, cylinderLength);
	}
	
}



Coordinate diffuseParticleSphere(Coordinate pos, double delPosX, double delPosY, double radius, double cylinderLength)
{
    double initialTheta = acos(pos.z / radius);
    double initialPhi = atan2(pos.y, pos.x);
    
	/*delPosX = 0.0;
	delPosY = 0.01;*/
    
    double stepsize = sqrt(delPosX*delPosX + delPosY*delPosY);				
    
    // hankbesser on github
    double theta = stepsize / radius;
    double phi = atan2(delPosY, delPosX);
	/*double phi = pi*0.5;
	
	phi = acos((cos(theta) - 1.0) / (tan(initialTheta)*sin(theta)));
	
	while(phi > pi)
	{
		phi = phi - 2.0*pi;
	}
	
	while(phi < -pi)
	{
		phi = phi + 2.0*pi;
	}*/
	
	//std::cout << phi << std::endl;
    
    // this is rotated position assuming that the particle started at the north pole (0,0,radius)
    struct Coordinate final_pos3D;
    final_pos3D.x = radius*cos(phi)*sin(theta);
    final_pos3D.y = radius*sin(phi)*sin(theta);
    final_pos3D.z = radius*cos(theta);
    
    struct Coordinate rotK;
    rotK.x = 0.0;
    rotK.y = 1.0;
    rotK.z = 0.0;
    
	
    // rotate by initialTheta
    double currDotProduct = rotK*final_pos3D;
    struct Coordinate currCrossProduct = rotK.crossProduct(final_pos3D);
	
	//Coordinate testPos = Coordinate {cos(initialTheta)*final_pos3D.x + sin(initialTheta)*final_pos3D.z, final_pos3D.y, -sin(initialTheta)*final_pos3D.x + cos(initialTheta)*final_pos3D.z};
    
    final_pos3D = final_pos3D*cos(initialTheta) + currCrossProduct*sin(initialTheta) + rotK*currDotProduct*(1.0 - cos(initialTheta));
	
	//std::cout << testPos.x << " " << testPos.y << " " << testPos.z << std::endl;
	//std::cout << final_pos3D.x << " " << final_pos3D.y << " " << final_pos3D.z << std::endl << std::endl;
    
    rotK.x = 0.0;
    rotK.y = 0.0;
    rotK.z = 1.0;
    
    // rotate by initialPhi
    currDotProduct = rotK*final_pos3D;
    currCrossProduct = rotK.crossProduct(final_pos3D);
    
	//testPos = Coordinate {cos(initialPhi)*final_pos3D.x - sin(initialPhi)*final_pos3D.y, sin(initialPhi)*final_pos3D.x + cos(initialPhi)*final_pos3D.y, final_pos3D.z};
	
    final_pos3D = final_pos3D*cos(initialPhi) + currCrossProduct*sin(initialPhi) + rotK*currDotProduct*(1.0 - cos(initialPhi));
	
	//std::cout << testPos.x << " " << testPos.y << " " << testPos.z << std::endl;
	//std::cout << final_pos3D.x << " " << final_pos3D.y << " " << final_pos3D.z << std::endl << std::endl;
	
	
	//std::cout << pos.x << " " << pos.y << " " << pos.z << std::endl;
	//std::cout << final_pos3D.x << " " << final_pos3D.y << " " << final_pos3D.z << std::endl << std::endl;
	
	//std::cout << stepsize << " " << calcArcLengthDistance(final_pos3D, pos, radius) << std::endl;
	
	// now need to see if this would diffuse the particle off the hemisphere onto the cylinder part
	// cylinder starts at z = 0
	
	Coordinate r0 = pos.getUnitCoord();
	Coordinate r1 = final_pos3D.getUnitCoord();
	
	double d = r0 * r1;
	
	double beta = 1.0 / sqrt(1.0 - d*d);
	double alpha = -d * beta;
	
	Coordinate w1 = alpha * r0 + beta * r1;
	
	
	
	double atan_value = -atan(r0.z / w1.z);
	
	double tm = atan_value;
	
	if(tm < 0.0)
	{
		tm = tm + pi;
	}
	
	//std::cout << r0.z << " " << w1.z << " " << tm << " " << (cos(tm)*r0 + sin(tm)*w1).z << std::endl;
	
	double d_acos = acos(d);
	
	//std::cout << tm << " " << d << " " << d_acos << " " << alpha << " " << beta << std::endl;
	
	//return final_pos3D;
	
	
	
	if(tm > 1E-10 && tm < d_acos)
	{
		// it crossed the boundary to the cylinder
		double arclength = tm * radius;
		
		double remainingLengthCylinder = stepsize - arclength;
		
		Coordinate tangentVectorAtBoundary = -sin(tm)*r0 + cos(tm)*w1;
		
		//std::cout << tangentVectorAtBoundary.x << " " << tangentVectorAtBoundary.y << " " << tangentVectorAtBoundary.z << " " << tangentVectorAtBoundary.getMagnitude() << std::endl;
		
		Coordinate startingPosCylinder = radius*(cos(tm)*r0 + sin(tm)*w1);
		
		//std::cout << startingPosCylinder.z << std::endl;
		
		
		//double deltaPhi = atan2(direction.y, direction.x);
		double dirX = sqrt(tangentVectorAtBoundary.x*tangentVectorAtBoundary.x + tangentVectorAtBoundary.y*tangentVectorAtBoundary.y);
		double dirZ = tangentVectorAtBoundary.z;
		
		double cyl_delX = remainingLengthCylinder*dirX;
		double cyl_delY = remainingLengthCylinder*dirZ;
		
		// now diffuse on cylinder starting at startingPosCylinder, going in direction tangentVectorAtBoundary for length remainingLengthCylinder
		final_pos3D = diffuseParticleCylinder(startingPosCylinder, cyl_delX, cyl_delY, radius, cylinderLength);
		
	}
    
    return final_pos3D;
}

//Coordinate diffuseParticleCylinder(Coordinate startingPos, Coordinate direction, double stepLength, double radius)
Coordinate diffuseParticleCylinder(Coordinate pos, double delPosX, double delPosY, double radius, double cylinderLength)
{	
	//return startingPos;
    double initialX = (atan2(pos.y, pos.x) + pi) / (2.0*pi) * (2.0*pi*radius);
	double initialZ = pos.z;
	
	//std::cout << initialX << " " << initialZ << std::endl;
	
	
	/*double dirZ = direction.z;
	//double deltaPhi = atan2(direction.y, direction.x);
	double dirX = sqrt(direction.x*direction.x + direction.y*direction.y);
	
	//std::cout << sqrt(dirZ*dirZ + dirX*dirX) << std::endl;*/
	
	double finalX = initialX + delPosX;
	double finalZ = initialZ + delPosY;
	
	
	//std::cout << sqrt((finalX-initialX)*(finalX-initialX) + (finalZ-initialZ)*(finalZ-initialZ)) << " " << stepLength << std::endl;
	while(finalX > 2*pi*radius)
	{
		finalX = finalX - 2*pi*radius;
	}
	
	double finalAngle = finalX / (2.0*pi*radius) * (2.0*pi);	
	
	Coordinate pos_3d = Coordinate {-radius * cos(finalAngle), -radius * sin(finalAngle), finalZ};
	
	//std::cout << pos.x << " " << pos.y << " " << pos.z << std::endl;
	//std::cout << pos_3d.x << " " << pos_3d.y << " " << pos_3d.z << std::endl << std::endl;
	
	//return pos_3d;
	
	
	
	// now check if pos is below reflecting boundary
	if(pos_3d.z < -cylinderLength)
	{
		pos_3d.z = -2.0*cylinderLength - pos_3d.z;
		
		//std::cout << pos_3d.z << std::endl;
		//dirZ = -dirZ;
	}
	
	//now need to check if pos_3d z is above boundary with hemisphere
	if(pos_3d.z > 0.0)
	{
		//it is on the sphere, need to back it up to position at boundary		
		double stepBackLength_Z = finalZ;
		
		double stepBackX = stepBackLength_Z / delPosY * delPosX;
		double stepBackY = stepBackLength_Z / delPosY * delPosY;
		
		double finalSphereX = finalX - stepBackX;
		double finalSphereZ = finalZ - stepBackY;
		
		while(finalSphereX > 2*pi*radius)
		{
			finalSphereX = finalSphereX - 2*pi*radius;
		}
		
		while(finalSphereX < 0)
		{
			finalSphereX = finalSphereX + 2*pi*radius;
		}
		
		double finalCylinderAngle = finalSphereX / (2.0*pi*radius) * (2.0*pi);	
		
		
		//now convert to 3d coordinates
		Coordinate boundary_pos_3d = Coordinate {-radius * cos(finalCylinderAngle), -radius * sin(finalCylinderAngle), finalSphereZ};
		
		
		
		
		
		
		/*double stepBackLength = sqrt(stepBackX*stepBackX + stepBackY*stepBackY);
		
		//Coordinate tangent_3d = Coordinate {-radius * cos(stepBackX/radius), -radius * sin(stepBackX/radius), stepBackY};
		
		Coordinate normal = boundary_pos_3d.getUnitCoord();
		Coordinate z = {0.0, 0.0, 1.0};
		Coordinate planar_x = z.crossProduct(normal);
		
		
		//Coordinate tangent_3d = pos_3d - boundary_pos_3d;
		
		Coordinate tangent_3d = stepBackY*z + stepBackX*planar_x;
		
		//std::cout << delPosX << " " << delPosY << std::endl;
		//std::cout << tangent_3d.x << " " << tangent_3d.y << " " << tangent_3d.z << std::endl << std::endl;
		//exit(0);
		
		//std::cout << tangent_3d.z << std::endl;
		
		pos_3d = diffuseParticleSphere_fromCyl(boundary_pos_3d, tangent_3d, stepBackLength, radius);*/
		
		
		
		
		
		
		
		double tmpAngle = atan2(stepBackY, stepBackX);
		
		double initialTheta = acos(boundary_pos_3d.z / radius);
		
		double step_size = sqrt(stepBackX*stepBackX + stepBackY*stepBackY);
		
		double step_theta = step_size / radius;
		
		double tmpPhi = acos((cos(step_theta) - 1.0) / (tan(initialTheta)*sin(step_theta)));
		
		tmpPhi = tmpPhi + tmpAngle;
		
		while(tmpPhi > pi)
		{
			tmpPhi = tmpPhi - 2.0*pi;
		}
		
		while(tmpPhi < -pi)
		{
			tmpPhi = tmpPhi + 2.0*pi;
		}
		
		double realStepX = step_size * cos(tmpPhi);
		double realStepY = step_size * sin(tmpPhi);
		
		//std::cout << tmpPhi << " " << realStepX << " " << realStepY << " " << stepBackX << " " << stepBackY << std::endl;
		
		//std::cout << boundary_pos_3d.x << " " << boundary_pos_3d.y << " " << boundary_pos_3d.z << std::endl;
		
		//now diffuse this on circle
		Coordinate pos_3d = diffuseParticleSphere(boundary_pos_3d, realStepX, realStepY, radius, cylinderLength);
		
	}
	
	return pos_3d;
	
}


/*Coordinate diffuseParticleSphere_fromCyl(Coordinate pos, Coordinate tangent_3d, double length, double radius)
{
	double angle = length / radius;
	
	//std::cout << pos.x << " " << pos.y << " " << pos.z << std::endl;
	
	Coordinate final_pos = radius*(cos(angle)*pos.getUnitCoord() + sin(angle)*tangent_3d.getUnitCoord());
	
	//std::cout << final_pos.x << " " << final_pos.y << " " << final_pos.z << std::endl << std::endl;
	
	return final_pos;	
}*/

Coordinate randomPointSurface(double radius, double cylinderLength, double randomNums[3])
{
	double areaHemisphere = 0.5 * 4.0 * pi * radius*radius;
	double areaCylinder = 2.0*pi*radius * cylinderLength;
	
	double areaTotal = areaHemisphere + areaCylinder;
	
	// pick either hemisphere or cylinder portion based on area of these regions
	if(areaCylinder / areaTotal < randomNums[0])
	{
		// put point on hemisphere
		double phiRnd = 2.0*pi*randomNums[1];
        double thetaRnd = acos(randomNums[2]);
            
        return Coordinate{radius*cos(phiRnd)*sin(thetaRnd), radius*sin(phiRnd)*sin(thetaRnd), radius*cos(thetaRnd)};
		
	}
	else
	{
		// put point on cylinder
		double phiRnd = 2.0*pi*randomNums[1];
        double zRnd = -cylinderLength * randomNums[2];
            
        return Coordinate{radius*cos(phiRnd), radius*sin(phiRnd), zRnd};
	}
}

Coordinate randomPointCylinder(double radius, double cylinderLength, double randomNums[3])
{
	double areaHemisphere = 0.5 * 4.0 * pi * radius*radius;
	double areaCylinder = 2.0*pi*radius * cylinderLength;
	
	double areaTotal = areaHemisphere + areaCylinder;

	{
		// put point on cylinder
		double phiRnd = 2.0*pi*randomNums[1];
        double zRnd = -cylinderLength * randomNums[2];
            
        return Coordinate{radius*cos(phiRnd), radius*sin(phiRnd), zRnd};
	}
}

double getFinalArcLengthAddArea(double s0, double AreaAdded, double radius)
{
	double quarterCircum = pi*radius*0.5;
	
	if(s0 > quarterCircum)
	{
		// both points will be on cylinder since adding area
		return AreaAdded / (2.0*pi*radius) + s0;
	}
	else
	{
		// first point is on sphere
		
		double acosArgument = AreaAdded / (2.0*pi*radius*radius) + (1.0 - cos(s0 / radius));
		
		if(acosArgument < 1.0)
		{
			if(acosArgument <= 0.0)
			{
				return 0.0;
			}
			
			// then final point will also be on sphere
			return radius * acos(1.0 - acosArgument);
		}
		else
		{
			// then final point will be on cylinder
			return acosArgument*radius + (pi*0.5 - 1.0)*radius;
		}
	}	
}



double appx_erf_imaginary_returnReal(double x, double y, int max_n)
{
    double const_val = 0.0;
    
    if(x != 0)
	{
        const_val = exp(-x*x)/(2.0*pi*x)*(sin(2.0*x*y));
	}
    else
	{
        const_val = y / pi;
	}
    
    //max_n = 10
    
    double summation = 0.0;
	
	for(int n = 1; n < max_n; n++)
	{
        //fn = 2.0*x-2.0*x*math.cosh(n*y)*math.cos(2.0*x*y) + n*math.sinh(n*y)*math.sin(2.0*x*y)
        double gn = 2.0*x*cosh(n*y)*sin(2.0*x*y) + n*sinh(n*y)*cos(2.0*x*y);
        
        double prefactor = exp(-0.25*n*n) / (n*n + 4*x*x);
        
        double val_in_sum = prefactor * (gn);
        
        summation += val_in_sum;
	}
        
    return const_val + 2.0 / pi * exp(-x*x)*summation;
}
	
double avgDisplacementAtArclength(double x)
{
	double disp = 0.0;
	
	double xval = (x - 0.05)/0.1;
	int xswitch = (int)xval;
	if(xval < 0.0)
	{
		xswitch = -1;
	}
	switch(xswitch) {
		case -1:
			disp = 0.0032416760778997966*x + 8.567746110826917e-22;
			break;
		case 0:
			disp = 0.0029986207888715575*x + 1.2152764451411864e-05;
			break;
		case 1:
			disp = 0.00372095189240205*x + -9.619690107816216e-05;
			break;
		case 2:
			disp = 0.0033613831389950125*x + -6.304712726402989e-06;
			break;
		case 3:
			disp = 0.003387725833400645*x + -1.552465576837376e-05;
			break;
		case 4:
			disp = 0.0025429366834025095*x + 0.00036463046173078773;
			break;
		case 5:
			disp = 0.0021025028377034123*x + 0.000606869076865291;
			break;
		case 6:
			disp = 0.0019364666482259423*x + 0.0007147926000256453;
			break;
		case 7:
			disp = 0.0011658225590549475*x + 0.0012927756669038914;
			break;
		case 8:
			disp = 0.0014769056211919817*x + 0.0010283550640874132;
			break;
		case 9:
			disp = 0.0008499669582840446*x + 0.0016239467938499537;
			break;
		case 10:
			disp = -1.2110225816278455e-05*x + 0.002529127837155293;
			break;
		case 11:
			disp = -0.00010746707609490588*x + 0.0026387882149757147;
			break;
		case 12:
			disp = -0.0007136029020452348*x + 0.003396457997413627;
			break;
		case 13:
			disp = -0.001147983766328732*x + 0.003982872164196348;
			break;
		case 14:
			disp = -0.0013442172727054804*x + 0.004267410748442633;
			break;
		case 15:
			disp = -0.0015190819290718672*x + 0.004538450965810533;
			break;
		case 16:
			disp = -0.0016425818064065816*x + 0.00474222576341281;
			break;
		case 17:
			disp = -0.0016738259999153366*x + 0.004796903102053131;
			break;
		case 18:
			disp = -0.0017818424655581238*x + 0.004996733563492289;
			break;
		case 19:
			disp = -0.001754580578066302*x + 0.004943572882883236;
			break;
		case 20:
			disp = -0.0017511560409477574*x + 0.004936552581790219;
			break;
		case 21:
			disp = -0.0016615086580877136*x + 0.004743810708641125;
			break;
		case 22:
			disp = -0.0016377732342938544*x + 0.004690406005104944;
			break;
		case 23:
			disp = -0.0015589245588133042*x + 0.004505111617725649;
			break;
		case 24:
			disp = -0.0013622157432266867*x + 0.004023175019538436;
			break;
		case 25:
			disp = -0.0012848136633894635*x + 0.003825799715953518;
			break;
		case 26:
			disp = -0.001186792017676669*x + 0.00356604235481461;
			break;
		case 27:
			disp = -0.0010377304741342052*x + 0.003156123110072835;
			break;
		case 28:
			disp = -0.0009561981835930621*x + 0.0029237560820305772;
			break;
		case 29:
			disp = -0.0007725661085796497*x + 0.0023820414607410106;
			break;
		case 30:
			disp = -0.0006300967022673964*x + 0.0019475097714886371;
			break;
		case 31:
			disp = -0.0005075236717357815*x + 0.0015614047253140492;
			break;
		case 32:
			disp = -0.00035895113001251083*x + 0.0010785439647134198;
			break;
		case 33:
			disp = -0.0002668143223903539*x + 0.0007698856591791944;
			break;
		case 34:
			disp = -0.00020304862923710141*x + 0.0005498940178004738;
			break;
		case 35:
			disp = -0.00013503708291501433*x + 0.0003084530283570644;
			break;
		case 36:
			disp = -7.549209922880212e-05*x + 9.111383790238981e-05;
			break;
		case 37:
			disp = -3.413910160362498e-05*x + -6.395990319202459e-05;
			break;
		case 38:
			disp = 1.0573204549001701e-05*x + -0.00023610228187963722;
			break;
		case 39:
			disp = 4.514160919373606e-05*x + -0.00037264748022633806;
			break;
		case 40:
			disp = 7.293048664080089e-05*x + -0.0004851924338869505;
			break;
		case 41:
			disp = 8.463809354776643e-05*x + -0.0005337790025508574;
			break;
		case 42:
			disp = 0.00010984051664463937*x + -0.0006408893007125678;
			break;
		case 43:
			disp = 0.00011635088250418322*x + -0.0006692093922015836;
			break;
		case 44:
			disp = 0.00011811050681831445*x + -0.0006770397203994674;
			break;
		case 45:
			disp = 0.0001340012662205909*x + -0.0007493426756798252;
			break;
		case 46:
			disp = 0.00013691029225067338*x + -0.0007628696467197086;
			break;
		case 47:
			disp = 0.0001289896881501288*x + -0.0007252467772421224;
			break;
		case 48:
			disp = 0.00013565684245982898*x + -0.0007575824756441681;
			break;
		case 49:
			disp = 0.0001309892251631302*x + -0.000734477770025509;
			break;
		case 50:
			disp = 0.000124351877627342*x + -0.0007009591649697785;
			break;
		case 51:
			disp = 0.0001199478294999091*x + -0.0006782783171134986;
			break;
		case 52:
			disp = 0.00011835194367191442*x + -0.0006698999165165263;
			break;
		case 53:
			disp = 0.00011782349629181938*x + -0.000667072723033018;
			break;
		case 54:
			disp = 0.00011611741464083499*x + -0.0006577745780351532;
			break;
		default:
			disp = 0.0;
	}
	
	return disp;
}

Coordinate getFinalPositionAddArea(Coordinate pos, double A_exo, double time_term, double radius, double cylinderLength)
{
	// this will remain fixed throughtout this function
	double initialPhi = atan2(pos.y, pos.x);
	
	// fit to sphere simulations
	double sigma_exo = 0.848*2.0/3.0;
	double sigma_endo = 1.838*2.0/3.0;
	
	double alpha = 1.0;
	
	if(pos.z < 0.0)
	{
		//initial pos is on cylinder
		double s0 = pi*radius*0.5 - pos.z;
		
		
		double s_final = s0 + time_term*avgDisplacementAtArclength(s0);
		
		if(s_final < 0.0)
		{
			s_final = 0.0;
		}
		
		double final_z = -(s_final - pi*radius*0.5);
		
		Coordinate pos_3d = Coordinate {pos.x, pos.y, final_z};
		
		// now check if pos is below reflecting boundary
		if(pos_3d.z < -cylinderLength)
		{
			pos_3d.z = -2.0*cylinderLength - pos_3d.z;
		}
		
		return pos_3d;
	}
	else
	{
		// initial pos in on sphere
		double initialTheta = acos(pos.z / radius);
		double s0 = radius*initialTheta;
		
		
		
		double s_final = s0 + time_term*avgDisplacementAtArclength(s0);
		
		if(s_final < 0.0)
		{
			s_final = 0.0;
		}
		
		if(s_final > pi*radius*0.5)
		{
			// then the final point is on the cylinder
			double final_z = -(s_final - pi*radius*0.5);
		
			Coordinate pos_3d = Coordinate {pos.x, pos.y, final_z};
			
			// now check if pos is below reflecting boundary
			if(pos_3d.z < -cylinderLength)
			{
				pos_3d.z = -2.0*cylinderLength - pos_3d.z;
			}
			
			return pos_3d;
		}
		else
		{
			// then the final point is still on the sphere
			double finalTheta = s_final / radius;
			
			Coordinate pos_3d = {radius*cos(initialPhi)*sin(finalTheta), radius*sin(initialPhi)*sin(finalTheta), radius*cos(finalTheta)};
			return pos_3d;			
		}
			
	}
	
}




using namespace std;

int main(int argc, char** argv)
{
    // uses physics definition of theta and phi, polar angle is theta, azimuthal angle is phi
	clock_t t1,t2;
	t1=clock();
    
    double ratio = 0.001;
    
    double lambda_ST = 9.6;
    
    double RsRt_frac = 5.0;
    
    double lambda_STT = 5.0;
    
    //double D_D = 4E-4;
    
   // double r_D =  0.03; 
    
   // double r_T = 0.005;
    
    //double D_T = 4E-4;
    
    //double D_ST = D_T;
    
    double D_D = 0.1;
    double D_T = 0.1;
    double D_ST = 0.0025;
    
    
    double r_D =  0.17; 
    
    double r_T = 0.015;
    

    double k_ST = 0.015;
	
	double k_hydro = 0.35*0.5;
    
    
    std::string positionFileName = "";
        
    if(argc >= 2)
    {
        for(int i = 1; i < argc; i++)
        {
            switch(i)
            {
                case 1:
                    r_D = atof(argv[1]);
                    cout << "r_D: " << r_D << endl;
                    break;
                case 2:
                    D_D = atof(argv[2]);
                    cout << "D_D: " << D_D << endl;
                    break;
                case 3:
                    positionFileName = argv[3];
                    break;
				
                default:
                    cout << "Ignoring additional arguments beyond 3rd" << endl;
            }
        }
    }
    
    //double k_D = 100.0 / 3000.0;
    
   
    
    
    
    //cout << positionFileName << endl;
    
    
    if(r_D < r_T)
    {
        r_T = r_D;
    }
    /*else
    {
        r_T = 0.015;
    }*/
	
	//r_T = 0.0;
    
    
    if(D_D < D_T)
    {
        D_T = D_D;
    }
    /*else
    {
        D_T = 0.1;
    }*/
    
    if(D_T < D_ST)
    {
        D_ST = D_T;
    }
    /*else
    {
        D_ST = 0.0025;
    }*/
    
    
     // approximation of kD to keep N_gdp,internal constant
    double k_D = 1.0/30.0 * (1000*r_T+2000*r_D) / (1000*0.005+2000*0.03);
    std::cout << "r_T: " << r_T << ", k_D: " << k_D << std::endl;
    
    double D_S = 0.0025;
    
    //double k_D = 100.0 / 3000.0;
    double k_S = 10.0 * ratio;

    
    
    double D_GAP = 1E-3;
    double k_GAP = 100.0 / 3000.0;
    double r_GAP = 1E-3;
    r_GAP = 0.0;
   
    
    
    double r_S = 0.1*RsRt_frac;
    

    double exo_rate = 40.88014155/60.0;
    double endo_rate = 146.1504417/60.0;
    

    double snapshotTimestep = 1.0 / 1000.0;
	
	double initialRadius = 1.8;
	double cylinderLength = 3.0;
    
    
	double initialSphereArea = 4.0*pi*initialRadius*initialRadius;
	
	double initialSpherocylinderArea = initialSphereArea*0.5 + 2.0*pi*initialRadius*cylinderLength;
    
	//multiply all numbers of particles by this ratio
    double SpherocylSphereAreaRatio = initialSpherocylinderArea / initialSphereArea;
	
	
    double r_ST = 0.1*RsRt_frac;
    

    
    int totalNumberScd2Scd1 = (int)round(450*0.5 * SpherocylSphereAreaRatio);
	
	//https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
	std::random_device rd;
    std::mt19937 gen(time(0));
    //std::mt19937 gen(0);
    
   
	
	//double cylinderLengthX = 2*pi*initialRadius;
	//double cylinderLengthZ = 9.0-2.0*initialRadius;
    
    double b_val_height = 0.71;
    
    int totalNumberOfCdc42 = (int)round(20000*0.25*SpherocylSphereAreaRatio);
    
    
    double diffusion_prefactorGDP = sqrt(2.0 * D_D);
    double diffusion_prefactorGTP = sqrt(2.0 * D_T);
    double diffusion_prefactorScd1Scd2 = sqrt(2.0 * D_S);
    double diffusion_prefactorScd1Scd2GTP = sqrt(2.0 * D_ST);
    double diffusion_prefactorGAP = sqrt(2.0 * D_GAP);
    
    int maxNumTypes = 6+2;
    
    std::vector <double> diffusionPrefactors(maxNumTypes);
    diffusionPrefactors[0] = diffusion_prefactorGDP;
    diffusionPrefactors[1] = diffusion_prefactorGTP;
    diffusionPrefactors[2] = -10.0;
    diffusionPrefactors[3] = diffusion_prefactorScd1Scd2;
    diffusionPrefactors[4] = -10.0;
    diffusionPrefactors[5] = diffusion_prefactorScd1Scd2GTP;
    diffusionPrefactors[6] = diffusion_prefactorGAP;
    diffusionPrefactors[7] = -10.0;
    
    /*diffusionPrefactors[0] = 0.0;
    diffusionPrefactors[1] = 0.0;
    diffusionPrefactors[3] = 0.0;
    diffusionPrefactors[5] = 0.0;*/ 
    
    
    // set number of particles on surface equal to steady state value without exo/endocytosis changing the number of particles on the surface
    int membraneNumberOfCdc42 = (int)(totalNumberOfCdc42 * 0.5);
    int cytoNumberOfCdc42 = totalNumberOfCdc42 - membraneNumberOfCdc42;
    
    int membraneNumberOfScd2Scd1 = (int)(totalNumberScd2Scd1 * 0.5);
    int cytoNumberOfScd2Scd1 = totalNumberScd2Scd1 - membraneNumberOfScd2Scd1;
    
    int totalNumberGAP = (int)round(totalNumberScd2Scd1*10);
    
    int membraneNumberOfGAP = (int)(totalNumberGAP * 0.5);
    int cytoNumberOfGAP = totalNumberGAP - membraneNumberOfGAP;
    
    //cout << totalNumberGAP << " " << membraneNumberOfGAP << " " << cytoNumberOfGAP << endl;
    //exit(0);
    
    if(cytoNumberOfCdc42 < 0)
    {
        cout << "Internal Number of particles is: " << cytoNumberOfCdc42 << ", exiting" << endl;
        exit(0);
    }
    

    cout << membraneNumberOfCdc42 << " " << cytoNumberOfCdc42 << endl;
	
    std::uniform_real_distribution<> dis(0.0, 1.0);
    std::normal_distribution<> gaussianDis1(0.0, 1.0);
    
    
    
    std::vector <std::vector <struct Particle>> particles(maxNumTypes);
    
    
    bool readin_state = true;
    
    if(readin_state && positionFileName != "")
    {
        //std::string positionFileName = "3dParamSetD_rD_0.10000-frame4999.xyz";
        
        ifstream inputfile (positionFileName);

        
        if (inputfile.is_open())
        {
            int linecount = 0;
            
            int currFilamentType = -1;
            int newFilamentTag = -1;
            
            std::string line;
            
            int numMembraneParticles;
            
            std::vector <bool> usedIndices (totalNumberOfCdc42+totalNumberScd2Scd1+totalNumberGAP, false);
            
        
            while(std::getline(inputfile, line))
            {
                std::istringstream iss(line);
                
                if(linecount == 0)
                {
                    iss >> numMembraneParticles;
                }
                else if(linecount == 1)
                {
                    std::string tempString, secondString, thirdString;
                    iss >> tempString >> secondString >> thirdString;
                    
                    int pos = tempString.find("=");
					
					int pos_comma = tempString.find(",");
                    
                    std::string timeString = tempString.substr(pos+1,pos_comma);
					
					double currSimTime = std::stod(timeString);
                    
					std::cout << thirdString << std::endl;
					
					cylinderLength = std::stod(thirdString);
                    
                    //exit(0);
                }
                else if(linecount > 1)
                {
                    int tag, type;
                    double x,y,z;
                    
                    iss >> type >> tag >> x >> y >> z;
                    
                    Coordinate c = {x, y, z};
                    
                    usedIndices[tag] = true;
                    
                    particles[type].emplace_back(c, tag, 0.0, 0.0);
                }
                
                linecount++;
            }
                                
            cytoNumberOfCdc42 = totalNumberOfCdc42 - particles[0].size() - particles[1].size() - particles[5].size();
            cytoNumberOfScd2Scd1 = totalNumberScd2Scd1 - particles[3].size() - particles[5].size();
            cytoNumberOfGAP = totalNumberGAP - particles[6].size();

            for(int i = 0; i < cytoNumberOfCdc42; i++)
            {
                Coordinate c = {0.0, 0.0, 0.0};
                
                int currIndex = -1;
                for(int j = 0; j < totalNumberOfCdc42; j++)
                {
                    if(usedIndices[j] == false)
                    {
                        currIndex = j;
                        usedIndices[j] = true;
                        break;
                    }
                }
                
                if(currIndex == -1)
                {
                    cout << "not enough indexes, Cdc42" << endl;
                }
                
                particles[2].emplace_back(c, currIndex, 0.0, 0.0);
            }
            for(int i = 0; i < cytoNumberOfScd2Scd1; i++)
            {
                Coordinate c = {0.0, 0.0, 0.0};
                
                int currIndex = -1;
                for(int j = totalNumberOfCdc42; j < usedIndices.size(); j++)
                {
                    if(usedIndices[j] == false)
                    {
                        currIndex = j;
                        usedIndices[j] = true;
                        break;
                    }
                }
                
                if(currIndex == -1)
                {
                    cout << "not enough indexes, Scd2/1" << endl;
                }
                
                particles[4].emplace_back(c, currIndex, 0.0, 0.0);
            }
            
            for(int i = 0; i < particles[5].size(); i++)
            {
                int secondIndex = -1;
                for(int j = totalNumberOfCdc42; j < usedIndices.size(); j++)
                {
                    if(usedIndices[j] == false)
                    {
                        secondIndex = j;
                        usedIndices[j] = true;
                        break;
                    }
                }
                
                if(secondIndex == -1)
                {
                    cout << "not enough indexes, type 5: " << i << endl;
                }
                
                particles[5][i].secondname = secondIndex;
            }
            
            
            for(int i = 0; i < cytoNumberOfGAP; i++)
            {
                Coordinate c = {0.0, 0.0, 0.0};
                
                int currIndex = -1;
                for(int j = totalNumberOfCdc42; j < usedIndices.size(); j++)
                {
                    if(usedIndices[j] == false)
                    {
                        currIndex = j;
                        usedIndices[j] = true;
                        break;
                    }
                }
                
                if(currIndex == -1)
                {
                    cout << "not enough indexes, GAP" << endl;
                }
                
                particles[7].emplace_back(c, currIndex, 0.0, 0.0);
            }
            
            for(int i = 0; i < usedIndices.size(); i++)
            {
                if(usedIndices[i] == false)
                {
                    cout << "missing particle " << i << endl;
                }
            }
            //exit(0);
        }
    }	
	else
    {
        /*std::vector <unsigned int> names(membraneNumberOfCdc42);
        std::vector <int> types(membraneNumberOfCdc42);
        std::vector <double> timesLastGlobalUpdate(membraneNumberOfCdc42);
        std::vector <double> posTimes(membraneNumberOfCdc42);*/
        
        for(int i = 0; i < membraneNumberOfCdc42; i++)
        {
            //double phiRnd = 2.0*pi*dis(gen);
            //double thetaRnd = acos(dis(gen));
			
			double rndNumbers [3] = {dis(gen), dis(gen), dis(gen)};
            
            //Coordinate c = {initialRadius*cos(phiRnd)*sin(thetaRnd), initialRadius*sin(phiRnd)*sin(thetaRnd), initialRadius*cos(thetaRnd)};
			Coordinate c = randomPointSurface(initialRadius, cylinderLength, rndNumbers);
            
            particles[0].emplace_back(c, i, 0.0, 0.0);
        }
        
        for(int i = 0; i < cytoNumberOfCdc42; i++)
        {
            Coordinate c = {0.0, 0.0, 0.0};
            
            particles[2].emplace_back(c, i+membraneNumberOfCdc42, 0.0, 0.0);
        }
        
        for(int i = 0; i < membraneNumberOfScd2Scd1; i++)
        {
            //double phiRnd = 2.0*pi*dis(gen);
            //double thetaRnd = acos(dis(gen));
			
			double rndNumbers [3] = {dis(gen), dis(gen), dis(gen)};
            
            //Coordinate c = {initialRadius*cos(phiRnd)*sin(thetaRnd), initialRadius*sin(phiRnd)*sin(thetaRnd), initialRadius*cos(thetaRnd)};
			Coordinate c = randomPointSurface(initialRadius, cylinderLength, rndNumbers);
            
            particles[3].emplace_back(c, i+membraneNumberOfCdc42+cytoNumberOfCdc42, 0.0, 0.0);
        }
        
        for(int i = 0; i < cytoNumberOfScd2Scd1; i++)
        {
            Coordinate c = {0.0, 0.0, 0.0};
            
            particles[4].emplace_back(c, i+membraneNumberOfCdc42+cytoNumberOfCdc42+membraneNumberOfScd2Scd1, 0.0, 0.0);
        }
        
        for(int i = 0; i < membraneNumberOfGAP; i++)
        {
           //double phiRnd = 2.0*pi*dis(gen);
            //double thetaRnd = acos(dis(gen));
			
			double rndNumbers [3] = {dis(gen), dis(gen), dis(gen)};
            
            //Coordinate c = {initialRadius*cos(phiRnd)*sin(thetaRnd), initialRadius*sin(phiRnd)*sin(thetaRnd), initialRadius*cos(thetaRnd)};
			Coordinate c = randomPointSurface(initialRadius, cylinderLength, rndNumbers);
            
            //c = {0.0, 0.0, initialRadius};
            
            particles[6].emplace_back(c, i+membraneNumberOfCdc42+cytoNumberOfCdc42+membraneNumberOfScd2Scd1+cytoNumberOfScd2Scd1, 0.0, 0.0);
        }
        for(int i = 0; i < cytoNumberOfGAP; i++)
        {
            Coordinate c = {0.0, 0.0, 0.0};
            
            particles[7].emplace_back(c, i+membraneNumberOfCdc42+cytoNumberOfCdc42+membraneNumberOfScd2Scd1+cytoNumberOfScd2Scd1+membraneNumberOfGAP, 0.0, 0.0);
        }
    }
        
    /*std::queue <unsigned int> recycledNames;
    for(int j = 0; j < cytoNumberOfCdc42; j++)
    {
        recycledNames.push(j + membraneNumberOfCdc42);
    }
    
    unsigned int namesMax = membraneNumberOfCdc42 + cytoNumberOfCdc42;*/
	
	double currTime = 0.0;
	double finalTime = 5000.0;
    //double currPosTime = 0.0;
	
    std::ostringstream out;
    out.precision(5);
    out << std::fixed << r_D;
    
    std::ostringstream out2;
    out2.precision(5);
    out2 << std::fixed << D_D;
    
    std::ostringstream out3;
    out3.precision(5);
    out3 << std::fixed << k_hydro;
    
    std::ostringstream out4;
    out4.precision(5);
    out4 << std::fixed << r_T;
    
	//std::string appendString = "density_kon_" + out.str() + "_r_D_" + out2.str() + "_D_" + out3.str() + "_endo-rate_" + out4.str();
	std::string appendString = "_rD_" + out.str() + "_DD_" + out2.str() + "_Run1b";
    
	
	//double del_t = finalTime;

	double nextSnapshotTime = 0.0;
	
	std::ofstream xyzFile ("3d" + appendString + ".xyz");
    
    std::ofstream numEachStateFile ("NumEachState" + appendString + ".txt");
	
	int currFrame = 0;
	
	Coordinate exclusionCenter = {0.0, -initialRadius, 0.0};
	
	double phiRnd = 2.0*pi*dis(gen);
	double thetaRnd = acos(2*dis(gen)-1);
		
	exclusionCenter = {initialRadius*cos(phiRnd)*sin(thetaRnd), initialRadius*sin(phiRnd)*sin(thetaRnd), initialRadius*cos(thetaRnd)};
	
    int numRates = 5;
	std::vector <double> rates(numRates);

	std::vector <double> Rki(numRates);

	std::vector <int> countEvents(numRates);
    
    std::ofstream occupancyFile ("3doccupancyHist" + appendString + ".dat");
    std::ofstream exoPositionsFile ("exoPositions" + appendString + ".dat");
    std::ofstream endoPositionsFile ("endoPositions" + appendString + ".dat");
    
    
    double exo_radius = 50.0 / 1000.0;
    double endo_radius = 22.6 / 1000.0;
    
    double endo_exo_cutoff_multiplier = 2.0;
    
    double area_per_vesicle_exo = 4.0 * pi * exo_radius*exo_radius;
    double height_spherical_cap = area_per_vesicle_exo / 2.0 / pi / initialRadius;
    double a = sqrt(area_per_vesicle_exo / pi - height_spherical_cap * height_spherical_cap);
    double exo_surface_radius = initialRadius * atan(a / (initialRadius - height_spherical_cap));
    double exo_surface_cutoff = endo_exo_cutoff_multiplier * exo_surface_radius;
    
    double var_l_exo_compare = acos(area_per_vesicle_exo / (2*pi*initialRadius*initialRadius) - 1.0)*initialRadius;
    var_l_exo_compare = 2.0;
    
    cout << exo_surface_cutoff << endl;
    
    double area_per_vesicle_endo = 4.0 * pi * endo_radius*endo_radius;
    height_spherical_cap = area_per_vesicle_endo / 2.0 / pi / initialRadius;
    a = sqrt(area_per_vesicle_endo / pi - height_spherical_cap * height_spherical_cap);
    double endo_surface_radius = initialRadius * atan(a / (initialRadius - height_spherical_cap));
    double endo_surface_cutoff = endo_exo_cutoff_multiplier * endo_surface_radius;
    
    double var_l_endo_compare = acos(area_per_vesicle_endo / (2*pi*initialRadius*initialRadius) - 1.0)*initialRadius;
    var_l_endo_compare = 2.0;
    
    cout << endo_surface_cutoff << endl;

    double kapa = 1.0; 
    
	// encapsulate net exo/endo area addition rate into a single parameter (um^2/s)
	double AreaAddRate = (exo_rate*area_per_vesicle_exo - endo_rate*area_per_vesicle_endo)*1.0;
	std::cout << "AreaAddRate: " << AreaAddRate << std::endl;
	
    //exo_rate = 0.0;
    //endo_rate = 0.0;
	
    rates[0] = 0.0;
    rates[1] = 0.0;
    
    
    double alpha = 1.0;
    
    int frameModValue = (int)(1.0 / snapshotTimestep);
    
    //currTime is the current system time at all times, currPosTime is the time at which the current position of the particles is at

	while(currTime < finalTime)
	{	
		//rates[1] = 254.5227 / 60.0;
        
        // assumes that all internal particles are GDP, not GTP
		rates[2] = k_D * particles[2].size();
        rates[3] = k_S * particles[4].size();
        rates[4] = k_GAP * particles[7].size();
        
		for(int i = 0; i < numRates; i++)
		{
			if(i == 0)
				Rki[i] = rates[i];
			else
				Rki[i] = rates[i] + Rki[i-1];
		}
		
        // time until next reaction
		double del_t = -1.0 / Rki[numRates-1] * log(dis(gen));
		
		//std::cout << currTime << " " << del_t << std::endl;
		
		//del_t = 1.0;
		//std::cout << del_t << std::endl;
		
		// print configuration
		while(currTime + del_t >= nextSnapshotTime)
		{
			
			
            int tempOccupancyTip = 0;
            int tempOccupancyBack = 0;
    
			if(currFrame % frameModValue == 0)
            {
                int totNumParticles = 0;
                for(int i = 0; i < particles.size(); i++)
                {
                    if(diffusionPrefactors[i] >= 0.0)
                    {
                        totNumParticles += particles[i].size();
                    }
                }
                        
                xyzFile << totNumParticles+1 << endl;
                xyzFile << "t=" << currTime << ", cl: " << cylinderLength << endl;
                
                numEachStateFile << currTime << ",";
                
                for(int i = 0; i < particles.size(); i++)
                {
                    if(diffusionPrefactors[i] >= 0.0)
                    {
                        for(int p = 0; p < particles[i].size(); p++)
                        {                            
                            xyzFile << i << " " << particles[i][p].name << " " << particles[i][p].pos.x << " " << particles[i][p].pos.y << " " << particles[i][p].pos.z << endl;
                        }
                    }
                    
                    numEachStateFile << particles[i].size() << ",";
                }
                
                xyzFile << 10 << " " << -1 << " " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
                
                numEachStateFile << endl;

                //xyzFile << 1 << " " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
                //xyzFile << 2 << " " << 0.0 << " " << 0.0 << " " << initialRadius << endl;
                //xyzFile << 3 << " " << 0.0 << " " << 0.0 << " " << -initialRadius << endl;
            }
			
			currFrame++;
			
			cylinderLength = AreaAddRate / alpha *(nextSnapshotTime - currTime)/2.0/pi/initialRadius + cylinderLength;
            // update to nextSnapshotTime from currPosTime
            for(int i = 0; i < particles.size(); i++)
            {
                double curr_diffusionPrefactor = diffusionPrefactors[i];
                
                if(curr_diffusionPrefactor > 0.0)
                {
                    for(int p = particles[i].size()-1; p > -1; p--)
                    {
                        double time_term = sqrt(nextSnapshotTime - particles[i][p].posTime);
                        
                        double gaussDirectionX = gaussianDis1(gen);
                        double gaussDirectionY = gaussianDis1(gen);

						double delPosX = curr_diffusionPrefactor * time_term * gaussDirectionX;
						double delPosY = curr_diffusionPrefactor * time_term * gaussDirectionY;
						
						//std::cout << "time: " << time_term << " " << particles[i][p].posTime << " " <<  nextSnapshotTime << " " << del_t << std::endl;
						
						if(time_term > 0.0)
						{
							//Coordinate initPos = particles[i][p].pos;
							
							particles[i][p].pos = diffuseParticle(particles[i][p].pos, delPosX, delPosY, initialRadius, cylinderLength);
							
							//Coordinate intermediatePos = particles[i][p].pos;
							
							
							particles[i][p].pos = getFinalPositionAddArea(particles[i][p].pos, area_per_vesicle_exo, time_term*time_term, initialRadius, cylinderLength);
							
							//Coordinate finalPos = particles[i][p].pos;
							
							/*if(i == 6)
							{
								// this is a GAP particle
								std::cout << calcArcLengthDistance(intermediatePos, initPos, initialRadius) << " " << calcArcLengthDistance(finalPos, intermediatePos, initialRadius) << " " << time_term*time_term << std::endl;
							}*/
							
							if(particles[i][p].pos.z >= 0.0)
							{
								tempOccupancyTip += 1;
							}
							else
							{
								tempOccupancyBack += 1;
							}
							
							particles[i][p].posTime = nextSnapshotTime;
						}
                    }
                }
            }
            
            int countScd1Reaction = 0;
            int countScd1InCyto = particles[4].size();
            int countCdc42GTP = particles[1].size();
            
            for(int i = 0; i < particles.size(); i++)
            {
                //cout << i << endl;
                for(int p = particles[i].size()-1; p > -1; p--)
                {

                    if(i == 0)
                    {
                        // conversion between GDP and GTP via Scd1Scd2
                        Coordinate pos_i = particles[i][p].pos;
                        
                        bool convertedToGTP = false;
                        
                        for(int k = 0; k < particles[3].size(); k++)
                        {
                            double tempArcLength = calcArcLengthDistance(pos_i, particles[3][k].pos, initialRadius);
                            
                            if(tempArcLength < 0.05)
                            {
                                //lambda parameter for Doi method
                                double lambda_SD = 5.3;
                                //double probability = lambda*(nextSnapshotTime - particles[i][p].timeLastGlobalUpdate);
                                double probability = 1.0 - exp(-lambda_SD * (nextSnapshotTime - particles[i][p].timeLastGlobalUpdate));
                                
                                double randVal = dis(gen);
                                
                                if(randVal < probability)
                                {
                                    // switch this particle from GDP to GTP
                                    Particle tempParticle = particles[0][p];
                                    particles[0][p] = particles[0].back();
                                    particles[0].pop_back();
                                    
                                    particles[1].push_back(tempParticle);
                                    convertedToGTP = true;
                                    break;
                                }
                            }
                        }
                        
                        if(convertedToGTP == false)
                        {
                            //conversion from GDP to GTP due to scd1scd2gtp complex
                            for(int k = 0; k < particles[5].size(); k++)
                            {
                                double tempArcLength = calcArcLengthDistance(pos_i, particles[5][k].pos, initialRadius);
                                
                                if(tempArcLength < 0.05)
                                {
                                    double probability = 1.0 - exp(-lambda_STT * (nextSnapshotTime - particles[i][p].timeLastGlobalUpdate));
                                    
                                    double randVal = dis(gen);
                                    
                                    if(randVal < probability)
                                    {
                                        // switch this particle from GDP to GTP
                                        Particle tempParticle = particles[0][p];
                                        particles[0][p] = particles[0].back();
                                        particles[0].pop_back();
                                        
                                        particles[1].push_back(tempParticle);
                                        
                                        break;
                                    }
                                }
                            }
                        }
                        
                    }
                    else if(i == 1)
                    {
                        bool convertedToGDP = false;
                        // conversion between GTP and GDP                        
                        double currProb = 1.0 - exp(-k_hydro * (nextSnapshotTime - particles[i][p].timeLastGlobalUpdate));
                        
                        double randomUniform = dis(gen);
                        
                        if(randomUniform < currProb)
                        {
                            // switch this particle from GTP to GDP
                            Particle tempParticle = particles[1][p];
                            particles[1][p] = particles[1].back();
                            particles[1].pop_back();
                            
                            particles[0].push_back(tempParticle);
                            
                            convertedToGDP = true;
                        }
                        
                        // hydrolosis by GAP
                        if(convertedToGDP == false)
                        {
                            Coordinate pos_i = particles[i][p].pos;
                            
                            for(int k = 0; k < particles[6].size(); k++)
                            {
                                double tempArcLength = calcArcLengthDistance(pos_i, particles[6][k].pos, initialRadius);
                                
                                if(tempArcLength < 0.05)
                                {
                                    //lambda parameter for Doi method
                                    double lambda_GAP = 0.1*100.0;
                                    //double probability = lambda*(nextSnapshotTime - particles[i][p].timeLastGlobalUpdate);
                                    double probability = 1.0 - exp(-lambda_GAP * (nextSnapshotTime - particles[i][p].timeLastGlobalUpdate));
                                    
                                    double randVal = dis(gen);
                                    
                                    if(randVal < probability)
                                    {
                                        Particle tempParticle = particles[1][p];
                                        particles[1][p] = particles[1].back();
                                        particles[1].pop_back();
                                        
                                        particles[0].push_back(tempParticle);
                                        
                                        convertedToGDP = true;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    else if(i == 3)
                    {
                        // conversion from Scd2/Scd1 and GTP to Scd2/Scd1/Cdc42GTP
                        Coordinate pos_i = particles[i][p].pos;
                        
                        for(int k = particles[1].size()-1; k > -1; k--)
                        {
                            double tempArcLength = calcArcLengthDistance(pos_i, particles[1][k].pos, initialRadius);
                            
                            if(tempArcLength < 0.05)
                            {
                                //lambda parameter for Doi method and Ramierez
                                //double lambda_ST = 9.6;
                                //double probability = lambda*(nextSnapshotTime - particles[i][p].timeLastGlobalUpdate);
                                double probability = 1.0 - exp(-lambda_ST * (nextSnapshotTime - particles[i][p].timeLastGlobalUpdate));
                                
                                double randVal = dis(gen);
                                
                                if(randVal < probability)
                                {
									// NEED TO CHECK FOR SPHEROCYLINDER!!!
                                    // switch Scd1Scd2 and GTP to Scd1Scd2GTP compound
                                    Particle tempParticle = particles[1][k];
                                    tempParticle.secondname = particles[3][p].name;
                                    
                                    Coordinate final_pos3D = 0.5*(particles[3][p].pos - particles[1][k].pos) + particles[1][k].pos;
									
									if(final_pos3D.z > 0.0)
									{
										// final position is on sphere surface
										final_pos3D = final_pos3D.getUnitCoord();
										final_pos3D = initialRadius * final_pos3D;
									}
									else
									{
										// final position is on cylinder surface
										double initialPhi = atan2(final_pos3D.y, final_pos3D.x);
										
										final_pos3D = Coordinate{initialRadius*cos(initialPhi), initialRadius*sin(initialPhi), final_pos3D.z};										
									}
						
									/*std::cout << particles[3][p].pos.x << " " << particles[3][p].pos.y << " " << particles[3][p].pos.z << std::endl;
                                    std::cout << particles[1][k].pos.x << " " << particles[1][k].pos.y << " " << particles[1][k].pos.z << std::endl;
									std::cout << final_pos3D.x << " " << final_pos3D.y << " " << final_pos3D.z << std::endl;
                                    
									std::cout << calcArcLengthDistance(particles[3][p].pos, final_pos3D, initialRadius) << " " << calcArcLengthDistance(particles[1][k].pos, final_pos3D, initialRadius) << std::endl;*/
                                    particles[3][p] = particles[3].back();
                                    particles[3].pop_back();
                                    
                                    
                                    particles[1][k] = particles[1].back();
                                    particles[1].pop_back();
                                    
                                    
                                    tempParticle.pos = final_pos3D;
                            
                                    particles[5].push_back(tempParticle);
                                    
                                    break;
                                }
                            }
                        }
                    }
                    else if(i == 4)
                    {
                        // conversion from Scd2/Scd1(cytoplasm) and GTP to Scd2/Scd1/Cdc42GTP
                        Coordinate pos_i = particles[i][p].pos;
                        
                        for(int k = particles[1].size()-1; k > -1; k--)
                        {
                            //double tempArcLength = calcArcLengthDistance(pos_i, particles[1][k].pos, initialRadius);
                            
                            //if(tempArcLength < 0.05)
                            {
                                //lambda parameter for Doi method and Ramierez
                                //double lambda = 256.0;
                                //double probability = lambda*(nextSnapshotTime - particles[i][p].timeLastGlobalUpdate);
                                //double probability = 1.0 - exp(-lambda * (nextSnapshotTime - particles[i][p].timeLastGlobalUpdate));
                                
                                double randVal = dis(gen);
                                
                                double currProb = 1.0 - exp(-k_ST * (nextSnapshotTime - particles[i][p].timeLastGlobalUpdate));
                            
                                if(randVal < currProb)
                                {
                                    countScd1Reaction++;
                                    // switch Scd1Scd2 and GTP to Scd1Scd2GTP compound
                                    
                                    Particle tempParticle = particles[1][k];
                                    tempParticle.secondname = particles[4][p].name;
                                    
                                    particles[4][p] = particles[4].back();
                                    particles[4].pop_back();
                                    
                                    particles[1][k] = particles[1].back();
                                    particles[1].pop_back();
                                    
                                    particles[5].push_back(tempParticle);
                                    
                                    break;
                                }
                            }
                        }
                    }
                    else if(i == 5)
                    {
                        // conversion from Scd2/Scd1/Cdc42GTP to Scd2/Scd1 and GTP
                        double currProb = 1.0 - exp(-r_ST * (nextSnapshotTime - particles[i][p].timeLastGlobalUpdate));
                        
                        double randomUniform = dis(gen);
                        
                        if(randomUniform < currProb)
                        {
                            Particle initialParticle = particles[i][p];
                            particles[i][p] = particles[i].back();
                            particles[i].pop_back();
                            
							double phi = 2.0*pi*dis(gen);
							double unbinding_rad = 0.055;
							
							
							double rndX = unbinding_rad*0.5*cos(phi);
							double rndY = unbinding_rad*0.5*sin(phi);
							
							Coordinate final_pos3d = diffuseParticle(initialParticle.pos, rndX, rndY, initialRadius, cylinderLength);
							particles[3].emplace_back(final_pos3d, initialParticle.secondname, initialParticle.timeLastGlobalUpdate, initialParticle.posTime);
							
							rndX = -rndX;
							rndY = -rndY;
							Coordinate final_pos3d_second = diffuseParticle(initialParticle.pos, rndX, rndY, initialRadius, cylinderLength);
							particles[1].emplace_back(final_pos3d_second, initialParticle.name, initialParticle.timeLastGlobalUpdate, initialParticle.posTime);
							
							
							/*std::cout << final_pos3d.x << " " << final_pos3d.y << " " << final_pos3d.z << std::endl;
							std::cout << final_pos3d_second.x << " " << final_pos3d_second.y << " " << final_pos3d_second.z << std::endl;
							
							std::cout << calcArcLengthDistance(final_pos3d_second, final_pos3d, initialRadius) << std::endl << std::endl;*/
                        }
                    }
                }
            }
            
            //cout << countScd1InCyto << " " << countScd1Reaction << " " << countCdc42GTP << endl;
            
            for(int i = 0; i < particles.size(); i++)
            {
                //cout << "two " << i << endl;
                for(int p = particles[i].size()-1; p > -1; p--)
                {
                    bool removedParticle = false;
                    // put koff simulations here instead since self-consistancy issues if in the continuous time part
                    
                    if(i == 0)
                    {
                        double currProb = 1.0 - exp(-r_D * (nextSnapshotTime - particles[i][p].timeLastGlobalUpdate));
                        
                        double randomUniform = dis(gen);
                        
                        if(randomUniform < currProb)
                        {
                            Particle tempParticle = particles[i][p];
                            particles[i][p] = particles[i].back();
                            particles[i].pop_back();
                            
                            particles[2].push_back(tempParticle);
                                    
                            cytoNumberOfCdc42 += 1;
                            //countEvents[3] += 1;
                            removedParticle = true;
                        }
                    }
                    else if(i == 1)
                    {
                        double currProb = 1.0 - exp(-r_T * (nextSnapshotTime - particles[i][p].timeLastGlobalUpdate));
                        
                        double randomUniform = dis(gen);
                        
                        if(randomUniform < currProb)
                        {
                            Particle tempParticle = particles[i][p];
                            particles[i][p] = particles[i].back();
                            particles[i].pop_back();
                            
                            particles[2].push_back(tempParticle);
                                    
                            cytoNumberOfCdc42 += 1;
                            //countEvents[3] += 1;
                            removedParticle = true;
                        }
                    }
                    else if(i == 3)
                    {
                        // koff for Scd2/Scd1
                        double currProb = 1.0 - exp(-r_S * (nextSnapshotTime - particles[i][p].timeLastGlobalUpdate));
                        
                        double randomUniform = dis(gen);
                        
                        if(randomUniform < currProb)
                        {
                            Particle tempParticle = particles[i][p];
                            particles[i][p] = particles[i].back();
                            particles[i].pop_back();
                            
                            particles[4].push_back(tempParticle);
                                    
                            cytoNumberOfCdc42 += 1;
                            //countEvents[3] += 1;
                            removedParticle = true;
                        }
                    }
                    else if(i == 6)
                    {
                        // koff for GAP
                        double currProb = 1.0 - exp(-r_GAP * (nextSnapshotTime - particles[i][p].timeLastGlobalUpdate));
                        
                        double randomUniform = dis(gen);
                        
                        if(randomUniform < currProb)
                        {
                            Particle tempParticle = particles[i][p];
                            particles[i][p] = particles[i].back();
                            particles[i].pop_back();
                            
                            particles[7].push_back(tempParticle);
                                    
                            //cytoNumberOfCdc42 += 1;
                            //countEvents[3] += 1;
                            removedParticle = true;
                        }
                    }
                    
                    if(removedParticle == false)
                    {
                        particles[i][p].timeLastGlobalUpdate = nextSnapshotTime;
                    }
                }
            }
                
            
            
            // update del_t to account for time particles diffused in this global update
            del_t = del_t - (nextSnapshotTime - currTime);
            
            currTime = nextSnapshotTime;
            nextSnapshotTime += snapshotTimestep;
            
            if(tempOccupancyBack > 0)
            {
                // ok to just divide the counts rather than densities since the area is constant       
                occupancyFile << currTime << " " << (double)tempOccupancyTip / (double)tempOccupancyBack << endl;
            }
		}
		
        // find out which event to select
		double randU = dis(gen);
		double randUQk = randU * Rki[numRates-1];
		
		int eventIndex = 0;
		
		if(randUQk > Rki[0])
		{
			for(int i = 1; i < numRates; i++)
			{
				if(randUQk > Rki[i-1] && randUQk < Rki[i])
					eventIndex = i;
			}
		}
        
        
        //only really need to update all particles if the event type is exo or endocytosis
		//update position of particles (diffusion)
		cylinderLength = AreaAddRate / alpha *(del_t)/2.0/pi/initialRadius + cylinderLength;
        if(eventIndex <= 1)
        {
            // update system to currTime + del_t from currPosTime
            for(int i = 0; i < particles.size(); i++)
            {
                double curr_diffusionPrefactor = diffusionPrefactors[i];
                
                if(curr_diffusionPrefactor > 0.0)
                {
                    for(int p = particles[i].size()-1; p > -1; p--)
                    {
                        double time_term = sqrt(del_t + currTime - particles[i][p].posTime);
                        
                        double gaussDirectionX = gaussianDis1(gen);
                        double gaussDirectionY = gaussianDis1(gen);
						
						double delPosX = curr_diffusionPrefactor * time_term * gaussDirectionX;
						double delPosY = curr_diffusionPrefactor * time_term * gaussDirectionY;
	
                        particles[i][p].pos = diffuseParticle(particles[i][p].pos, delPosX, delPosY, initialRadius, cylinderLength);
						particles[i][p].pos = getFinalPositionAddArea(particles[i][p].pos, area_per_vesicle_exo, time_term*time_term, initialRadius, cylinderLength);
    
                        //particles[i][p].pos = diffuseParticleSphere(time_term, curr_diffusionPrefactor, particles[i][p].pos, gaussDirectionX, gaussDirectionY, initialRadius);
                        
                        particles[i][p].posTime = del_t + currTime;
                    }
                }
            }
        }
		
		//THIS ALL NEEDS CHANGING SPHEROCYLINDER!!!
        /*if(eventIndex == 0)
		{
			//exocytosis
			double exclusionRadius = exo_surface_radius;
            
            double radiusRatio = exo_radius / initialRadius;            
			
            double phiRnd, thetaRnd;
            Coordinate c;
            double randomUniform = 1.0;
            
            double sigma_distribution = 0.01;
            double inv_sigma_distribution = 1.0 / sigma_distribution;
            double inv_sigma_distribution_sq = inv_sigma_distribution*inv_sigma_distribution;
            
            //cout << "exocytosis attempt" << endl;
            //int attemptCount = 0;
            
            if(particles[1].size() > 0.0)
            {
                double prob_value = 0.0;
                do
                {
                    // get random position					
					double rndNumbers[] = [dis(gen), dis(gen), dis(gen)];
					c = randomPointSurface(radius, cylinderLength, rndNumbers);
					
                    
                    prob_value = 0.0;
                    // determine probability to add particle at position c by measuring arc length distances (Cdc42-GTP)
                    for(int p = 0; p < particles[1].size(); p++)
                    {
                        Coordinate temp_pos = particles[1][p].pos;
                        
                        double dist = calcArcLengthDistance(c, temp_pos, initialRadius);
                        
                        prob_value += exp(-0.5*dist*dist*inv_sigma_distribution_sq);
                    }
                    
                    prob_value = prob_value / (double)particles[1].size();
                        
                    
                    randomUniform = dis(gen);
                    
                    //cout << "exo attempt count: " << attemptCount << " " << prob_value << endl;
                    //attemptCount++;
                    
                    if(prob_value >= 1.0)
                        cout << prob_value << " " << randomUniform << endl;
                }
                while(randomUniform > prob_value);
            }
            else
            {
				double rndNumbers[] = [dis(gen), dis(gen), dis(gen)];
				c = randomPointSurface(radius, cylinderLength, rndNumbers);	
            }
            
            
			exclusionCenter = c;
            
            // need to find max Cdc42-GTP concentration point and rotate exclusionCenter based around this point (north pole -> max Cdc42-GTP point)
            
            exoPositionsFile << exclusionCenter.x << " " << exclusionCenter.y << " " << exclusionCenter.z << endl;
			
			Coordinate v1 = exclusionCenter.getUnitCoord();
			
			//push particles out of exclusion zone
            for(int i = 0; i < particles.size(); i++)
            {
                if(diffusionPrefactors[i] >= 0.0)
                {
                    for(int p = particles[i].size()-1; p > -1; p--)
                    {                    
                        if(exclusionRadius > pi*initialRadius)
                            exclusionRadius = pi*initialRadius;
                        
                        struct Coordinate disp = particles[i][p].pos - exclusionCenter;
                        
                        // straight line distance between pos[p] and exclusion center
                        double var_d = sqrt(disp.x*disp.x + disp.y*disp.y + disp.z*disp.z);
                        
                        if(var_d > 0.0)
                        {
                        
                            // arc length distance between pos[p] and exclusion center using law of cosines
                            double var_l = initialRadius * acos(1.0 - var_d*var_d / 2.0 / (initialRadius*initialRadius));
                            
                            if(var_l < var_l_exo_compare)
                            {
                                Coordinate v2 = particles[i][p].pos - exclusionCenter;
                                
                                // projection of v2 onto refrence plane with normal vector v1
                                v2 = v2 - (v2*v1)*v1;
                                
                                v2 = v2.getUnitCoord();
                                
                                // angle made between exclusion center and final position
                                // kapa is how steep the dropoff of the push is (should go to zero if var_l = exo_surface_cutoff
                                // temp_exclusion_radius is additional distance away from the exocytotic center that the particle should move
                                //double var_t = (exclusionRadius) / initialRadius;
                                //double var_t = (exclusionRadius + var_l) / initialRadius;
                                
                                
                                double RSphereNew = sqrt(initialRadius*initialRadius + area_per_vesicle_exo / (4.0 * pi));
                                
                                double acosParameter = 1.0 - initialRadius*initialRadius / (RSphereNew*RSphereNew) * (1.0 - cos(var_l / initialRadius)) - area_per_vesicle_exo / (2.0 * pi * RSphereNew*RSphereNew * alpha);
                                if(acosParameter < -1.0)
                                {
                                   // cout << "acosParameter < 0.0: " << acosParameter << endl;
                                    acosParameter = -1.0;
                                }
                                else if(acosParameter > 1.0)
                                {
                                    //cout << "acosParameter > pi: " << acosParameter << endl;
                                    acosParameter = 1.0;
                                }
                                
                                double delR = RSphereNew * acos(acosParameter) - var_l;
                                double var_t = (var_l + delR) / initialRadius;
                                
                                if(var_t > pi)
                                    var_t = pi;
                                
                                Coordinate finalPos = initialRadius*cos(var_t)*v1 + initialRadius*sin(var_t)*v2;
                                
                                //cout << var_l << " " << var_l+delR << endl;
                                particles[i][p].pos = finalPos;
                            }
                        }
                    }
                }
            }
            
            
			countEvents[0] += 1;
		}*/
		/*else if(eventIndex == 1)
		{
			//endocytosis
			double inclusionRadius = endo_surface_radius;
			
            
            
            double phiRnd, thetaRnd;
            Coordinate c;
            double randomUniform = 1.0;
            
            double sigma_distribution = 1.51;
            double inv_sigma_distribution = 1.0 / sigma_distribution;
            double inv_sigma_distribution_sq = inv_sigma_distribution*inv_sigma_distribution;
            
            if(particles[1].size() > 0.0)
            {
                double prob_value = 0.0;
                do
                {
                    // get random position					
					double rndNumbers[] = [dis(gen), dis(gen), dis(gen)];
					c = randomPointSurface(radius, cylinderLength, rndNumbers);
			
                    
                    prob_value = 0.0;
                    // determine probability to add particle at position c by measuring arc length distances (Cdc42-GTP)
                    for(int p = 0; p < particles[1].size(); p++)
                    {
                        Coordinate temp_pos = particles[1][p].pos;
                        
                        double dist = calcArcLengthDistance(c, temp_pos, initialRadius);
                        
                        prob_value += exp(-0.5*dist*dist*inv_sigma_distribution_sq);
                    }
                    
                    prob_value = prob_value / (double)particles[1].size();
                        
                    
                    randomUniform = dis(gen);
                    
                    if(prob_value >= 1.0)
                        cout << prob_value << " " << randomUniform << endl;
                }
                while(randomUniform > prob_value);
            }
            else
            {
				double rndNumbers[] = [dis(gen), dis(gen), dis(gen)];
				c = randomPointSurface(radius, cylinderLength, rndNumbers);
            }
            
            
			Coordinate inclusionCenter = c;
            
      
            
            endoPositionsFile << inclusionCenter.x << " " << inclusionCenter.y << " " << inclusionCenter.z << endl;
			
			Coordinate v1 = inclusionCenter.getUnitCoord();
			
			int proteinRemoveCount = 0;
			
			//pull particles into inclusion zone
            for(int i = 0; i < particles.size(); i++)
            {
                if(diffusionPrefactors[i] >= 0.0)
                {
                    for(int p = particles[i].size()-1; p > -1; p--)
                    {
                        if(inclusionRadius > pi*initialRadius)
                            inclusionRadius = pi*initialRadius;
                        
                        struct Coordinate disp = particles[i][p].pos - inclusionCenter;
                        
                        double var_d = sqrt(disp.x*disp.x + disp.y*disp.y + disp.z*disp.z);
                        
                        if(var_d > 0.0)
                        {
                        
                            double var_l = initialRadius * acos(1.0 - var_d*var_d / 2.0 / (initialRadius*initialRadius));
                            
                            if(var_l < inclusionRadius)
                            {
                                particles[i][p].pos = inclusionCenter;
                            }
                            else if(var_l < var_l_endo_compare)
                            //else
                            {
                                Coordinate v2 = particles[i][p].pos - inclusionCenter;
                                
                                v2 = v2 - (v2*v1)*v1;
                                
                                v2 = v2.getUnitCoord();
                                
                                // angle made between exclusion center and final position
                                //double var_t = (var_l - inclusionRadius) / initialRadius;
                                
                                double RSphereNew = sqrt(initialRadius*initialRadius - area_per_vesicle_endo / (4.0 * pi));
                                
                                double acosParameter = 1.0 - initialRadius*initialRadius / (RSphereNew*RSphereNew) * (1.0 - cos(var_l / initialRadius)) + area_per_vesicle_endo / (2.0 * pi * RSphereNew*RSphereNew * alpha);
                                if(acosParameter < -1.0)
                                {
                                   // cout << "acosParameter < 0.0: " << acosParameter << endl;
                                    acosParameter = -1.0;
                                }
                                else if(acosParameter > 1.0)
                                {
                                    //cout << "acosParameter > pi: " << acosParameter << endl;
                                    acosParameter = 1.0;
                                }
                                
                                double delR = var_l - RSphereNew * acos(acosParameter);
                                double var_t = (var_l - delR) / initialRadius;
                                
                                if(var_t > pi)
                                    var_t = pi;
                                
                                Coordinate finalPos = initialRadius*cos(var_t)*v1 + initialRadius*sin(var_t)*v2;
                                particles[i][p].pos = finalPos;
                            }		
                        }
                    }
                }
            }
			
			cytoNumberOfCdc42 += proteinRemoveCount;
			
			countEvents[1] += 1;
		}*/
		
		if(eventIndex == 2)
		{
            //cout << "Cdc42 " << particles[2].size()  << endl;
			//Cdc42 proteinAdd
			/*double phiRnd = 2.0*pi*dis(gen);
			double thetaRnd = acos(2*dis(gen)-1);
			
			Coordinate c = {initialRadius*cos(phiRnd)*sin(thetaRnd), initialRadius*sin(phiRnd)*sin(thetaRnd), initialRadius*cos(thetaRnd)};*/
			double rndNumbers [3] = {dis(gen), dis(gen), dis(gen)};
			Coordinate c = randomPointSurface(initialRadius, cylinderLength, rndNumbers);

            if(particles[2].size() >= 1)
            {
                int randomIndex = (int)(particles[2].size()*dis(gen));
                
                Particle tempParticle = particles[2][randomIndex];
                particles[2][randomIndex] = particles[2].back();
                particles[2].pop_back();
                
                tempParticle.pos = c;
                tempParticle.timeLastGlobalUpdate = currTime + del_t;
                tempParticle.posTime = currTime + del_t;
                
                particles[0].push_back(tempParticle);
                
                cytoNumberOfCdc42 -= 1;
                
                countEvents[2] += 1;
            }
		}
        else if(eventIndex == 3)
        {
            //cout << "Scd1Scd2 " << particles[4].size()  << endl;
            //Scd1Scd2 proteinAdd
			/*double phiRnd = 2.0*pi*dis(gen);
			double thetaRnd = acos(2*dis(gen)-1);
			
			Coordinate c = {initialRadius*cos(phiRnd)*sin(thetaRnd), initialRadius*sin(phiRnd)*sin(thetaRnd), initialRadius*cos(thetaRnd)};*/
			double rndNumbers [3] = {dis(gen), dis(gen), dis(gen)};
			Coordinate c = randomPointSurface(initialRadius, cylinderLength, rndNumbers);
			
            //c = {0.0, 0.0, initialRadius};
            
            if(particles[4].size() >= 1)
            {
                int randomIndex = (int)(particles[4].size()*dis(gen));
                
                
                Particle tempParticle = particles[4][randomIndex];
                particles[4][randomIndex] = particles[4].back();
                particles[4].pop_back();
                
                tempParticle.pos = c;
                tempParticle.timeLastGlobalUpdate = currTime + del_t;
                tempParticle.posTime = currTime + del_t;
                
                particles[3].push_back(tempParticle);
                
                cytoNumberOfCdc42 -= 1;
                
                countEvents[3] += 1;
            }
        }
        else if(eventIndex == 4)
        {
            //cout << "GAP add: " << particles[7].size() << " " << particles[6].size() << endl;
            //GAP proteinAdd
			/*double phiRnd = 2.0*pi*dis(gen);
			double thetaRnd = acos(2*dis(gen)-1);

			Coordinate c = {initialRadius*cos(phiRnd)*sin(thetaRnd), initialRadius*sin(phiRnd)*sin(thetaRnd), initialRadius*cos(thetaRnd)};*/
			
			double rndNumbers [3] = {dis(gen), dis(gen), dis(gen)};
			Coordinate c = randomPointSurface(initialRadius, cylinderLength, rndNumbers);
            //c = {0.0, 0.0, initialRadius};
            
            if(particles[7].size() >= 1)
            {
                int randomIndex = (int)(particles[7].size()*dis(gen));
                
                
                Particle tempParticle = particles[7][randomIndex];
                particles[7][randomIndex] = particles[7].back();
                particles[7].pop_back();
                
                tempParticle.pos = c;
                tempParticle.timeLastGlobalUpdate = currTime + del_t;
                tempParticle.posTime = currTime + del_t;
                
                particles[6].push_back(tempParticle);
                
                cytoNumberOfCdc42 -= 1;
                
                countEvents[4] += 1;
            }
        }
        
		/*else if(eventIndex == 3)
		{
			//proteinRemove
			int randInt = (int)(particles.size() * dis(gen));
			//randMembraneProtein = random.randint(0, len(particles)-1);
			int randMembraneProtein = randInt;
			
			particles[randMembraneProtein] = particles.back();
			particles.pop_back();
            
            recycledNames.push(names[randMembraneProtein]);
            
            names[randMembraneProtein] = names.back();
            names.pop_back();
            
            types[randMembraneProtein] = types.back();
            types.pop_back();
					
			cytoNumberOfCdc42 += 1;
			countEvents[3] += 1;
		}*/
        
        currTime = currTime + del_t;
	}
	
	t2=clock();
	float diff ((float)t2-(float)t1);
	float seconds = diff / CLOCKS_PER_SEC;
	
    xyzFile.close();
    occupancyFile.close();
    exoPositionsFile.close();
    endoPositionsFile.close();
    numEachStateFile.close();
    
	//cout<<"Total runtime (s): " << seconds<<endl;
	
	/*for(int i = 0; i < countEvents.size(); i++)
	{
		cout << i << " " << countEvents[i] << endl;
	}*/
}
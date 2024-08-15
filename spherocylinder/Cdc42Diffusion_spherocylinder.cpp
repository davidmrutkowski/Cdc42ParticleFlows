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

/* 
   This code diffuses (diffuseParticle), reacts (lines 1336-1662; 1692-1769),
   and exerts the average effects of exo/endocytosis as a function of arclength 
   from the half-spherocylinder tip (getFinalPositionAddArea, avgDisplacementAtArclength). 
   Values in avgDisplacementAtArclength function are determined from sphere simulations at
   either WT or 1ritC conditions. The python code in this folder 
   (FlowsFromNPPatch_functionalFit.py) pulls in an xyz file to determine these values.
   on structs of type Particle restricted to the surface of a half-spherocylinder.
   The spherocylinder length (cylinderLength) increases as the simulation proceeds in
   order to reflect the growing cell suface.   
   The simulations proceeds in discrete timesteps of dt or (for Particle addition to
   the surface) del_t. The simulation runs until finalTime taking xyz snapshots every
   snapshotTime. 
   If given as a command line argument, positionFileName (an xyz frame)
   can be used to restart the simulation from a previous state rather  than initializing
   the Particles to random positions on the half-spherocylinder surface.
   
   Particle types refer to the following species:
   0: Cdc42-GDP
   1: Cdc42-GTP
   2: Cdc42-GDP (cytoplasmic)
   3: Scd2/1
   4: Scd2/1 (cytoplasmic)
   5: Scd2/1/Cdc42-GTP
   6: sGAP
   7: sGAP (cytoplasmic)
   
   This code was compiled with the C++17 standard using the GNU compiler.
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

// struct used to store 3D positions and vectors
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

//struct used to store information about each particle on the surface
struct Particle
{
    Coordinate pos;
    int name;
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

// calculates straight line distance between posi and posj
double calcDistance(Coordinate posi, Coordinate posj)
{
    Coordinate rij = posj - posi;
    
    return rij.getMagnitude();
}

// calculates arclength distance for two particles on a sphere surface
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

// calculates arclength distance for two positions on the cylinder surface
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

// calculates the arclength between two particles on the half-spherocylinder by projection
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
    
// function definitions for use by diffuseParticle
Coordinate diffuseParticleSphere(Coordinate pos, double delPosX, double delPosY, double radius, double cylinderLength);
Coordinate diffuseParticleCylinder(Coordinate pos, double delPosX, double delPosY, double radius, double cylinderLength);
Coordinate diffuseParticleSphere_fromCyl(Coordinate pos, Coordinate tangent_3d, double length, double radius);

// general diffuse particle for both surfaces, calls appropriate function based on where particle starts;
// if the particle z position is greater than 0.0 then it is on the hemisphere otherwise it is on the cylinder
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

// diffuse a particle on the hemisphere surface starting at pos; can call diffuseParticleCylinder if particle crosses
// the boundary between the hemisphere and then cylinder
Coordinate diffuseParticleSphere(Coordinate pos, double delPosX, double delPosY, double radius, double cylinderLength)
{
    double initialTheta = acos(pos.z / radius);
    double initialPhi = atan2(pos.y, pos.x);
    
    double stepsize = sqrt(delPosX*delPosX + delPosY*delPosY);              
    
    // hankbesser on github
    double theta = stepsize / radius;
    double phi = atan2(delPosY, delPosX);
    
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
    
    final_pos3D = final_pos3D*cos(initialTheta) + currCrossProduct*sin(initialTheta) + rotK*currDotProduct*(1.0 - cos(initialTheta));
    
    rotK.x = 0.0;
    rotK.y = 0.0;
    rotK.z = 1.0;
    
    // rotate by initialPhi
    currDotProduct = rotK*final_pos3D;
    currCrossProduct = rotK.crossProduct(final_pos3D);
    
    final_pos3D = final_pos3D*cos(initialPhi) + currCrossProduct*sin(initialPhi) + rotK*currDotProduct*(1.0 - cos(initialPhi));
    
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
    double d_acos = acos(d);
    
    if(tm > 1E-10 && tm < d_acos)
    {
        // particle crossed the boundary to the cylinder
        double arclength = tm * radius;
        
        double remainingLengthCylinder = stepsize - arclength;
        
        Coordinate tangentVectorAtBoundary = -sin(tm)*r0 + cos(tm)*w1;
        
        Coordinate startingPosCylinder = radius*(cos(tm)*r0 + sin(tm)*w1);

        double dirX = sqrt(tangentVectorAtBoundary.x*tangentVectorAtBoundary.x + tangentVectorAtBoundary.y*tangentVectorAtBoundary.y);
        double dirZ = tangentVectorAtBoundary.z;
        
        double cyl_delX = remainingLengthCylinder*dirX;
        double cyl_delY = remainingLengthCylinder*dirZ;
        
        // now diffuse on cylinder starting at startingPosCylinder, going in direction tangentVectorAtBoundary for length remainingLengthCylinder
        final_pos3D = diffuseParticleCylinder(startingPosCylinder, cyl_delX, cyl_delY, radius, cylinderLength);
        
    }
    
    return final_pos3D;
}

// diffuse a particle on the cylinder surface starting at pos; can call diffuseParticleSphere if particle crosses
// the boundary between the cylinder and then hemisphere
Coordinate diffuseParticleCylinder(Coordinate pos, double delPosX, double delPosY, double radius, double cylinderLength)
{   
    double initialX = (atan2(pos.y, pos.x) + pi) / (2.0*pi) * (2.0*pi*radius);
    double initialZ = pos.z;
    
    double finalX = initialX + delPosX;
    double finalZ = initialZ + delPosY;
    
    while(finalX > 2*pi*radius)
    {
        finalX = finalX - 2*pi*radius;
    }
    
    double finalAngle = finalX / (2.0*pi*radius) * (2.0*pi);    
    
    Coordinate pos_3d = Coordinate {-radius * cos(finalAngle), -radius * sin(finalAngle), finalZ};
    
    // now check if pos is below reflecting boundary
    if(pos_3d.z < -cylinderLength)
    {
        pos_3d.z = -2.0*cylinderLength - pos_3d.z;
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
        
        //now diffuse this on circle
        Coordinate pos_3d = diffuseParticleSphere(boundary_pos_3d, realStepX, realStepY, radius, cylinderLength);
        
    }
    
    return pos_3d;
    
}

// generate a random point on the surface of the half-spherocylinder
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

// returns a random point on a cylinder surface (with top edge at z=0)
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

// determines final arclength starting at arclength s0 and adding area AreaAdded
// to a half-spherocylinder with radius radius
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

// values in this function are taken from average displacement measurements from Sphere code
// for WT condition (run python code FlowsFromNPPatch_functionalFit.py to get values)
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

// gets final position for a particle starting at pos due to net effect of exo/endocytosis as
// determiend by avgDisplacementAtArclength for the starting arclength of the pos
Coordinate getFinalPositionAddArea(Coordinate pos, double time_term, double radius, double cylinderLength)
{
    // this will remain fixed throughtout this function
    double initialPhi = atan2(pos.y, pos.x);
    
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
    double pi = 3.14159265359;
    double initialRadius = 1.8;
    double cylinderLength = 3.0;
    double initialSphereArea = 4.0*pi*initialRadius*initialRadius;
    double initialSpherocylinderArea = initialSphereArea*0.5 + 2.0*pi*initialRadius*cylinderLength;
    
    //multiply all numbers of particles by this ratio
    double SpherocylSphereAreaRatio = initialSpherocylinderArea / initialSphereArea;    
    
    //Table S2 Parameters
    double dt = 1.0 / 1000.0;
    double D_T = 0.1;
    double D_D = 0.1;
    double D_S = 0.0025;
    double D_ST = 0.0025;
    double D_GAP = 1E-3;

    double r_T = 0.015;
    double r_D =  0.17;
    double r_S = 0.5;
    double r_ST = 0.5;
    double r_GAP = 0.0;
    double k_hydro = 0.175;
    double k_S = 0.01;
    double k_ST = 0.015;
    // rho is hardcoded in reactions below
    // sigma, radius at which to place to particles that unbind
    double unbinding_rad = 0.055;
    
    //lambda parameters for Doi method
    double lambda_SD = 5.3;
    double lambda_ST = 9.6;
    double lambda_STT = 5.0;
    double lambda_GAP = 10.0;
    
    int totalNumberOfCdc42 = (int)round(20000*0.25*SpherocylSphereAreaRatio);
    int totalNumberScd2Scd1 = (int)round(450*0.5 * SpherocylSphereAreaRatio);
    int totalNumberGAP = (int)round(totalNumberScd2Scd1*10);
    
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
    
    if(r_D < r_T)
    {
        r_T = r_D;
    }
    if(D_D < D_T)
    {
        D_T = D_D;
    }
    if(D_T < D_ST)
    {
        D_ST = D_T;
    }
    
    // approximation of kD to keep N_gdp,internal constant
    double k_D = 1.0/30.0 * (1000*r_T+2000*r_D) / (1000*0.005+2000*0.03);
    std::cout << "r_T: " << r_T << ", k_D: " << k_D << std::endl;
    
    // on rate of sGAP
    double k_GAP = 100.0 / 3000.0;
                              

    //Table S3 parameters
    double exo_radius = 50.0 / 1000.0;
    double area_per_vesicle_exo = 4.0 * pi * exo_radius*exo_radius;
    
    double endo_radius = 22.6 / 1000.0;
    double area_per_vesicle_endo = 4.0 * pi * endo_radius*endo_radius;

    double exo_rate = 40.88014155/60.0;
    double endo_rate = 146.1504417/60.0;
    double w_exo = 0.01;
    double w_endo = 1.51;
    double alpha = 0.5;
    //gamma implicitly 1.0
    // R_cutoff parameter
    double var_l_exo_compare = 2.0;
    double var_l_endo_compare = 2.0;
    
    
    // time related variables
    double currTime = 0.0;
    double finalTime = 5000.0;
    
    double nextSnapshotTime = 0.0;
    
    //time at which to add another snapshot to xyz output file
    double snapshotTime = 10.0;
    int frameModValue = (int)(snapshotTime / dt);
    
    //https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
    std::random_device rd;
    std::mt19937 gen(time(0));

    double diffusion_prefactorGDP = sqrt(2.0 * D_D);
    double diffusion_prefactorGTP = sqrt(2.0 * D_T);
    double diffusion_prefactorScd1Scd2 = sqrt(2.0 * D_S);
    double diffusion_prefactorScd1Scd2GTP = sqrt(2.0 * D_ST);
    double diffusion_prefactorGAP = sqrt(2.0 * D_GAP);
    
    int maxNumTypes = 8;
    
    std::vector <double> diffusionPrefactors(maxNumTypes);
    diffusionPrefactors[0] = diffusion_prefactorGDP;
    diffusionPrefactors[1] = diffusion_prefactorGTP;
    diffusionPrefactors[2] = -10.0;
    diffusionPrefactors[3] = diffusion_prefactorScd1Scd2;
    diffusionPrefactors[4] = -10.0;
    diffusionPrefactors[5] = diffusion_prefactorScd1Scd2GTP;
    diffusionPrefactors[6] = diffusion_prefactorGAP;
    diffusionPrefactors[7] = -10.0;
    
    
    // set number of particles on surface equal to steady state value without exo/endocytosis changing the number of particles on the surface
    int membraneNumberOfCdc42 = (int)(totalNumberOfCdc42 * 0.5);
    int cytoNumberOfCdc42 = totalNumberOfCdc42 - membraneNumberOfCdc42;
    
    int membraneNumberOfScd2Scd1 = (int)(totalNumberScd2Scd1 * 0.5);
    int cytoNumberOfScd2Scd1 = totalNumberScd2Scd1 - membraneNumberOfScd2Scd1;
    
    int membraneNumberOfGAP = (int)(totalNumberGAP * 0.5);
    int cytoNumberOfGAP = totalNumberGAP - membraneNumberOfGAP;
    
    if(cytoNumberOfCdc42 < 0)
    {
        cout << "Internal Number of particles is: " << cytoNumberOfCdc42 << ", exiting" << endl;
        exit(0);
    }
    

    cout << membraneNumberOfCdc42 << " " << cytoNumberOfCdc42 << endl;
    
    std::uniform_real_distribution<> dis(0.0, 1.0);
    std::normal_distribution<> gaussianDis1(0.0, 1.0);
    
    std::vector <std::vector <struct Particle>> particles(maxNumTypes);
    
    // initializing Particle positions either from a previous xyz frame (positionFileName)
    // or randomly placing them on the sphere surface                            
    if(positionFileName != "")
    {
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
                    
                    //read in the cylinderLength from the xyz file
                    cylinderLength = std::stod(thirdString);
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
        }
    }   
    else
    {   
        // random initialization of the particles on the surface of the spherocylinder
        // sGAP is initialized on the cylinder part of the spherocylinder only
        for(int i = 0; i < membraneNumberOfCdc42; i++)
        {           
            double rndNumbers [3] = {dis(gen), dis(gen), dis(gen)};
            
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
            double rndNumbers [3] = {dis(gen), dis(gen), dis(gen)};
            
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
            double rndNumbers [3] = {dis(gen), dis(gen), dis(gen)};
            
            Coordinate c = randomPointCylinder(initialRadius, cylinderLength, rndNumbers);

            particles[6].emplace_back(c, i+membraneNumberOfCdc42+cytoNumberOfCdc42+membraneNumberOfScd2Scd1+cytoNumberOfScd2Scd1, 0.0, 0.0);
        }
        for(int i = 0; i < cytoNumberOfGAP; i++)
        {
            Coordinate c = {0.0, 0.0, 0.0};
            
            particles[7].emplace_back(c, i+membraneNumberOfCdc42+cytoNumberOfCdc42+membraneNumberOfScd2Scd1+cytoNumberOfScd2Scd1+membraneNumberOfGAP, 0.0, 0.0);
        }
    }
    
    std::ostringstream out;
    out.precision(5);
    out << std::fixed << r_D;
    
    std::ostringstream out2;
    out2.precision(5);
    out2 << std::fixed << D_D;
    
    
    std::string appendString = "_rD_" + out.str() + "_DD_" + out2.str() + "_Run1";
    
    std::ofstream xyzFile ("3d" + appendString + ".xyz");
    
    std::ofstream numEachStateFile ("NumEachState" + appendString + ".txt");
    
    int currFrame = 0;

    
    int numRates = 5;
    std::vector <double> rates(numRates);

    std::vector <double> Rki(numRates);

    std::vector <int> countEvents(numRates);
    
    std::ofstream occupancyFile ("3doccupancyHist" + appendString + ".dat");
    std::ofstream exoPositionsFile ("exoPositions" + appendString + ".dat");
    std::ofstream endoPositionsFile ("endoPositions" + appendString + ".dat");
    

    // encapsulate net exo/endo area addition rate into a single parameter (um^2/s)
    double AreaAddRate = (exo_rate*area_per_vesicle_exo - endo_rate*area_per_vesicle_endo)*1.0;
    std::cout << "AreaAddRate: " << AreaAddRate << std::endl;
    
    // Gillespie algorithm parameters for exo/endocytosis, both are zero since
    // for spherocylinder simulations discrete events are not implemented,
    // instead avgDisplacementAtArclength is used to implement the net effects
    // of exo/endocytosis
    rates[0] = 0.0;
    rates[1] = 0.0;
    
    //main simulation loop    
    //currTime is the current system time at all times, currPosTime is the time at which the current position of the particles is at
    while(currTime < finalTime)
    {   
        // Gillespie rates than depend on particle amount on surface
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
        
        // advance to next snapshot time while currTime + next Gillspie event time is not yet reached
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
            }
            
            currFrame++;
            
            //update cylinder length due to additional area added to surface between nextSnapshotTime and currTime
            cylinderLength = AreaAddRate *(nextSnapshotTime - currTime)/2.0/pi/initialRadius + cylinderLength;
            
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
                        
                        if(time_term > 0.0)
                        {
                            particles[i][p].pos = diffuseParticle(particles[i][p].pos, delPosX, delPosY, initialRadius, cylinderLength);
                            
                            particles[i][p].pos = getFinalPositionAddArea(particles[i][p].pos, time_term*time_term, initialRadius, cylinderLength);
                            
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
            
            // go through reactions based on discrete timestep
            for(int i = 0; i < particles.size(); i++)
            {
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
                                double probability = 1.0 - exp(-lambda_ST * (nextSnapshotTime - particles[i][p].timeLastGlobalUpdate));
                                
                                double randVal = dis(gen);
                                
                                if(randVal < probability)
                                {
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
                        }
                    }
                }
            }
            
            //additional koff unbinding reactions based on discrete timestep
            for(int i = 0; i < particles.size(); i++)
            {
                for(int p = particles[i].size()-1; p > -1; p--)
                {
                    bool removedParticle = false;
                    // put koff simulations here instead since self-consistancy issues if in the continuous time part
                    
                    if(i == 0)
                    {
                        // koff for Cdc42-GDP
                        double currProb = 1.0 - exp(-r_D * (nextSnapshotTime - particles[i][p].timeLastGlobalUpdate));
                        
                        double randomUniform = dis(gen);
                        
                        if(randomUniform < currProb)
                        {
                            Particle tempParticle = particles[i][p];
                            particles[i][p] = particles[i].back();
                            particles[i].pop_back();
                            
                            particles[2].push_back(tempParticle);
                                    
                            cytoNumberOfCdc42 += 1;
                            removedParticle = true;
                        }
                    }
                    else if(i == 1)
                    {
                        // koff for Cdc42-GTP
                        double currProb = 1.0 - exp(-r_T * (nextSnapshotTime - particles[i][p].timeLastGlobalUpdate));
                        
                        double randomUniform = dis(gen);
                        
                        if(randomUniform < currProb)
                        {
                            Particle tempParticle = particles[i][p];
                            particles[i][p] = particles[i].back();
                            particles[i].pop_back();
                            
                            particles[2].push_back(tempParticle);
                                    
                            cytoNumberOfCdc42 += 1;
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
            nextSnapshotTime += dt;
            
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
        
        if(eventIndex == 2)
		{
			//Cdc42 proteinAdd
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
            //Scd1Scd2 proteinAdd

			double rndNumbers [3] = {dis(gen), dis(gen), dis(gen)};
			Coordinate c = randomPointSurface(initialRadius, cylinderLength, rndNumbers);
            
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
            //GAP proteinAdd
			
			double rndNumbers [3] = {dis(gen), dis(gen), dis(gen)};
			Coordinate c = randomPointCylinder(initialRadius, cylinderLength, rndNumbers);
            
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

        currTime = currTime + del_t;
    }
    
    xyzFile.close();
    occupancyFile.close();
    exoPositionsFile.close();
    endoPositionsFile.close();
    numEachStateFile.close();
}
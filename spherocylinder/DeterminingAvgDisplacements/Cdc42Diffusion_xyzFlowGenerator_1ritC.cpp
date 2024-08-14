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
 
/* Modified version of Cdc42Diffusion code that generates averaged displacement vectors
   associated with discrete points on the surface of the sphere as determined by ipack.3.33002.txt
   (given at http://neilsloane.com/icosahedral.codes/ ) for a given configuration of Cdc42-GTP as
   defined in "Spherocylinder simulation" section of the manuscript. No particles diffuse or move 
   in this code; only displacement vectors are calculated and averaged.
   Type 1 - Cdc42-GTP
   Type 6 - Tracer particles
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
    
    Coordinate totalDisplacement;
    
    int secondname;
    
    Particle(Coordinate c, int n, double t, double pt, Coordinate td)
        : pos(c)
        , name(n)
        , timeLastGlobalUpdate(t)
        , posTime(pt)
        , totalDisplacement(td)
        {}
};

// periodic wrap for positions varying between 0 and length
double periodicWrap(double pos, double length)
{
    pos = pos - length * 0.5;
    return pos - length * round(pos / length) + length*0.5;
}
    
double calcDistance(double posi, double posj, double length)
{
    double dist = posi - posj;
    return dist;
}

double calcDistance(Coordinate posi, Coordinate posj)
{
    Coordinate rij = posj - posi;
    
    return rij.getMagnitude();
}

double calcArcLengthDistance(Coordinate posi, Coordinate posj, double sphereRadius)
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

Coordinate diffuseParticle(double time_term, double prefactor, Coordinate pos, double gaussDirectionX, double gaussDirectionY, double radius)
{
    double initialTheta = acos(pos.z / radius);
    double initialPhi = atan2(pos.y, pos.x);
 
    
    double delPosX, delPosY;
    
    
    delPosX = prefactor * time_term * gaussDirectionX;
    delPosY = prefactor * time_term * gaussDirectionY;
    
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
    
    return final_pos3D;
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
    
    
    
    
    
    double r_ST = 0.1*RsRt_frac;
    

    
    int totalNumberScd2Scd1 = (int)(450*0.5);
    
    //https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
    std::random_device rd;
    std::mt19937 gen(time(0));
    //std::mt19937 gen(0);
    
    double pi = 3.14159265359;
    double initialRadius = 1.8;
    double cylinderLengthX = 2*pi*initialRadius;
    double cylinderLengthZ = 9.0-2.0*initialRadius;
    
    double b_val_height = 0.71;
    
    int totalNumberOfCdc42 = (int)(20000*0.25);
    
    
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
    
    int totalNumberGAP = totalNumberScd2Scd1*10;
    
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
    std::normal_distribution<> gaussianDis2(0.0, 1.13170366);
    
    
    
    std::vector <std::vector <struct Particle>> particles(maxNumTypes);
    
    
    bool readin_state = true;
    
    
    {
        /*std::vector <unsigned int> names(membraneNumberOfCdc42);
        std::vector <int> types(membraneNumberOfCdc42);
        std::vector <double> timesLastGlobalUpdate(membraneNumberOfCdc42);
        std::vector <double> posTimes(membraneNumberOfCdc42);*/
        
        // define two centers around which to place Gaussian distributed particles
        Coordinate center_one = Coordinate {initialRadius, 0.0, 0.0};
        Coordinate center_two = Coordinate {0.0, initialRadius, 0.0};
        
        Coordinate center_three = initialRadius*(center_one + center_two).getUnitCoord();
        
        Coordinate center_north = Coordinate {0.0, 0.0, initialRadius};
        
        double sigma = 1.0;
            
        for(int i = 0; i < membraneNumberOfCdc42; i++)
        {
            double phiRnd = 2.0*pi*dis(gen);
            double thetaRnd = acos(2*dis(gen)-1);
            
            Coordinate c = {initialRadius*cos(phiRnd)*sin(thetaRnd), initialRadius*sin(phiRnd)*sin(thetaRnd), initialRadius*cos(thetaRnd)};
            
            double gaussDirectionX = gaussianDis2(gen);
            double gaussDirectionY = gaussianDis2(gen);
            
            Coordinate final_pos = diffuseParticle(1.0, 1.0, center_north, gaussDirectionX, gaussDirectionY, initialRadius);
            
            particles[1].emplace_back(final_pos, i, 0.0, 0.0, Coordinate{0.0, 0.0, 0.0});
        }
        
        /*for(int i = 0; i < membraneNumberOfCdc42; i++)
        {
            double phiRnd = 2.0*pi*dis(gen);
            double thetaRnd = acos(2*dis(gen)-1);
            
            Coordinate c = {initialRadius*cos(phiRnd)*sin(thetaRnd), initialRadius*sin(phiRnd)*sin(thetaRnd), initialRadius*cos(thetaRnd)};
            
            double gaussDirectionX = gaussianDis1(gen);
            double gaussDirectionY = gaussianDis1(gen);
            
            Coordinate final_pos = diffuseParticle(1.0, 0.3, center_three, gaussDirectionX, gaussDirectionY, initialRadius);
            
            particles[1].emplace_back(final_pos, i, 0.0, 0.0, Coordinate{0.0, 0.0, 0.0});
        }*/
        
       
        //std::string positionFileName = "3dParamSetD_rD_0.10000-frame4999.xyz";
        
        ifstream closepackfile ("ipack.3.33002.txt");

        
        if (closepackfile.is_open())
        {
            std::string line;
            
            int linecount = 0;

            double x,y,z;
            
            while(std::getline(closepackfile, line))
            {
                std::istringstream iss(line);
                
                
                
                if(linecount % 3 == 0)
                {
                    iss >> x;
                }
                else if(linecount % 3 == 1)
                {
                    iss >> y;
                }
                else if(linecount % 3 == 2)
                {
                    iss >> z;
                    
                    particles[6].emplace_back(Coordinate {x*initialRadius, y*initialRadius, z*initialRadius}, membraneNumberOfCdc42*2.0 + particles[6].size(), 0.0, 0.0, Coordinate{0.0, 0.0, 0.0});
                }
                
                linecount += 1;
            }
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
    std::string appendString = "_rD_" + out.str() + "_DD_" + out2.str() + "_singleNPPatch_correctedWidth";
    
    
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
    
    //exo_rate = 0.0;
    //endo_rate = 0.0;
    
    rates[0] = exo_rate;
    rates[1] = endo_rate;
    
    
    double alpha = 0.5;
    
    int frameModValue = (int)(10.0 / snapshotTimestep);
    
    //currTime is the current system time at all times, currPosTime is the time at which the current position of the particles is at

    int eventCount = 0;
    int numAttempts = 10000;
    
    while(eventCount < numAttempts)
    {   
        //int eventIndex = 0;
        
        // find out which event to select
        double randU = dis(gen);
        double randUQk = randU * Rki[numRates-1];
        
        int eventIndex = 0;
        
        for(int i = 0; i < numRates; i++)
        {
            if(i == 0)
                Rki[i] = rates[i];
            else
                Rki[i] = rates[i] + Rki[i-1];
        }
        
        
        if(randUQk > Rki[0])
        {
            for(int i = 1; i < numRates; i++)
            {
                if(randUQk > Rki[i-1] && randUQk < Rki[i])
                    eventIndex = i;
            }
        }
        
        std::cout << eventCount << " " << eventIndex << std::endl;
       
        if(eventIndex == 0)
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
                    phiRnd = 2.0*pi*dis(gen);
                    thetaRnd = acos(2*dis(gen)-1);
                    
                    c = {initialRadius*cos(phiRnd)*sin(thetaRnd), initialRadius*sin(phiRnd)*sin(thetaRnd), initialRadius*cos(thetaRnd)};
                    
                    prob_value = 0.0;
                    // determine probability to add particle at position c by measuring arc length distances (Cdc42-GTP)
                    for(int p = 0; p < particles[1].size(); p++)
                    {
                        Coordinate temp_pos = particles[1][p].pos;
                        
                        double dist = calcArcLengthDistance(c, temp_pos, initialRadius);
                        
                        prob_value += exp(-0.5*dist*dist*inv_sigma_distribution_sq);
                    }
                    
                    prob_value = prob_value / ((double)particles[1].size()*0.01);
                        
                    
                    randomUniform = dis(gen);
                    
                    //cout << "exo attempt count: " << attemptCount << " " << prob_value << endl;
                    //attemptCount++;
                    
                    if(prob_value >= 1.0)
                    {
                        std::cout << "Prob. value too large exocytosis" << std::endl;
                        cout << prob_value << " " << randomUniform << endl;
                        exit(0);
                    }
                }
                while(randomUniform > prob_value);
            }
            else
            {
                phiRnd = 2.0*pi*dis(gen);
                thetaRnd = acos(2*dis(gen)-1);
                
                c = {initialRadius*cos(phiRnd)*sin(thetaRnd), initialRadius*sin(phiRnd)*sin(thetaRnd), initialRadius*cos(thetaRnd)};
            }
            
            
            exclusionCenter = c;
            
            
            exoPositionsFile << exclusionCenter.x << " " << exclusionCenter.y << " " << exclusionCenter.z << endl;
            
            Coordinate v1 = exclusionCenter.getUnitCoord();
            
            // record displacement vectors on gridppoitnts (GAP) of sphere
            for(int i = 6; i < 7; i++)
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
                            /*particles[p] = particles.back();
                            particles.pop_back();
                            
                            recycledNames.push(names[p]);
                            
                            names[p] = names.back();
                            names.pop_back();
                            
                            cytoNumberOfCdc42 += 1;*/

                            Coordinate v2 = particles[i][p].pos - exclusionCenter;
                            
                            // projection of v2 onto refrence plane with normal vector v1
                            v2 = v2 - (v2*v1)*v1;
                            
                            v2 = v2.getUnitCoord();
                            
                            // angle made between exclusion center and final position
                            // kapa is how steep the dropoff of the push is (should go to zero if var_l = exo_surface_cutoff
                            // temp_exclusion_radius is additional distance away from the exocytotic center that the particle should move
                            //double var_t = (exclusionRadius) / initialRadius;
                            //double var_t = (exclusionRadius + var_l) / initialRadius;
                            /*double gamma = (1.0 / (var_l + kapa) - 1.0 / (exo_surface_cutoff + kapa)) / (1.0/kapa - 1.0 / (exo_surface_cutoff + kapa));
                            double temp_exclusion_radius = exclusionRadius * gamma;
                            double var_t = (temp_exclusion_radius + var_l) / initialRadius;*/
                            
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
                            
                            particles[i][p].totalDisplacement = particles[i][p].totalDisplacement + (finalPos - particles[i][p].pos);
                        }
                    }
                }
            }
            
            
            countEvents[0] += 1;
        }
        else if(eventIndex == 1)
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
                    phiRnd = 2.0*pi*dis(gen);
                    thetaRnd = acos(2*dis(gen)-1);
                    
                    c = {initialRadius*cos(phiRnd)*sin(thetaRnd), initialRadius*sin(phiRnd)*sin(thetaRnd), initialRadius*cos(thetaRnd)};
                    
                    prob_value = 0.0;
                    // determine probability to add particle at position c by measuring arc length distances (Cdc42-GTP)
                    for(int p = 0; p < particles[1].size(); p++)
                    {
                        Coordinate temp_pos = particles[1][p].pos;
                        
                        double dist = calcArcLengthDistance(c, temp_pos, initialRadius);
                        
                        prob_value += exp(-0.5*dist*dist*inv_sigma_distribution_sq);
                    }
                    
                    prob_value = prob_value / (double)(particles[1].size());
                        
                    
                    randomUniform = dis(gen);
                    
                    if(prob_value >= 1.0)
                        cout << prob_value << " " << randomUniform << endl;
                }
                while(randomUniform > prob_value);
            }
            else
            {
                phiRnd = 2.0*pi*dis(gen);
                thetaRnd = acos(2*dis(gen)-1);
                
                c = {initialRadius*cos(phiRnd)*sin(thetaRnd), initialRadius*sin(phiRnd)*sin(thetaRnd), initialRadius*cos(thetaRnd)};
            }
            
            
            Coordinate inclusionCenter = c;
            
       
            
            
            
            
            endoPositionsFile << inclusionCenter.x << " " << inclusionCenter.y << " " << inclusionCenter.z << endl;
            
            Coordinate v1 = inclusionCenter.getUnitCoord();
            
            int proteinRemoveCount = 0;
            
            //pull particles into inclusion zone
            for(int i = 6; i < 7; i++)
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
                            particles[i][p].totalDisplacement = particles[i][p].totalDisplacement + (inclusionCenter - particles[i][p].pos);
                            //particles[i][p].pos = inclusionCenter;
                            /*particles[p] = particles.back();
                            particles.pop_back();
                            
                            recycledNames.push(names[p]);
                            
                            names[p] = names.back();
                            names.pop_back();
                            
                            proteinRemoveCount++;*/
                        }
                        else if(var_l < var_l_endo_compare)
                        //else
                        {
                            Coordinate v2 = particles[i][p].pos - inclusionCenter;
                            
                            v2 = v2 - (v2*v1)*v1;
                            
                            v2 = v2.getUnitCoord();
                            
                            // angle made between exclusion center and final position
                            //double var_t = (var_l - inclusionRadius) / initialRadius;
                           /* double gamma = (1.0 / (var_l - inclusionRadius + kapa) - 1.0 / (endo_surface_cutoff - inclusionRadius + kapa)) / (1.0/kapa - 1.0 / (endo_surface_cutoff - inclusionRadius + kapa));
                            double temp_inclusion_radius = inclusionRadius * gamma;
                            double var_t = (var_l - temp_inclusion_radius) / initialRadius;*/
                            
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
                            
                            particles[i][p].totalDisplacement = particles[i][p].totalDisplacement + (finalPos - particles[i][p].pos);
                        }       
                    }
                }
            }
            
            countEvents[1] += 1;
        }
        
        eventCount += 1;
    }

    xyzFile << particles[6].size() + particles[1].size() << endl << endl;
    
    //now write xyz file
    
    
    for(int p = 0; p < particles[1].size(); p++)
     {
         xyzFile << 1 << " " << p << " " << particles[1][p].pos.x << " " <<  particles[1][p].pos.y << " " << particles[1][p].pos.z << " " << particles[1][p].totalDisplacement.x / (double)numAttempts << " " << particles[1][p].totalDisplacement.y / (double)numAttempts << " " << particles[1][p].totalDisplacement.z / (double)numAttempts << std::endl;
     }
     
     for(int p = 0; p < particles[6].size(); p++)
     {
         xyzFile << 6 << " " << p << " " << particles[6][p].pos.x << " " <<  particles[6][p].pos.y << " " << particles[6][p].pos.z << " " << particles[6][p].totalDisplacement.x / (double)numAttempts << " " << particles[6][p].totalDisplacement.y / (double)numAttempts << " " << particles[6][p].totalDisplacement.z / (double)numAttempts << std::endl;
     }
    
    
    
    xyzFile.close();
    occupancyFile.close();
    exoPositionsFile.close();
    endoPositionsFile.close();
    numEachStateFile.close();
    
    
}
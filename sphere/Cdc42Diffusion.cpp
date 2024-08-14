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
   This code diffuses (diffuseParticle), reacts (lines 719-1079; 1422-1504),
   and exerts the effects of instantaneous exo/endocytosis (lines 1173-1420)
   on structs of type Particle restricted to the surface of a sphere. 
   The simulations proceeds in discrete timesteps of dt or (for exo/endocytosis
   or Particle addition to the surface) del_t. The simulation runs until finalTime
   taking xyz snapshots every snapshotTime. 
   If given as a command line argument, positionFileName (an xyz frame) can be used
   to restart the simulation from a previous state rather  than initializing the
   Particles to random positions on the sphere surface.
   Uses physics definition of theta and phi, polar angle is theta, azimuthal angle is phi
   
   Particle types refer to the following species:
   0: Cdc42-GDP
   1: Cdc42-GTP
   2: Cdc42-GDP (cytoplasmic)
   3: Scd2/1
   4: Scd2/1 (cytoplasmic)
   5: Scd2/1/Cdc42-GTP
   6: sGAP
   7: sGAP (cytoplasmic)
   
   This code was originally compiled with the C++17 standard using the GNU compiler.
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

// calculate the straight line diestance    
double calcDistance(Coordinate posi, Coordinate posj)
{
    Coordinate rij = posj - posi;
    
    return rij.getMagnitude();
}

//calculate distance between posi and posj on sphere surface with radius sphereRadius
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

// diffuse a particle starting at position pos by 
// 1) temporarily placing particle at north pole
// 2) displacing by gaussDirectionX and gaussDirectionY in the tangent plane
// 3) conserving this distance and projected down to the sphere surface 
// 4) finally rotating the particle back by its initial theta and phi angles from the north pole 
Coordinate diffuseParticle(double time_term, double prefactor, Coordinate pos, double gaussDirectionX, double gaussDirectionY, double radius)
{
    double initialTheta = acos(pos.z / radius);
    double initialPhi = atan2(pos.y, pos.x);
 
    double delPosX, delPosY;
    
    delPosX = prefactor * time_term * gaussDirectionX;
    delPosY = prefactor * time_term * gaussDirectionY;
    
    double stepsize = sqrt(delPosX*delPosX + delPosY*delPosY);              
    
    // inspired by hankbesser on github
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
    double pi = 3.14159265359;
    double initialRadius = 1.8;
    
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
    
    int totalNumberOfCdc42 = 5000;
    int totalNumberScd2Scd1 = 225;
    int totalNumberGAP = totalNumberScd2Scd1*10;
    
    // command line read-in overwrites default parameters
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

    // implementation of three conditions listed in "Simulations with varying Cdc42 mobility" section
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
    
    // calculation of endo_surface_radius, the arclength of the hemispherical cap on the sphere surface
    // that gives an equivalent area as area_per_vesicle_endo, particles closer than this distance
    // to the endocytosis event are placed at the center of endocytosis event
    double height_spherical_cap = area_per_vesicle_endo / 2.0 / pi / initialRadius;
    double a = sqrt(area_per_vesicle_endo / pi - height_spherical_cap * height_spherical_cap);
    double endo_surface_radius = initialRadius * atan(a / (initialRadius - height_spherical_cap));
    
    // time related variables
    double currTime = 0.0;
    double finalTime = 5000.0;
    
    double nextSnapshotTime = 0.0;
    
    //time at which to add another snapshot to xyz output file
    double snapshotTime = 10.0;
    int frameModValue = (int)(snapshotTime / dt);

    int currFrame = 0;
    
    
    //https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
    std::random_device rd;
    std::mt19937 gen(time(0));
   
    double diffusion_prefactorGDP = sqrt(2.0 * D_D);
    double diffusion_prefactorGTP = sqrt(2.0 * D_T);
    double diffusion_prefactorScd1Scd2 = sqrt(2.0 * D_S);
    double diffusion_prefactorScd1Scd2GTP = sqrt(2.0 * D_ST);
    double diffusion_prefactorGAP = sqrt(2.0 * D_GAP);
    
    int maxNumTypes = 9;
    
    //vector of diffusion prefactors for different components to use in diffuseParticle function
    //-10.0 indicates a particle that is cytoplasmic and does not explicitly diffuse
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
        // reading in particle positions from positionFileName xyz file
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
                    std::string tempString;
                    iss >> tempString;
                    
                    int pos = tempString.find("=");
                    
                    std::string timeString = tempString.substr(pos+1,tempString.length());
                    
                    double currSimTime = std::stod(timeString);
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
        // random initialization of the particles on the surface of the sphere
        for(int i = 0; i < membraneNumberOfCdc42; i++)
        {
            double phiRnd = 2.0*pi*dis(gen);
            double thetaRnd = acos(2*dis(gen)-1);
            
            Coordinate c = {initialRadius*cos(phiRnd)*sin(thetaRnd), initialRadius*sin(phiRnd)*sin(thetaRnd), initialRadius*cos(thetaRnd)};
            
            particles[0].emplace_back(c, i, 0.0, 0.0);
        }
        
        for(int i = 0; i < cytoNumberOfCdc42; i++)
        {
            Coordinate c = {0.0, 0.0, 0.0};
            
            particles[2].emplace_back(c, i+membraneNumberOfCdc42, 0.0, 0.0);
        }
        
        for(int i = 0; i < membraneNumberOfScd2Scd1; i++)
        {
            double phiRnd = 2.0*pi*dis(gen);
            double thetaRnd = acos(2*dis(gen)-1);
            
            Coordinate c = {initialRadius*cos(phiRnd)*sin(thetaRnd), initialRadius*sin(phiRnd)*sin(thetaRnd), initialRadius*cos(thetaRnd)};
            
            particles[3].emplace_back(c, i+membraneNumberOfCdc42+cytoNumberOfCdc42, 0.0, 0.0);
        }
        
        for(int i = 0; i < cytoNumberOfScd2Scd1; i++)
        {
            Coordinate c = {0.0, 0.0, 0.0};
            
            particles[4].emplace_back(c, i+membraneNumberOfCdc42+cytoNumberOfCdc42+membraneNumberOfScd2Scd1, 0.0, 0.0);
        }
        
        for(int i = 0; i < membraneNumberOfGAP; i++)
        {
            double phiRnd = 2.0*pi*dis(gen);
            double thetaRnd = acos(2*dis(gen)-1);
            
            Coordinate c = {initialRadius*cos(phiRnd)*sin(thetaRnd), initialRadius*sin(phiRnd)*sin(thetaRnd), initialRadius*cos(thetaRnd)};
            
            particles[6].emplace_back(c, i+membraneNumberOfCdc42+cytoNumberOfCdc42+membraneNumberOfScd2Scd1+cytoNumberOfScd2Scd1, 0.0, 0.0);
        }
        for(int i = 0; i < cytoNumberOfGAP; i++)
        {
            Coordinate c = {0.0, 0.0, 0.0};
            
            particles[7].emplace_back(c, i+membraneNumberOfCdc42+cytoNumberOfCdc42+membraneNumberOfScd2Scd1+cytoNumberOfScd2Scd1+membraneNumberOfGAP, 0.0, 0.0);
        }
    }
    
    //output file setup
    std::ostringstream out;
    out.precision(5);
    out << std::fixed << r_D;
    
    std::ostringstream out2;
    out2.precision(5);
    out2 << std::fixed << D_D;
    
    std::string appendString = "_rD_" + out.str() + "_DD_" + out2.str() + "_Run1";

    std::ofstream xyzFile ("3d" + appendString + ".xyz");
    
    std::ofstream numEachStateFile ("NumEachState" + appendString + ".txt");
    
    std::ofstream occupancyFile ("3doccupancyHist" + appendString + ".dat");
    std::ofstream exoPositionsFile ("exoPositions" + appendString + ".dat");
    std::ofstream endoPositionsFile ("endoPositions" + appendString + ".dat");

    
    // setup of Gillespie algorithm
    int numRates = 5;
    std::vector <double> rates(numRates);

    std::vector <double> Rki(numRates);

    std::vector <int> countEvents(numRates);
    
    //Gillespie exo and endocytosis rates
    rates[0] = exo_rate;
    rates[1] = endo_rate;
    
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
        
        // time until next reaction from Gillespie algorithm
        double del_t = -1.0 / Rki[numRates-1] * log(dis(gen));
        
        // advance to next snapshot time while currTime + next Gillspie event time is not yet reached
        while(currTime + del_t >= nextSnapshotTime)
        {
            currFrame++;
            
            int tempOccupancyTip = 0;
            int tempOccupancyBack = 0;
    
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

                        particles[i][p].pos = diffuseParticle(time_term, curr_diffusionPrefactor, particles[i][p].pos, gaussDirectionX, gaussDirectionY, initialRadius);
                        
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
                                    final_pos3D = final_pos3D.getUnitCoord();
                                    final_pos3D = initialRadius * final_pos3D;
                                    
                                    
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
                            
                            // inspired by hankbesser on github
                            double theta = unbinding_rad *0.5 / initialRadius;
                            double phi = 2.0*pi*dis(gen);
                            
                            double initialTheta = acos(initialParticle.pos.z / initialRadius);
                            double initialPhi = atan2(initialParticle.pos.y, initialParticle.pos.x);
                    
                            // this is rotated position assuming that the particle started at the north pole (0,0,initialradius)
                            struct Coordinate final_pos3D;
                            final_pos3D.x = initialRadius*cos(phi)*sin(theta);
                            final_pos3D.y = initialRadius*sin(phi)*sin(theta);
                            final_pos3D.z = initialRadius*cos(theta);
                            
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
                            
                            particles[3].emplace_back(final_pos3D, initialParticle.secondname, initialParticle.timeLastGlobalUpdate, initialParticle.posTime);
                            
                            
                            theta = -theta;
                    
                            // this is rotated position assuming that the particle started at the north pole (0,0,initialradius)
                            final_pos3D.x = initialRadius*cos(phi)*sin(theta);
                            final_pos3D.y = initialRadius*sin(phi)*sin(theta);
                            final_pos3D.z = initialRadius*cos(theta);
                            
                            rotK.x = 0.0;
                            rotK.y = 1.0;
                            rotK.z = 0.0;
                            
                            // rotate by initialTheta
                            currDotProduct = rotK*final_pos3D;
                            currCrossProduct = rotK.crossProduct(final_pos3D);
                            
                            final_pos3D = final_pos3D*cos(initialTheta) + currCrossProduct*sin(initialTheta) + rotK*currDotProduct*(1.0 - cos(initialTheta));
                            
                            rotK.x = 0.0;
                            rotK.y = 0.0;
                            rotK.z = 1.0;
                            
                            // rotate by initialPhi
                            currDotProduct = rotK*final_pos3D;
                            currCrossProduct = rotK.crossProduct(final_pos3D);
                            
                            final_pos3D = final_pos3D*cos(initialPhi) + currCrossProduct*sin(initialPhi) + rotK*currDotProduct*(1.0 - cos(initialPhi));
                            
                            particles[1].emplace_back(final_pos3D, initialParticle.name, initialParticle.timeLastGlobalUpdate, initialParticle.posTime);
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
                            
                            removedParticle = true;
                        }
                    }
                    
                    if(removedParticle == false)
                    {
                        particles[i][p].timeLastGlobalUpdate = nextSnapshotTime;
                    }
                }
            }
            
            // write current configuration of particles if it is time to do so based on snapshotTime
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
                xyzFile << "t=" << nextSnapshotTime << endl;
                
                numEachStateFile << nextSnapshotTime << ",";
                
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
                
                // add particle in xyz file to represent center of sphere
                xyzFile << 10 << " " << -1 << " " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
                
                numEachStateFile << endl;
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
        
        // now it is time to do a single Gillespie algorithm step
        // find out which event to select based for Gillespie algorithm
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
        
        
        //only need to update all particles if the event type is exo or endocytosis
        //update position of particles (diffusion)
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
    
                        particles[i][p].pos = diffuseParticle(time_term, curr_diffusionPrefactor, particles[i][p].pos, gaussDirectionX, gaussDirectionY, initialRadius);
                        
                        particles[i][p].posTime = del_t + currTime;
                    }
                }
            }
        }
        
        if(eventIndex == 0)
        {
            //exocytosis
            double phiRnd, thetaRnd;
            Coordinate c;
            double randomUniform = 1.0;
            
            double inv_sigma_distribution = 1.0 / w_exo;
            double inv_sigma_distribution_sq = inv_sigma_distribution*inv_sigma_distribution;
            
            // using positions of Cdc42-GTP particles to find a position for the exocytosis event
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
                    
                    prob_value = prob_value / (double)particles[1].size();
                        
                    
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
            
            
            Coordinate exclusionCenter = c;
            
            exoPositionsFile << exclusionCenter.x << " " << exclusionCenter.y << " " << exclusionCenter.z << endl;
            
            Coordinate v1 = exclusionCenter.getUnitCoord();
            
            //push particles out of exclusion zone according to equation in "Membrane flows" section
            for(int i = 0; i < particles.size(); i++)
            {
                //only move particles that are on the surface
                if(diffusionPrefactors[i] >= 0.0)
                {
                    for(int p = particles[i].size()-1; p > -1; p--)
                    {                                            
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
                                // temp_exclusion_radius is additional distance away from the exocytotic center that the particle should move                                
                                double RSphereNew = sqrt(initialRadius*initialRadius + area_per_vesicle_exo / (4.0 * pi));
                                
                                double acosParameter = 1.0 - initialRadius*initialRadius / (RSphereNew*RSphereNew) * (1.0 - cos(var_l / initialRadius)) - area_per_vesicle_exo / (2.0 * pi * RSphereNew*RSphereNew * alpha);
                                if(acosParameter < -1.0)
                                {
                                    acosParameter = -1.0;
                                }
                                else if(acosParameter > 1.0)
                                {
                                    acosParameter = 1.0;
                                }
                                
                                double delR = RSphereNew * acos(acosParameter) - var_l;
                                double var_t = (var_l + delR) / initialRadius;
                                
                                if(var_t > pi)
                                    var_t = pi;
                                
                                Coordinate finalPos = initialRadius*cos(var_t)*v1 + initialRadius*sin(var_t)*v2;

                                particles[i][p].pos = finalPos;
                            }
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
            
            double inv_sigma_distribution = 1.0 / w_endo;
            double inv_sigma_distribution_sq = inv_sigma_distribution*inv_sigma_distribution;
            
            // using positions of Cdc42-GTP particles to find a position for the endocytosis event
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
                    
                    prob_value = prob_value / (double)particles[1].size();
                        
                    
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
            
            //pull particles into inclusion zone according to equation in "Membrane flows" section
            for(int i = 0; i < particles.size(); i++)
            {
                //only move particles that are on the surface
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
                                // place particle at center of endocytosis event if it is within the arclength of the hemispherical cap
                                // associated with the endocytosis vesicle area
                                particles[i][p].pos = inclusionCenter;
                            }
                            else if(var_l < var_l_endo_compare)
                            {
                                Coordinate v2 = particles[i][p].pos - inclusionCenter;
                                
                                v2 = v2 - (v2*v1)*v1;
                                
                                v2 = v2.getUnitCoord();
                                
                                // angle made between exclusion center and final position
                                
                                double RSphereNew = sqrt(initialRadius*initialRadius - area_per_vesicle_endo / (4.0 * pi));
                                
                                double acosParameter = 1.0 - initialRadius*initialRadius / (RSphereNew*RSphereNew) * (1.0 - cos(var_l / initialRadius)) + area_per_vesicle_endo / (2.0 * pi * RSphereNew*RSphereNew * alpha);
                                if(acosParameter < -1.0)
                                {
                                    acosParameter = -1.0;
                                }
                                else if(acosParameter > 1.0)
                                {
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
        }
        
        else if(eventIndex == 2)
        {
            //Cdc42 proteinAdd
            double phiRnd = 2.0*pi*dis(gen);
            double thetaRnd = acos(2*dis(gen)-1);
            
            Coordinate c = {initialRadius*cos(phiRnd)*sin(thetaRnd), initialRadius*sin(phiRnd)*sin(thetaRnd), initialRadius*cos(thetaRnd)};

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
            double phiRnd = 2.0*pi*dis(gen);
            double thetaRnd = acos(2*dis(gen)-1);
            
            Coordinate c = {initialRadius*cos(phiRnd)*sin(thetaRnd), initialRadius*sin(phiRnd)*sin(thetaRnd), initialRadius*cos(thetaRnd)};
            
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
            double phiRnd = 2.0*pi*dis(gen);
            double thetaRnd = acos(2*dis(gen)-1);
            
            Coordinate c = {initialRadius*cos(phiRnd)*sin(thetaRnd), initialRadius*sin(phiRnd)*sin(thetaRnd), initialRadius*cos(thetaRnd)};
            
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
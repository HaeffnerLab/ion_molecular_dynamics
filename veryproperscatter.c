#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>


typedef int bool;
enum { false, true };

#define DBL_MAX 1e37

#define Npartmax 2000    // dimension of arrays ( Npart < Npartmax )
#define MAXSWAPS 100000


typedef struct{
  double x, y, z;
} position;

typedef struct{
  double vx, vy, vz;
} velocity;

typedef struct{
  double px, py, pz;
} momentum;

typedef struct{
  double ax, ay, az;
} acceleration;

int Npart = 5;  // No of particles
double mass = 6.6551077e-26; // in kg
position R[Npartmax];       // positions
velocity V[Npartmax];       // velocities
acceleration acc[Npartmax]; // accelerations
FILE *fp_traj, *fp_traj_xyz, *fp_traj_swap, *fp_traj_readout;
position Rp[Npartmax], Rnew[Npartmax];
double dt, tmin, tmax;      // time grid
double trap_freq_x, trap_freq_y, trap_freq_z;
double boltzmann = 1.3806503e-23;
double hbar = 1.05457148e-34;
double initialTemp = 0.0005; // Kelvin
double decaytime = 0.0000000071;//0.0000000071; // decay time for working transition in Ca -> 7.1ns
double detuning_parameter = -.5;
double detuning;
double s_beam = 1.0; // saturation parameter
double laser_direction_temp[] = {1.0, 1.0, 1.0}; //laser comes in the x direction only // laser comes in at 45 degree angle
double laser_direction[] = {0.0, 0.0, 0.0};
double gamma_param = 140845070; //1/(7.1ns) //22416189.2; // (1/(2*Pi*7.1ns) in hertz // 1.319e8;
double laser_velocity = 0.024877240126316254; // hbar*k/mass
double rho_ee;
double timeStepFactor;
double theta_initial, phi_initial, theta_recoil, phi_recoil;
//int *currentlyExcitedPart = NULL; // array that keeps track of excited particles (1 == excited, 0 == ground)
int *decayState = NULL;
double *excitedTimer = NULL; // array to keep track of the atom's internal timer during absorption (prior to emmission)
double *lifetimes = NULL;           //  decay Times for getting average natural lifetime
double *time_in_excited_state = NULL;
double *time_in_ground_state = NULL;
double t; // official simulation time
double initialTime, finalTime;
double gamma_lasers; // sum of all the calculated gamma_lasers, used for averaging

// swap stuff
int swapCounter = 0;
int swapsToMelt = 2000;
int swapCatalog[MAXSWAPS][2]; // 5 swaps, 1st column initial, 2nd final
int swapCatalogAbs[MAXSWAPS][2]; // 5 swaps, 1st column initial, 2nd final
int ionTracking[Npartmax]; // keep track of an ion's position relative to others in the chain (first column is initial ion, second is its relative position)

// micromotion - these values correspond to trap frequencies of around 1MHz, 1Mhz and 1KHz for x y and z (see mathematica notebook) s

double V_1 = 45.55362938;//14.5;//25.0;
double V_0 = 41.94377859;//21.1;//84.3;//13.8;
double U_0 = 16.39851231;//3.745;//.4;
double R_0 = 0.0005;
double Z_0 = 0.01;
double kappa = 0.5;
double trap_drive_omega;
double electron_charge = 1.602e-19;

// used for averaging and plotting.
int binCount;      // a counter to keep track of the number of RF periods
int rfCount;       // a counter to keep track of time steps within an RF period
double bin;        // number of sets of RF periods CRUCIAL
double numRFsteps; // number of time steps in one RF period CRUCIAL
double *emit_counter = NULL;
double *total_emit_counter = NULL;
double *total_check_counter = NULL;
double *velocitiesX = NULL;
double *velocitiesY = NULL;
double *velocitiesZ = NULL;
double *velocitiesSquaredX = NULL;
double *velocitiesSquaredY = NULL;
double *velocitiesSquaredZ = NULL;
double *tempavg = NULL;
double *tempavgx = NULL;
double *tempavgy = NULL;
double *tempavgz = NULL;
double *fluoravg = NULL;
double *totalAvg = NULL;

bool startReadout = false;

// Can't define these globally since you can't run a function in a global scope, so assigning by pointers instead
// these are set inside 'initial_conditions', make sure to do calculations with (*variable) not (variable)
double *lattice_power;
double *damping;
double *k;
double *recoil_velocity = NULL;
double *k_laser = NULL;

// specific experiment

int laser_counter = 0;
                         // variables for averaging and plottings
void setup_averaging_variables()
{
    
    int i;
    binCount = 0;
    bin = 100.0; // NUMBER OF SETS OF RF STEPS
    numRFsteps = 100.0; // NUMBER OF RF STEPS
    rfCount = 0;
    emit_counter = malloc(Npart*sizeof(double));
    total_emit_counter = malloc(Npart*sizeof(double));
    total_check_counter = malloc(Npart*sizeof(double));
    velocitiesX = malloc(Npart*sizeof(double));
    velocitiesY = malloc(Npart*sizeof(double));
    velocitiesZ = malloc(Npart*sizeof(double));
    velocitiesSquaredX = malloc(Npart*sizeof(double));
    velocitiesSquaredY = malloc(Npart*sizeof(double));
    velocitiesSquaredZ = malloc(Npart*sizeof(double));
    tempavg = malloc(Npart*sizeof(double));
    tempavgx = malloc(Npart*sizeof(double));
    tempavgy = malloc(Npart*sizeof(double));
    tempavgz = malloc(Npart*sizeof(double));

    fluoravg = malloc(Npart*sizeof(double));
    totalAvg = malloc(Npart*sizeof(double));


    for (i = 0; i<Npart; i++)
    {
	    emit_counter[i] = 0.0;
            total_emit_counter[i] = 0.0;
	    total_check_counter[i] = 0.0;
	    velocitiesX[i] = 0.0;
	    velocitiesY[i] = 0.0;
	    velocitiesZ[i] = 0.0;
	    velocitiesSquaredX[i] = 0.0;
	    velocitiesSquaredY[i] = 0.0;
	    velocitiesSquaredZ[i] = 0.0;
	    tempavg[i] = 0.0;
            tempavgx[i] = 0.0;
            tempavgy[i] = 0.0;
            tempavgz[i] = 0.0;
	    fluoravg[i] = 0.0;
	    totalAvg[i] = 0.0;
     }
}

                            // output routine
void output_x_v()
{
  int i;

  for ( i=0 ; i<Npart ; i++ )
    {
      fprintf( stderr, "Position & velocity:  %d %.20f %.20f %.20f %.20f %.20f %.20f\n", 
           i, R[i].x, R[i].y, R[i].z, V[i].vx, V[i].vy, V[i].vz ); 
      double temp;
      temp = (V[i].vx*V[i].vx +  V[i].vy*V[i].vy + V[i].vz*V[i].vz)*mass/(2*boltzmann);
      printf("Temperature: %.30f\n", temp);
    }
}

                            // record trajectory in movie format!
void record_trajectories_xyz( )
{
  int i;
  double x, y, z;

  fprintf(fp_traj_xyz, "%d\n", Npart);
  fprintf(fp_traj_xyz, "Mark Kokish Test XYZ\n");
  
  for ( i=0 ; i<Npart ; i++ )
    {
      x = R[i].x*1000000.0;
      y = R[i].y*1000000.0;
      z = R[i].z*1000000.0;
      if (i == 1 || i==4)
      {
        fprintf(fp_traj_xyz, "Cr %.30f %.30f %.30f\n", x, y, z);
      }
      else
      {
	fprintf(fp_traj_xyz, "V %.30f %.30f %.30f\n", x, y, z);
      }
 
    }
}


                            // record trajectory in file
void record_trajectories( )
{
  int i;
  double xtemp, ytemp;

  for ( i=0 ; i<Npart ; i++ )
  {
    velocitiesX[i] += V[i].vx;
    velocitiesY[i] += V[i].vy;
    velocitiesZ[i] += V[i].vz;
  }

  if (rfCount == numRFsteps)
  {
    for ( i=0 ; i<Npart ; i++ )
    {

      velocitiesSquaredX[i] += (velocitiesX[i]/numRFsteps)*(velocitiesX[i]/numRFsteps);
      velocitiesSquaredY[i] += (velocitiesY[i]/numRFsteps)*(velocitiesY[i]/numRFsteps);
      velocitiesSquaredZ[i] += (velocitiesZ[i]/numRFsteps)*(velocitiesZ[i]/numRFsteps);
      velocitiesX[i] = 0.0;
      velocitiesY[i] = 0.0;
      velocitiesZ[i] = 0.0;
    }
 
    rfCount = 0;
    binCount++;
//    record_trajectories_xyz(); // record movie every RF cycle!!
  }

  if (binCount == bin)
  {
    if (startReadout == false)
    {
      fprintf(fp_traj, "%.30f ", t*1e3); // milliseconds
    }
    else
    {
      fprintf(fp_traj_readout, "%.30f ", t*1e3); // milliseconds
    }
    for ( i=0 ; i<Npart ; i++ )
    {

	tempavg[i] = 2.0*mass*(velocitiesSquaredX[i]/bin+velocitiesSquaredY[i]/bin+velocitiesSquaredZ[i]/bin)/(3.0*hbar*gamma_param);
        tempavgx[i] = 2.0*mass*(velocitiesSquaredX[i]/bin)/(3.0*hbar*gamma_param);
        tempavgy[i] = 2.0*mass*(velocitiesSquaredY[i]/bin)/(3.0*hbar*gamma_param);
        tempavgz[i] = 2.0*mass*(velocitiesSquaredZ[i]/bin)/(3.0*hbar*gamma_param);
	totalAvg[i] += tempavg[i];
	fluoravg[i] = ((emit_counter[i]/(bin*numRFsteps*dt))/(1e6)); // MegaCounts per sec
	
        if (startReadout == false)
        {
	  fprintf(fp_traj, " %.30f %.30f %.30f %.30f %.30f %.30f %.30f %.30f %.30f %.30f %.30f ", tempavg[i], fluoravg[i], tempavgx[i], tempavgy[i], tempavgz[i], R[i].x, R[i].y, R[i].z, V[i].vx, V[i].vy, V[i].vz);
        }
        else
        {
          fprintf(fp_traj_readout, " %.30f ", fluoravg[i]);
        }
//	fprintf(fp_traj, " %.30f %.30f ", tempavg[i], fluoravg[i]);
	emit_counter[i] = 0.0;
	velocitiesSquaredX[i] = 0.0;
	velocitiesSquaredY[i] = 0.0;
	velocitiesSquaredZ[i] = 0.0;
	binCount = 0;
    }	
    if (startReadout == false)
    {
      fprintf(fp_traj, "\n");
    }
    else
    {
      fprintf(fp_traj_readout, "\n");
    }
  }

  

}

void record_trajectories_readout( )
{
  int i;
  double xtemp, ytemp;

  fprintf(fp_traj_readout, "%.30f ", t*1e3); // milliseconds
  
  for ( i=0 ; i<Npart ; i++ )
  {

    velocitiesX[i] = V[i].vx;
    velocitiesY[i] = V[i].vy;
    velocitiesZ[i] = V[i].vz;

    velocitiesSquaredX[i] = (velocitiesX[i])*(velocitiesX[i]);
    velocitiesSquaredY[i] = (velocitiesY[i])*(velocitiesY[i]);
    velocitiesSquaredZ[i] = (velocitiesZ[i])*(velocitiesZ[i]);
    velocitiesX[i] = 0.0;
    velocitiesY[i] = 0.0;
    velocitiesZ[i] = 0.0;

    tempavg[i] = 2.0*mass*(velocitiesSquaredX[i]+velocitiesSquaredY[i]+velocitiesSquaredZ[i])/(3.0*hbar*gamma_param);
    totalAvg[i] += tempavg[i];
    fluoravg[i] = ((emit_counter[i]/(dt))/(1e6)); // MegaCounts per sec
	
    fprintf(fp_traj_readout, " %.30f %.30f ", tempavg[i], fluoravg[i]);
    emit_counter[i] = 0.0;
    velocitiesSquaredX[i] = 0.0;
    velocitiesSquaredY[i] = 0.0;
    velocitiesSquaredZ[i] = 0.0;
  }	
  fprintf(fp_traj_readout, "\n");



}


void initial_conditions( )
{
  srand(time(NULL)); // random setup
  double aspect_ratio, ee, epsilon_0, length_scale, energy_scale, initialVelocity;
  int i, j;

  setup_averaging_variables();

  double initialZPositions[] = {-0.00003592197686982309,-0.00001694384644835323, -0.00000000000003984558,0.00001694384636012078,0.00003592197676604132};//{-0.00012400338910260807, -0.00010032968053235605, -0.00008071905149605598, -0.00006305777988679104, -0.00004652649444710114, -0.00003068542210321371, -0.00001524865657422565, -0.00000000000000001103, 0.00001524865657420344, 0.00003068542210319132, 0.00004652649444707858, 0.00006305777988676827, 0.00008071905149603312, 0.00010032968053233323, 0.00012400338910258609};
 
  double omega_para[] = {483.438,478.828,47.786};//{158.280,166.347,15.617};

  trap_drive_omega = 2.0*M_PI*30000000.0; // 2.0*M_PI*15000000
  
//  k = malloc(sizeof(double));
//  *k = 2.0*M_PI/(2.0*3.075430571317475407e-04/(1.0+sqrt(5.0)));
  
  aspect_ratio = omega_para[2]/omega_para[0];
  ee = (1.60217646e-19)*(1.60217646e-19);
  epsilon_0 = 8.8542e-12;

  length_scale = pow(((ee*aspect_ratio)*(ee*aspect_ratio))/(4.0*M_PI*epsilon_0*mass*(omega_para[2]*1000.0)*(omega_para[2]*1000.0)), (1.0/3.0));

//  energy_scale = ee/(4.0*M_PI*epsilon_0*length_scale);
  lattice_power = malloc(sizeof(double));//*0.020*energy_scale; // lattice power is off for now  
  *lattice_power = 0.0;
  damping = malloc(sizeof(double));
  *damping = 0.0; //1.0e-19;  //1e-18 is overdamped // damping is off for now
   
//  trap_freq_x = 2.0*M_PI*1000.0*omega_para[0];
//  trap_freq_y = 2.0*M_PI*1000.0*omega_para[1];
//  trap_freq_z = 2.0*M_PI*1000.0*omega_para[2];
  
//  currentlyExcitedPart = malloc(Npart*sizeof(int));
  decayState = malloc(Npart*sizeof(int));
  excitedTimer = malloc(Npart*sizeof(double));
  lifetimes = malloc(Npart*sizeof(double));
  time_in_excited_state = malloc(Npart*sizeof(double));
  time_in_ground_state = malloc(Npart*sizeof(double));
  
  // laser stuff
  
  double laser_direction_temp_magnitude;
  
  laser_direction_temp_magnitude = sqrt(laser_direction_temp[0]*laser_direction_temp[0] + laser_direction_temp[1]*laser_direction_temp[1] + laser_direction_temp[2]*laser_direction_temp[2]); 
  
  // establish proper sizes for laser_direction and k_laser
 
  recoil_velocity = malloc(3*sizeof(double));
  k_laser = malloc(3*sizeof(double));
  for (j = 0; j<3; j++)
  {
    laser_direction[j] = laser_direction_temp[j]/laser_direction_temp_magnitude;
    recoil_velocity[j] = laser_velocity*laser_direction[j];
    k_laser[j] = (M_PI*2.0/(396.847e-9))*laser_direction[j];
  }
 
  detuning = detuning_parameter*gamma_param;
  initialVelocity = sqrt(2.0*boltzmann*initialTemp/mass);
  
  for ( i=0 ; i<Npart ; i++ )
    {
      
                     // initial positions
      R[i].x = 0.0;
      R[i].y = 0.0;
      R[i].z = initialZPositions[i];

      // the ions are intially already labelled in order
      ionTracking[i] = i;


//      currentlyExcitedPart[i] = 0;
      decayState[i] = 0;
      excitedTimer[i] = 0.0;
      lifetimes[i] = 0.0;
      time_in_excited_state[i] = 0.0;
      time_in_ground_state[i] = 0.0;
      
   
                       // random velocities
                       // random number times something on the order of kbT
  
      theta_initial = acos(2.0*((double)rand()/(double)RAND_MAX)-1.0);
      phi_initial = 2.0*M_PI*((double)rand()/(double)RAND_MAX);
      V[i].vx = initialVelocity*sin(theta_initial)*cos(phi_initial);
      V[i].vy = initialVelocity*sin(theta_initial)*sin(phi_initial);
      V[i].vz = initialVelocity*cos(theta_initial);



                        // previous value (needed in Verlet step) Rp == "R Previous"
                        // Euler step
      Rp[i].x = R[i].x - dt*V[i].vx;
      Rp[i].y = R[i].y - dt*V[i].vy;
      Rp[i].z = R[i].z - dt*V[i].vz;
    }
    
    printf("First Ion Tracking : ");
    for (j = 0; j<Npart; j++)
    {
      printf("%d, ",  ionTracking[j]);
    }
    printf("\n");
    
//     R[7].x = 0.000058062228770773562360960063;
//     Rp[7].x = R[7].x - dt*V[7].vx;
    
  printf("Initial velocity: %.30f %.30f %.30f\n", V[0].vx, V[0].vy, V[0].vz );  
}


                            // force (acceleration) due to Coulomb potential
                            // on each body caused by
                            // all other bodies
                            
                            // All the forces! harmonic and Coulomb
void accelerations( int iprint )
{
  int i, j;
  double xij, yij, zij;
  double Fijx, Fijy, Fijz;
  double Fx[Npart], Fy[Npart], Fz[Npart];
  double factor;
  double magnitude;

  // micromotion
  
  for ( i=0 ; i<Npart ; i++ )
  {
  
    Fx[i] =  electron_charge*( (-4.0)*V_0*cos(trap_drive_omega*t)*R[i].x/(R_0*R_0) + kappa*U_0*R[i].x/(Z_0*Z_0) );// - (*damping)*V[i].vx;
    Fy[i] =  electron_charge*( (4.0)*V_1*cos(trap_drive_omega*t)*R[i].y/(R_0*R_0) + kappa*U_0*R[i].y/(Z_0*Z_0) );// - (*damping)*V[i].vy;
    Fz[i] =  electron_charge*( (-2.0)*kappa*U_0*R[i].z/(Z_0*Z_0) );// - (*damping)*V[i].vz + (*lattice_power)*(*k)*sin(R[i].z*(*k));
   
  }
  
  
  // harmonic pseudopotential
//   for ( i=0 ; i<Npart ; i++ )
//   {
//     Fx[i] = -R[i].x*mass*(trap_freq_x*trap_freq_x);// - (*damping)*V[i].vx;
// //     printf("Force on x: %.30f\n", Fx[i]);
// //     printf("Velocity in x: %.30f\n", V[i].vx);
//     Fy[i] = -R[i].y*mass*(trap_freq_y*trap_freq_y);// - (*damping)*V[i].vy;
//     Fz[i] = -R[i].z*mass*(trap_freq_z*trap_freq_z);// - (*damping)*V[i].vz + (*lattice_power)*(*k)*sin(R[i].z*(*k));
//   }
  
  // Coulomb
  
  for ( i=0 ; i<Npart-1 ; i++ )
    {
    //  printf("i: %d\n", i);
      for ( j= i+1 ; j<Npart ; j++ )
	{

	  magnitude = sqrt((R[i].x - R[j].x)*(R[i].x - R[j].x) + (R[i].y - R[j].y)*(R[i].y - R[j].y) + (R[i].z - R[j].z)*(R[i].z - R[j].z));
              xij = (R[i].x - R[j].x)/fabs(pow(magnitude, 3.0));
              yij = (R[i].y - R[j].y)/fabs(pow(magnitude, 3.0));
	      zij = (R[i].z - R[j].z)/fabs(pow(magnitude, 3.0));

	      Fx[i] = Fx[i] + 2.30e-28*xij;
	      Fy[i] = Fy[i] + 2.30e-28*yij;
	      Fz[i] = Fz[i] + 2.30e-28*zij;
              Fx[j] = Fx[j] - 2.30e-28*xij;
	      Fy[j] = Fy[j] - 2.30e-28*yij;
	      Fz[j] = Fz[j] - 2.30e-28*zij;
	      

	} // end j loop
    } // i loop 

   for ( i=0 ; i<Npart ; i++ )
    {
      acc[i].ax = Fx[i] / mass;
      acc[i].ay = Fy[i] / mass;
      acc[i].az = Fz[i] / mass;
    }

}


double get_excitation_probability( int ion )
{
  // calculate doppler
  double doppler;
  double probability;
  double gamma_laser;
  double s;

  doppler = ((k_laser[0])*V[ion].vx + (k_laser[1])*V[ion].vy  + (k_laser[2])*V[ion].vz); 
  s = s_beam/(pow((2.0*(detuning-doppler)/gamma_param), 2.0) + 1.0);
  gamma_laser = gamma_param*s/2.0;
  gamma_lasers += gamma_laser;
//  printf("%.30f\n", gamma_laser);
  probability = gamma_laser*dt;
  
  return probability;
}

void catalogIonSwaps()
{
  // This function finds out if ions have swapped places and catalogs which ion moved where

  int i;
  int j;
  int count;


  int ionTrackingTemp[Npart];


  for (i=0; i<Npart; i++)
  {
    count = 0;
    for (j = 0; j<Npart; j++)
    { 
        if (Rnew[i].z < Rnew[j].z)
        {
          count++;
        }
    }
    ionTrackingTemp[i] = Npart - count - 1;
  }

  for (i=0; i<Npart; i++)
  {
    if (ionTrackingTemp[i] != ionTracking[i])
    {
      swapCatalog[swapCounter][0] = ionTrackingTemp[i]; //i;
      swapCatalog[swapCounter][1] = ionTracking[i];//indexToRemove;   
//      printf("Swap! %d %.30f\n", swapCounter, t);
      fprintf(fp_traj_swap, "%.30f %d %d\n", t, swapCatalog[swapCounter][0], swapCatalog[swapCounter][1] );
      swapCounter++; 
    }
    
    ionTracking[i] = ionTrackingTemp[i];    
  }
//*/  
}

                            // Velocity Verlet time step
time_step_verlet( int iprint )
{
  int i;
  double emit_rand_temp;
  double excitation_prob; 
                            // accelerations at t
  accelerations( iprint );
  

                            // basic Verlet
  for ( i=0 ; i<Npart ; i++ )
  {
    // check if a particle is already excited or not
    if (decayState[i] == 0)
    { 
      // if we're here, the particle is not currently excited
      // check for excitation
      time_in_ground_state[i] += dt;
      total_check_counter[i]++;
      excitation_prob = get_excitation_probability(i);
      if ( (double)rand()/(double)RAND_MAX < excitation_prob )
      {
	// excite! change decay state
	decayState[i] = 1;
	// start the excited timer!
        excitedTimer[i] += dt;

	// absorb!
	Rnew[i].x = 2.0*R[i].x - Rp[i].x + recoil_velocity[0]*dt + acc[i].ax*dt*dt;
	Rnew[i].y = 2.0*R[i].y - Rp[i].y + recoil_velocity[1]*dt + acc[i].ay*dt*dt;
	Rnew[i].z = 2.0*R[i].z - Rp[i].z + recoil_velocity[2]*dt + acc[i].az*dt*dt;

/*
	else
	{
	  // absorb AND emit (time step greater than decay time)
	  emit_counter[i]++;
	  phi_recoil = acos(2.0*((double)rand()/(double)RAND_MAX)-1.0);
	  theta_recoil = M_PI*((double)rand()/(double)RAND_MAX);

	  Rnew[i].x = 2.0*R[i].x - Rp[i].x + ( (laser_velocity*sin(theta_recoil)*cos(phi_recoil) + recoil_velocity[0])*dt ) + acc[i].ax*dt*dt;
//	    printf("recoilvelocity and kick in x (should not happen): %.30f\n", recoil_velocity[0]+laser_velocity*sin(theta_recoil)*cos(phi_recoil));
	  Rnew[i].y = 2.0*R[i].y - Rp[i].y + ( (laser_velocity*sin(theta_recoil)*sin(phi_recoil) + recoil_velocity[1])*dt ) + acc[i].ay*dt*dt;
	  Rnew[i].z = 2.0*R[i].z - Rp[i].z + ( (laser_velocity*cos(theta_recoil) + recoil_velocity[2])*dt ) + acc[i].az*dt*dt;	    
	  
	  // reset decay timer
	  decaytimer = 0.0;
	}
*/	    
	V[i].vx = ( Rnew[i].x - Rp[i].x ) / ( 2.0*dt );
	V[i].vy = ( Rnew[i].y - Rp[i].y ) / ( 2.0*dt );
	V[i].vz = ( Rnew[i].z - Rp[i].z ) / ( 2.0*dt );
	 
      }
      else
      {
	// no excitation evolve normally
	
	Rnew[i].x = 2.0*R[i].x - Rp[i].x + acc[i].ax*dt*dt;
	Rnew[i].y = 2.0*R[i].y - Rp[i].y + acc[i].ay*dt*dt;
	Rnew[i].z = 2.0*R[i].z - Rp[i].z + acc[i].az*dt*dt;
	
	V[i].vx = ( Rnew[i].x - Rp[i].x ) / ( 2.0*dt );
	V[i].vy = ( Rnew[i].y - Rp[i].y ) / ( 2.0*dt );
	V[i].vz = ( Rnew[i].z - Rp[i].z ) / ( 2.0*dt );

      }
    }

    else
    {
      // if we're here, that means the ion is already excited!
      time_in_excited_state[i] += dt;
      emit_rand_temp = (double)rand()/(double)RAND_MAX;
      if (emit_rand_temp < (excitation_prob + gamma_param*dt) ) // check if it's time to emit for this ion
      {
	// yes! do a verlet step that includes the emission
	if (emit_rand_temp < gamma_param*dt)
        {
	  emit_counter[i]++;
          total_emit_counter[i]++;
	  lifetimes[i] += excitedTimer[i]; // add time spent in excited to the total lifetimes

  	  phi_recoil = acos(2.0*((double)rand()/(double)RAND_MAX)-1.0);
	  theta_recoil = 2.0*M_PI*((double)rand()/(double)RAND_MAX);

	  Rnew[i].x = 2.0*R[i].x - Rp[i].x + ( (laser_velocity*sin(theta_recoil)*cos(phi_recoil))*dt ) + acc[i].ax*dt*dt;
	  Rnew[i].y = 2.0*R[i].y - Rp[i].y + ( (laser_velocity*sin(theta_recoil)*sin(phi_recoil))*dt ) + acc[i].ay*dt*dt;
	  Rnew[i].z = 2.0*R[i].z - Rp[i].z + ( (laser_velocity*cos(theta_recoil))*dt ) + acc[i].az*dt*dt;
        }
        else
        {
          // stimulated emission means recoil in the direciton opposite the laser
	  Rnew[i].x = 2.0*R[i].x - Rp[i].x - recoil_velocity[0]*dt + acc[i].ax*dt*dt;
	  Rnew[i].y = 2.0*R[i].y - Rp[i].y - recoil_velocity[1]*dt + acc[i].ay*dt*dt;
	  Rnew[i].z = 2.0*R[i].z - Rp[i].z - recoil_velocity[2]*dt + acc[i].az*dt*dt;
	}

	V[i].vx = ( Rnew[i].x - Rp[i].x ) / ( 2.0*dt );
	V[i].vy = ( Rnew[i].y - Rp[i].y ) / ( 2.0*dt );
	V[i].vz = ( Rnew[i].z - Rp[i].z ) / ( 2.0*dt );
	
	decayState[i] = 0; // go back to ground state
//        lifetimes[i] += excitedTimer[i]; // add time spent in excited to the total lifetimes
        excitedTimer[i] = 0.0; // set the timer back to 0
      }
      else
      {
	excitedTimer[i] += dt;
	// if we're here, that means we did not emit, so just evolve normally
	Rnew[i].x = 2.0*R[i].x - Rp[i].x + acc[i].ax*dt*dt;
	Rnew[i].y = 2.0*R[i].y - Rp[i].y + acc[i].ay*dt*dt;
	Rnew[i].z = 2.0*R[i].z - Rp[i].z + acc[i].az*dt*dt;
	
	V[i].vx = ( Rnew[i].x - Rp[i].x ) / ( 2.0*dt );
	V[i].vy = ( Rnew[i].y - Rp[i].y ) / ( 2.0*dt );
	V[i].vz = ( Rnew[i].z - Rp[i].z ) / ( 2.0*dt );
	
      }
    }  
  }

    
  // here you now certainly have the old and new positions R and Rnew
  // catalog ion swapping      
  catalogIonSwaps();

                            // reset arrays
  for ( i=0 ; i<Npart ; i++ )
  {
    Rp[i].x = R[i].x;
    Rp[i].y = R[i].y;
    Rp[i].z = R[i].z;
    R[i].x = Rnew[i].x;
    R[i].y = Rnew[i].y;
    R[i].z = Rnew[i].z;
  }

}



int main( int argn, char * argv[] )
{
  int dbg_print, special, print_energy, swaps;
  double coolingTime, endDelayPeriod, delayTime, heatingTime, recoolTime;
  bool startHeating = false;
  bool startRecooling = false;

  char traj_filename[20];
  strcpy(traj_filename, "trajectories");
  strcat(traj_filename, argv[1]);
  char swaptimes_filename[20];
  strcpy(swaptimes_filename, "swaptimes");
  strcat(swaptimes_filename, argv[1]);
  char readout_filename[20];
  strcpy(readout_filename, "readout");
  strcat(readout_filename, argv[1]);
  initialTime = clock();
                        // extra debug print
                            // time grid
  tmin = 0.0;
  tmax =   .040100389084; // 2 million 0.010;
  dt =     .000000000333; // sec (for micromotion should be at least a 50th of the drive frequency. It's a 100th right now) 

  coolingTime = .00100;
  heatingTime = 0.019100389084;
  delayTime = atof(argv[2])*.0000001;
  endDelayPeriod = coolingTime + delayTime;
  recoolTime = endDelayPeriod + heatingTime;

  double percentDoneTimes[101];
  int i;
  for (i=0;i<101; i++)
  {
    percentDoneTimes[i] = tmax*0.01*(double)i;
  }

  timeStepFactor = 1.0/dt; // number of time steps per second
  
                        // information
  printf( "\n Npart: %d tmax %f\n", Npart, tmax ); 
                        // recording trajectories in file
  fp_traj = fopen( traj_filename, "w" );
  strcat(traj_filename, ".xyz");
  fp_traj_xyz = fopen( traj_filename, "w" );
  fp_traj_swap = fopen( swaptimes_filename, "w" );
  fp_traj_readout = fopen (readout_filename, "w");

                        // initial time
  t = tmin;
                        // initial positions & velocities
  initial_conditions( );
                        // print positions & velocities
  fprintf( stderr,
   "\nInitial positions, velocities & energy-momentum (time %f)\n", t );
  output_x_v();
                        // record positions
//  record_trajectories( ); 
  
  // Do one time step to initialize 
  
                       // loop over time
  int percentDoneCounter = 0;
  while ( t < tmax )
    {
                        // time step -- velocity-verlet
      time_step_verlet( dbg_print );
      t = t + dt;
      if (t >= percentDoneTimes[percentDoneCounter])
      {
        fprintf(stderr, "Percent done: %d\n", (int)((t/tmax)*100.0));
        percentDoneCounter++;
      }
                        // record positions
      record_trajectories( );
//      record_trajectories_xyz( );
      
      rfCount++;
      

      // this is for applying different pulses at different times, for example, to turn off laser, change detuning to something far away
//      laser_counter++;

      // experiment 1: just cool first, then heat continuously and monitor the fluorescence

      if (t < coolingTime)
      {
	detuning = -0.5*gamma_param;
      }
      else if (t > coolingTime && t < endDelayPeriod)
      {
	detuning = -500.5*gamma_param; // laser off
      }
      else if (t > endDelayPeriod && t < recoolTime)
      {
        if (startHeating == false)
        {
          detuning = 5.0*gamma_param;
          startHeating = true;
        } 
      }
      else
      {
        if (startRecooling == false)
        {
          double laser_direction_temp[] = {1.0, 1.0, 1.0};
          double laser_direction_temp_magnitude;
          double laser_direction[] = {0.0, 0.0, 0.0};
          laser_direction_temp_magnitude = sqrt(laser_direction_temp[0]*laser_direction_temp[0] + laser_direction_temp[1]*laser_direction_temp[1] + laser_direction_temp[2]*laser_direction_temp[2]); 
          int j;
          // establish proper sizes for laser_direction and k_laser
          for (j = 0; j<3; j++)
          {
            laser_direction[j] = laser_direction_temp[j]/laser_direction_temp_magnitude;
            recoil_velocity[j] = laser_velocity*laser_direction[j];
            k_laser[j] = (M_PI*2.0/(396.847e-9))*laser_direction[j];
          }
          detuning = 5.0*gamma_param;
          startRecooling = true;
          s_beam = 0.0;
          for (i = 0; i < Npart; i++)
          {
            emit_counter[i] = 0;
            fluoravg[i] = 0;
          }
//	  binCount = 0;
//	  rfCount = 0;

        } 
//	record_trajectories( );
//	rfCount++;
      }
//        record_trajectories_readout( );
//	rfCount++;

      // experiment 2: cool first, start heating, collect fluorescence during readout time! make sure you change tmax!

/*      
      if (t < coolingTime)
      {
	detuning = -0.5*gamma_param;
      }
      else if (t > coolingTime && t < timeOfReadout)
      {
	detuning = 0.5*gamma_param;
      }
      else
      {
      	detuning = 0.0*gamma_param;
	record_trajectories_readout( );
      }
*/

/*
     swaps = (swapCounter / 2);
      if (swaps >= swapsToMelt)
      {
        tmax = t;        
        break;
      }
*/      
    }
  printf("\nAverage energies and fluorescence: \n");
  for (i = 0; i < Npart; i++)
  {
    printf("Average energy for particle %d: %.30f\n", i, totalAvg[i]/((tmax - tmin)/(dt*bin)));
  }
  
  for (i = 0; i < Npart; i++)
  {
    printf("Emit count for particle %d: %.30f\n", i, total_emit_counter[i]);
//    printf("Excitation probability for particle %d: %.30f\n", i, total_emit_counter[i]/total_check_counter[i]);
//    printf("Average lifetime for particle %d: %.30f ns\n", i, 1.0e9*lifetimes[i]/total_emit_counter[i]); // total lifetimes divided by number of emissions (average lifetime)
    printf("Probability of particle %d being in excited state: %.30f\n", i, time_in_excited_state[i]/tmax);
    printf("Probability of particle %d being in ground state: %.30f\n", i, time_in_ground_state[i]/tmax);
    printf("Sum of ground and excited probabilities for particle %d: %.30f\n", i, (time_in_excited_state[i]/tmax) + (time_in_ground_state[i]/tmax) );
    printf("time*(linewidth)*Pexcited for particle %d (should equal emit count): %.30f\n", i, tmax*(gamma_param)*(time_in_excited_state[i]/tmax) ); 
    printf("%.30f\n", gamma_lasers/total_check_counter[i]);   
  }  
    printf("Rho_ee for our s_beam with 0 detuning: %.30f\n", (s_beam*0.5 / (1.0 + s_beam)) ); // weird cus there's no gamma??
  
                        // print positions & velocities
  output_x_v();
                        // close trajectories file
  fclose( fp_traj );
//  fclose( fp_traj_xyz );
  fclose( fp_traj_swap );
  fclose( fp_traj_readout );

  finalTime = clock();
  printf("Time to simulate %d particles for %.15f seconds: %.30f\n", Npart, tmax, (double)(finalTime - initialTime)/CLOCKS_PER_SEC);
  // All done, let's print the results!
  int count;
  printf("swapcount: %d \n", swapCounter);
/*
  for (count = 0; count<swapCounter; count++)
  {
    printf("Ion %d is now in position %d\n", swapCatalog[count][0], swapCatalog[count][1]);
  }
*/

}

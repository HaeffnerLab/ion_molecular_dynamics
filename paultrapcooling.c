#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

/* This code was adapted from another simulation - I can't seem to find the source. */

/********************************************************************************************************************************* 

If you are not familiar with C:

Note: certain values for some variables defined must be calculated. Because values cannot be calculated in the global scope,
the values are set here as NULL pointers and are set to the correct values in later functions such as 'initial_conditions'.

Make sure when doing calculations with these specific variables, they are done as: (*variable) not (variable)

*********************************************************************************************************************************/



/* allows you to use "true" and "false" as booleans */
typedef int bool;
enum { false, true };

#define DBL_MAX 1e37     // largest number a double can handle
#define Npartmax 2000    // Largest number of particles ( Npart < Npartmax )
#define MAXSWAPS 100000  // maximum number of swaps. Cataloging more swaps than this number will lead to a segmentation fault

/* structures that contain particle data */
typedef struct{
  double x, y, z;
} position;

typedef struct{
  double vx, vy, vz;
} velocity;

typedef struct{
  double ax, ay, az;
} acceleration;

int Npart = 5;                                                  // Number of particles
double mass = 6.6551077e-26;                                    // mass of calcium in kg
position R[Npartmax];                                           // positions
velocity V[Npartmax];                                           // velocities
acceleration acc[Npartmax];                                     // accelerations
FILE *fp_traj, *fp_traj_xyz, *fp_traj_swap, *fp_traj_readout;   // files
position Rp[Npartmax], Rnew[Npartmax];                          // previous and new positions (needed for Verlet algorithm)
double dt, tmin, tmax;                                          // time setup
double trap_freq_x, trap_freq_y, trap_freq_z;                   // for the pseudopotential
double boltzmann = 1.3806503e-23;
double hbar = 1.05457148e-34;
double initialTemp = 0.0005;                                    // in Kelvin - the temperature the ions will start at
double decaytime = 0.0000000071;                                // decay time for working transition in Ca -> 7.1ns
double detuning_parameter = -.5;                                // initial detuning of the laser
double detuning;                                                // detuning_parameter * gamma
double s_beam = 1.0;                                            // saturation parameter (arbitrarily set)
double laser_direction_temp[] = {1.0, 1.0, 1.0};                // laser comes in at 45 degree angle - temporary variable
double laser_direction[] = {0.0, 0.0, 0.0};                     // laser direction 
double gamma_param = 140845070;                                 // 1 / (7.1ns)
double laser_velocity = 0.024877240126316254;                   // hbar*k/mass
double theta_initial, phi_initial, theta_recoil, phi_recoil;    // used for creating a random direction
int *decayState = NULL;                                         // array to keep track of the ions' electronic states
double *excitedTimer = NULL;                                    // array to keep track of the ions' internal timer during absorption (prior to emmission)
double *time_in_excited_state = NULL;                           // array to keep track of the ions' time in the excited state
double *time_in_ground_state = NULL;                            // array to keep track of the ions' time in the excited state
double t;                                                       // official simulation time
double initialTime, finalTime;                                  // used for computing how long the simulation runs

/* Variables used to keep track of ion swapping */
int swapCounter = 0;                                            // count of the total number of swaps
int swapsToMelt = 2000;                                         // if implemented (in the main while loop), the simulation will stop when this number of swaps is reached       
int swapCatalog[MAXSWAPS][2];                                   // a 2 column array, each entry records which two ions swapped, ions are labeled: 0, 1, 2, ....
int ionTracking[Npartmax];                                      // keep track of an ion's position relative to others in the chain

/* micromotion - these values were calculated such that the desired trap frequencies are achieved - see mathematica notebook */

double V_1 = 45.55362938;                                       // voltage applied to the Y blades in volts
double V_0 = 41.94377859;                                       // voltage applied to the X blades in volts
double U_0 = 16.39851231;                                       // voltage applied to the endcaps (Z-direction) in volts
double R_0 = 0.0005;                                            // trap radius in the XY plane
double Z_0 = 0.01;                                              // trap length (Z-direction)
double kappa = 0.5;                                             // geometry factor
double trap_drive_omega;                                        // RF drive frequency
double electron_charge = 1.602e-19;                             // charge of an electron in coulombs

/* used for averaging and plotting. */
int binCount;                                                   // a counter to keep track of the number of RF periods
int rfCount;                                                    // a secondary counter for keeping track of time steps within an RF period (this is NOT enforced)
double bin;                                                     // number of sets of RF periods to average over
double numRFsteps;                                              // number of time steps in one RF period to average over (again, nothing enforces this to be exactly the number of time steps within an RF cycle)
double *emit_counter = NULL;                                    // a counter to keep track of the number of photon emissions within an averaging cycle
double *total_emit_counter = NULL;                              // a counter to keep track of the total number of photon emissions
double *velocitiesX = NULL;                                     // sum of velocities in the X direction within an averaging cycle
double *velocitiesY = NULL;                                     // sum of velocities in the X direction within an averaging cycle
double *velocitiesZ = NULL;                                     // sum of velocities in the Z direction within an averaging cycle
double *velocitiesSquaredX = NULL;                              // sum of squares of velocities in the X direction within an averaging cycle
double *velocitiesSquaredY = NULL;                              // sum of squares of velocities in the Y direction within an averaging cycle
double *velocitiesSquaredZ = NULL;                              // sum of squares of velocities in the Z direction within an averaging cycle
double *tempavg = NULL;                                         // average "temperature" within an averaging cycle (defined in units such that temperature = 1 when the Doppler temperature is reached - see below)
double *tempavgx = NULL;                                        // average "temperature" in the X direction
double *tempavgy = NULL;                                        // average "temperature" in the Y direction
double *tempavgz = NULL;                                        // average "temperature" in the Z direction
double *fluoravg = NULL;                                        // average fluorescence                                                
double *totalAvg = NULL;                                        // average of the temperature over the whole simulation (not always a useful value)

bool startReadout = false;                                      // used to signal when "readout" has begun, data is recording in a separate file called "readout" - note that this is purely for convenience, for example:
                                                                // if one wants to separate the heating data from the readout data (this makes data analysis quicker)

/* The following variables require calculation prior to assignment - see note above */
double *lattice_power;                                          // power of the optical lattice (not used currently)
double *damping;                                                // damping parameter
double *recoil_velocity = NULL;                                 // recoil velocity of the ions upon absorption
double *k_laser = NULL;                                         // wavevector of the laser


/* variables for averaging and plottings */
void setup_averaging_variables()
{
    
    int i;
    binCount = 0;
    bin = 100.0; // NUMBER OF SETS OF RF STEPS
    numRFsteps = 100.0; // NUMBER OF RF STEPS
    rfCount = 0;

    /* allocate memory for the arrays */
    emit_counter = malloc(Npart*sizeof(double));
    total_emit_counter = malloc(Npart*sizeof(double));
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

/* output routine */
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

/* record trajectory in .xyz format (for movies!) */
void record_trajectories_xyz( )
{
  int i;
  double x, y, z;

  fprintf(fp_traj_xyz, "%d\n", Npart);
  fprintf(fp_traj_xyz, "Ions in a Paul Trap with a Laser XYZ\n");
  
  for ( i=0 ; i<Npart ; i++ )
    {
      /* scale the positions so that the movements can be noticed in VMD */
      x = R[i].x*1000000.0;
      y = R[i].y*1000000.0;
      z = R[i].z*1000000.0;
      /* This part is specific to 5 ions, this ensures that the 2nd and 5th ion have different colors than the others by assigning them different atom names */
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


/* record trajectory in a file */
void record_trajectories( )
{
  int i;
  double xtemp, ytemp;

  /* add the velocities at every time step */
  for ( i=0 ; i<Npart ; i++ )
  {
    velocitiesX[i] += V[i].vx;
    velocitiesY[i] += V[i].vy;
    velocitiesZ[i] += V[i].vz;
  }

  /* if an RF period is reached - calculated the velocities squared and add them to the previous value */
  if (rfCount == numRFsteps)
  {
    for ( i=0 ; i<Npart ; i++ )
    {

      velocitiesSquaredX[i] += (velocitiesX[i]/numRFsteps)*(velocitiesX[i]/numRFsteps);
      velocitiesSquaredY[i] += (velocitiesY[i]/numRFsteps)*(velocitiesY[i]/numRFsteps);
      velocitiesSquaredZ[i] += (velocitiesZ[i]/numRFsteps)*(velocitiesZ[i]/numRFsteps);
      /* reset the velocities for the next RF cycle */
      velocitiesX[i] = 0.0;
      velocitiesY[i] = 0.0;
      velocitiesZ[i] = 0.0;
    }
 
    rfCount = 0; 
    binCount++; // update the number of RF cycles within this averaging period
//    record_trajectories_xyz(); // record movie every RF cycle!! uncomment this to record movies (this takes a lot of disk space)
  }

  /* if the number of RF cycles is reached - calculate all values to be recorded and reset all variables*/
  if (binCount == bin)
  {
    /* check if it's time to do "readout" */
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
          fprintf(fp_traj_readout, " %.30f ", fluoravg[i]); // readout only cares about the fluorescence
        }
        /* reset variables */
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

void initial_conditions( )
{
  srand(time(NULL)); // random setup
  double initialVelocity;
  int i, j;

  setup_averaging_variables();

  double initialZPositions[] = {-0.00003592197686982309,-0.00001694384644835323, -0.00000000000003984558,0.00001694384636012078,0.00003592197676604132};//{-0.00012400338910260807, -0.00010032968053235605, -0.00008071905149605598, -0.00006305777988679104, -0.00004652649444710114, -0.00003068542210321371, -0.00001524865657422565, -0.00000000000000001103, 0.00001524865657420344, 0.00003068542210319132, 0.00004652649444707858, 0.00006305777988676827, 0.00008071905149603312, 0.00010032968053233323, 0.00012400338910258609};
 

  trap_drive_omega = 2.0*M_PI*30000000.0; // 2.0*M_PI*15000000
  
  lattice_power = malloc(sizeof(double)); // lattice power is off for now  
  *lattice_power = 0.0;
  damping = malloc(sizeof(double));
  *damping = 0.0; //1.0e-19;  //1e-18 is overdamped // damping is off for now
   
//  trap_freq_x = 2.0*M_PI*1000000.0
//  trap_freq_y = 2.0*M_PI*1000000.0
//  trap_freq_z = 2.0*M_PI*100000.0
  
  decayState = malloc(Npart*sizeof(int));
  excitedTimer = malloc(Npart*sizeof(double));
  time_in_excited_state = malloc(Npart*sizeof(double));
  time_in_ground_state = malloc(Npart*sizeof(double));
  
  /* laser-specific */
  
  double laser_direction_temp_magnitude;
  laser_direction_temp_magnitude = sqrt(laser_direction_temp[0]*laser_direction_temp[0] + laser_direction_temp[1]*laser_direction_temp[1] + laser_direction_temp[2]*laser_direction_temp[2]); 
  /* establish proper sizes for laser_direction and k_laser */
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
      
      /* initial positions */
      R[i].x = 0.0;
      R[i].y = 0.0;
      R[i].z = initialZPositions[i];

      /* the ions are intially already labelled in order */
      ionTracking[i] = i;

      decayState[i] = 0;
      excitedTimer[i] = 0.0;
      time_in_excited_state[i] = 0.0;
      time_in_ground_state[i] = 0.0;
      
      /* random velocities normalized such that the initial temperature condition is met */
      theta_initial = acos(2.0*((double)rand()/(double)RAND_MAX)-1.0);
      phi_initial = 2.0*M_PI*((double)rand()/(double)RAND_MAX);
      V[i].vx = initialVelocity*sin(theta_initial)*cos(phi_initial);
      V[i].vy = initialVelocity*sin(theta_initial)*sin(phi_initial);
      V[i].vz = initialVelocity*cos(theta_initial);

      /* previous value (needed in Verlet step) Rp == "R Previous" - Euler step */
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
    
  printf("Initial velocity: %.30f %.30f %.30f\n", V[0].vx, V[0].vy, V[0].vz );  
}


/* all forces calculated here */
void accelerations()
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
  
/*  
  // harmonic pseudopotential
   for ( i=0 ; i<Npart ; i++ )
   {
     Fx[i] = -R[i].x*mass*(trap_freq_x*trap_freq_x);// - (*damping)*V[i].vx;
     Fy[i] = -R[i].y*mass*(trap_freq_y*trap_freq_y);// - (*damping)*V[i].vy;
     Fz[i] = -R[i].z*mass*(trap_freq_z*trap_freq_z);// - (*damping)*V[i].vz + (*lattice_power)*(*k)*sin(R[i].z*(*k));
   }
*/  
  
  /* Coulomb */  
  for ( i=0 ; i<Npart-1 ; i++ )
    {
      for ( j= i+1 ; j<Npart ; j++ )
	{

	  magnitude = sqrt((R[i].x - R[j].x)*(R[i].x - R[j].x) + (R[i].y - R[j].y)*(R[i].y - R[j].y) + (R[i].z - R[j].z)*(R[i].z - R[j].z));
              xij = (R[i].x - R[j].x)/fabs(pow(magnitude, 3.0));
              yij = (R[i].y - R[j].y)/fabs(pow(magnitude, 3.0));
	      zij = (R[i].z - R[j].z)/fabs(pow(magnitude, 3.0));
              
              /* 2.30e-28  is 4*Pi*epsilon_0 */
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
  double doppler;
  double probability;
  double gamma_laser;
  double s;

  doppler = ((k_laser[0])*V[ion].vx + (k_laser[1])*V[ion].vy  + (k_laser[2])*V[ion].vz); 
  s = s_beam/(pow((2.0*(detuning-doppler)/gamma_param), 2.0) + 1.0);
  gamma_laser = gamma_param*s/2.0;
  probability = gamma_laser*dt;
  
  return probability;
}

void catalogIonSwaps()
{
  /* This function finds out if ions have swapped places and catalogs which ion moved where - it does this by comparing the Z positions relative to each other and keeping track of the result */

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
      swapCatalog[swapCounter][0] = ionTrackingTemp[i];
      swapCatalog[swapCounter][1] = ionTracking[i];   
      fprintf(fp_traj_swap, "%.30f %d %d\n", t, swapCatalog[swapCounter][0], swapCatalog[swapCounter][1] );
      swapCounter++; 
    }
    
    ionTracking[i] = ionTrackingTemp[i];    
  }

}

/* Velocity Verlet time step */
time_step_verlet()
{
  int i;
  double emit_rand_temp;
  double excitation_prob; 
                  
  accelerations();           // accelerations at t
                  
  for ( i=0 ; i<Npart ; i++ )
  {
    // check if a particle is already excited or not
    if (decayState[i] == 0)
    { 
      /* if we're here, the particle is not currently excited - check for excitation */
      time_in_ground_state[i] += dt;
      excitation_prob = get_excitation_probability(i);
      if ( (double)rand()/(double)RAND_MAX < excitation_prob )
      {
	/* excite! change decay state */
	decayState[i] = 1;
	/* start the excited timer! */
        excitedTimer[i] += dt;

	// absorb! - (notice the recoil velocity in the Verlet step)
	Rnew[i].x = 2.0*R[i].x - Rp[i].x + recoil_velocity[0]*dt + acc[i].ax*dt*dt;
	Rnew[i].y = 2.0*R[i].y - Rp[i].y + recoil_velocity[1]*dt + acc[i].ay*dt*dt;
	Rnew[i].z = 2.0*R[i].z - Rp[i].z + recoil_velocity[2]*dt + acc[i].az*dt*dt;
   
	V[i].vx = ( Rnew[i].x - Rp[i].x ) / ( 2.0*dt );
	V[i].vy = ( Rnew[i].y - Rp[i].y ) / ( 2.0*dt );
	V[i].vz = ( Rnew[i].z - Rp[i].z ) / ( 2.0*dt );
	 
      }
      else
      {
	/* no excitation evolve normally */
	
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
      /* if we're here, that means the ion is already excited! */
      time_in_excited_state[i] += dt;
      emit_rand_temp = (double)rand()/(double)RAND_MAX;
      if (emit_rand_temp < (excitation_prob + gamma_param*dt) ) // check if it's time to emit for this ion
      {
	/* yes! do a Verlet step that includes the emission */
	if (emit_rand_temp < gamma_param*dt)
        {
	  emit_counter[i]++;
          total_emit_counter[i]++;

  	  phi_recoil = acos(2.0*((double)rand()/(double)RAND_MAX)-1.0);
	  theta_recoil = 2.0*M_PI*((double)rand()/(double)RAND_MAX);

	  Rnew[i].x = 2.0*R[i].x - Rp[i].x + ( (laser_velocity*sin(theta_recoil)*cos(phi_recoil))*dt ) + acc[i].ax*dt*dt;
	  Rnew[i].y = 2.0*R[i].y - Rp[i].y + ( (laser_velocity*sin(theta_recoil)*sin(phi_recoil))*dt ) + acc[i].ay*dt*dt;
	  Rnew[i].z = 2.0*R[i].z - Rp[i].z + ( (laser_velocity*cos(theta_recoil))*dt ) + acc[i].az*dt*dt;

        }
        else
        {
          /* stimulated emission - notice the recoil velocity in the Verlet Step */
	  Rnew[i].x = 2.0*R[i].x - Rp[i].x - recoil_velocity[0]*dt + acc[i].ax*dt*dt;
	  Rnew[i].y = 2.0*R[i].y - Rp[i].y - recoil_velocity[1]*dt + acc[i].ay*dt*dt;
	  Rnew[i].z = 2.0*R[i].z - Rp[i].z - recoil_velocity[2]*dt + acc[i].az*dt*dt;
	}

	V[i].vx = ( Rnew[i].x - Rp[i].x ) / ( 2.0*dt );
	V[i].vy = ( Rnew[i].y - Rp[i].y ) / ( 2.0*dt );
	V[i].vz = ( Rnew[i].z - Rp[i].z ) / ( 2.0*dt );
	
	decayState[i] = 0; // go back to ground state
        excitedTimer[i] = 0.0; // set the timer back to 0
      }
      else
      {
	excitedTimer[i] += dt;
	/* if we're here, that means we did not emit, so just evolve normally */
	Rnew[i].x = 2.0*R[i].x - Rp[i].x + acc[i].ax*dt*dt;
	Rnew[i].y = 2.0*R[i].y - Rp[i].y + acc[i].ay*dt*dt;
	Rnew[i].z = 2.0*R[i].z - Rp[i].z + acc[i].az*dt*dt;
	
	V[i].vx = ( Rnew[i].x - Rp[i].x ) / ( 2.0*dt );
	V[i].vy = ( Rnew[i].y - Rp[i].y ) / ( 2.0*dt );
	V[i].vz = ( Rnew[i].z - Rp[i].z ) / ( 2.0*dt );
	
      }
    }  
  }
   
  /* catalog ion swapping - the proceeding line can be commented out if one does not wish to record swaps */
  catalogIonSwaps();

  /* reset arrays */
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
  
  /* USAGE: "./paultrapcooling FILENUMBER DELAYTIME(in us) "
     example: "./paultrapcooling 0 0 " will create files with the suffix "0" */  

  int swaps;                                                                    // keeps track of the swaps
  double coolingTime, endDelayPeriod, delayTime, heatingTime, recoolTime;       // times throughout the "pulse" sequence
  bool startCooling = false;
  bool startHeatDelay = false;
  bool startHeating = false;
  bool startRecooling = false;

  /* filenames */
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

  tmin = 0.0;
  tmax =   .030700000000; 
  dt =     .000000000333;                                                       // sec (for micromotion should be at least a 50th of the drive frequency. It's a 100th right now) 

  coolingTime = .00100;                                                         // initial cooling time
  heatingTime = .01970;                                                         // heating time
  delayTime = atof(argv[2])*.0000001;
  endDelayPeriod = coolingTime + delayTime;                                     // a delay before heating is implemented in the command line
  recoolTime = endDelayPeriod + heatingTime;

  /* for printing the percentage of the simulation complete */
  double percentDoneTimes[101];
  int i;
  for (i=0;i<101; i++)
  {
    percentDoneTimes[i] = tmax*0.01*(double)i;
  }

  

  printf( "\n Npart: %d tmax %f\n", Npart, tmax ); 
  fp_traj = fopen( traj_filename, "w" );
  strcat(traj_filename, ".xyz");
  fp_traj_xyz = fopen( traj_filename, "w" );
  fp_traj_swap = fopen( swaptimes_filename, "w" );
  fp_traj_readout = fopen (readout_filename, "w");


  t = tmin;
  initial_conditions( );                                                        // set initial conditions
  /* print positions & velocities */
  fprintf( stderr, "\nInitial positions, velocities & energy-momentum (time %f)\n", t );
  output_x_v();

  /* main loop */ 
  int percentDoneCounter = 0;
  while ( t < tmax )
    {

      time_step_verlet();
      t = t + dt;
      if (t >= percentDoneTimes[percentDoneCounter])
      {
        fprintf(stderr, "Percent done: %d\n", (int)((t/tmax)*100.0));
        percentDoneCounter++;
      }

      record_trajectories( );                                                   // record data - averaging is taken care of within the function
      rfCount++;                                                                // update the number of steps within an RF cycle - this number MUST be incremented in order to record properly
      

      if (t < coolingTime)
      {
        if (startCooling == false)
        {
          detuning = -0.5*gamma_param;
          startCooling = true;
        } 
      }
      else if (t > coolingTime && t < endDelayPeriod)
      {
        if (startHeatDelay == false)
        {
          detuning = -1000000000000.0*gamma_param; // laser off
          startHeatDelay = true;
        } 
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
        if (startRecooling == false) // "recooling" refers to an experiment where the "readout" after heating is far enough red-detuned such that the ions cool - one can do an actual readout simply by adjust the detuning here to whatever one wants
        {
          /* readjust the laser - direction*/
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
          /* readjust the laser - detuning/saturation */
          detuning = -0.4*gamma_param;
          s_beam = 1.0;
          startRecooling = true;
          /* if one wants to record data ONLY during this time period uncomment the following: */
/*                                                                                                      <-- THIS LINE                                                                                             
          for (i = 0; i < Npart; i++)
          {
            emit_counter[i] = 0;
            fluoravg[i] = 0;
          }
	  binCount = 0;
	  rfCount = 0;
*/                                                                                                     //<-- THIS LINE
        } 
//	record_trajectories( );                                                                         <-- THIS LINE - note: comment out the same line above (in the while loop)
//	rfCount++;                                                                                      <-- THIS LINE - note: comment out the same line above (in the while loop)
      }

/* cut the simulation short if a desired number of swaps is reached */
/*
     swaps = (swapCounter / 2);
      if (swaps >= swapsToMelt)
      {
        tmax = t;        
        break;
      }
*/      
    }

  /* print relevent info */  

  printf("\nAverage energies and fluorescence: \n");
  for (i = 0; i < Npart; i++)
  {
    printf("Average energy for particle %d: %.30f\n", i, totalAvg[i]/((tmax - tmin)/(dt*bin)));
  }
  
  for (i = 0; i < Npart; i++)
  {
    printf("Emit count for particle %d: %.30f\n", i, total_emit_counter[i]);
    printf("Probability of particle %d being in excited state: %.30f\n", i, time_in_excited_state[i]/tmax);
    printf("Probability of particle %d being in ground state: %.30f\n", i, time_in_ground_state[i]/tmax);
    printf("Sum of ground and excited probabilities for particle %d: %.30f\n", i, (time_in_excited_state[i]/tmax) + (time_in_ground_state[i]/tmax) );
    printf("time*(linewidth)*Pexcited for particle %d (should equal emit count): %.30f\n", i, tmax*(gamma_param)*(time_in_excited_state[i]/tmax) ); 
  }  
    printf("Rho_ee for our s_beam with 0 detuning: %.30f\n", (s_beam*0.5 / (1.0 + s_beam)) ); // weird cus there's no gamma??
  

  output_x_v(); // print final pos/vel

  fclose( fp_traj );
  fclose( fp_traj_xyz );
  fclose( fp_traj_swap );
  fclose( fp_traj_readout );

  finalTime = clock();
  printf("Time to simulate %d particles for %.15f seconds: %.30f\n", Npart, tmax, (double)(finalTime - initialTime)/CLOCKS_PER_SEC);
  // All done, let's print the results!
  int count;
  printf("swapcount: %d \n", swapCounter);

}

//  dsmceq - Dilute gas simulation using DSMC algorithm
//  This version illustrates the approach to equilibrium

#include "NumMeth.h"
#include "SortList.h"

using namespace std;

double rand( long& seed );
int colider( Matrix& v, Matrix& crmax, double tau, long& seed,
			  Matrix& selxtra, double coeff, SortList& sD );
void sorter( Matrix& x, double L, SortList &sD );

int main() {

  // Initialize constants  (particle mass, diameter, etc.)
  const double pi = 3.141592654;
  const double boltz = 1.3806e-23;    // Boltzmann's constant (J/K)
  double mass = 6.63e-26;       // Mass of argon atom (kg)
  double diam = 3.66e-10;       // Effective diameter of argon atom (m)
  double T = 273;               // Temperature (K)
  double density = 1.78;        // Density of argon at STP (kg/m^3)
  double L = 1e-6;              // System size is one micron

  // My additions to compute number density and theoretical mean free path
  double n = density / mass;    // Number density (particles / m^3)
  double sigma = pi * diam * diam; // Collision cross-section
  double lambda_th = 1.0 / (sqrt(2.0) * n * sigma); // mean-free path
  // cout << "Theoretical mean free path: " << lambda_th << " m" << endl;

  // User input for nPart
  cout << "Enter number of simulation particles: "; 
  int npart; cin >> npart;
  double eff_num = density/mass*L*L*L/npart;
  cout << "Each particle represents " << eff_num << " atoms" << endl;

  // Assign random positions and velocities to particles
  long seed = 1;       // Initial seed for rand (DO NOT USE ZERO)
  double v_init = sqrt(3.0*boltz*T/mass);    // Initial speed
  Matrix x(npart), y(npart), z(npart), v(npart,3), vxI(npart);
  Matrix xU(npart), yU(npart), zU(npart); // “U” for “unwrapped”
  int i;
  for (i = 1; i <= npart; i++) {
    double rx = rand(seed);
    double ry = rand(seed);
    double rz = rand(seed);
    
    // Original x array (with periodic)
    x(i) = L * rx;
    y(i) = L * rx;
    z(i) = L * rx;
    
    // Unwrapped arrays (NO fmod here)
    xU(i) = L * rx;
    yU(i) = L * ry;
    zU(i) = L * rz;

	// My addition for adjusting isotropic initial condition
    int plusMinus1 = (1 - 2*((int)(2*rand(seed))));
    int plusMinus2 = (1 - 2*((int)(2*rand(seed))));
    int plusMinus3 = (1 - 2*((int)(2*rand(seed))));
    double factor = 1.0 / sqrt(3.0);  // Scale factor so that the total speed is v_init
    v(i,1) = plusMinus1 * factor * v_init;
    v(i,2) = plusMinus2 * factor * v_init;
    v(i,3) = plusMinus3 * factor * v_init;
    
    // Record the initial x velocity
    vxI(i) = v(i,1);
  }

  // Record inital particle speeds
  Matrix vmagI(npart);
  for( i=1; i<=npart; i++ )
	vmagI(i) = sqrt( v(i,1)*v(i,1) + v(i,2)*v(i,2) + v(i,3)*v(i,3) );

  // Initialize variables used for evaluating collisions
  int ncell = 50; //15                       // Number of cells
  double tau = 0.2*(L/ncell)/v_init;    // Set timestep tau
  Matrix vrmax(ncell), selxtra(ncell);
  vrmax.set(3*v_init);    // Estimated max rel. speed
  selxtra.set(0.0);       // Used by routine "colider"
  double coeff = 0.5*eff_num*pi*diam*diam*tau/(L*L*L/ncell);
  int coltot = 0;         // Count total collisions
  double totalDist = 0.0;

  // Declare object for lists used in sorting
  SortList sortData(ncell,npart);  

  // Loop for the desired number of time steps
  cout << "Enter total number of time steps: ";
  int istep, nstep; cin >> nstep;
  for( istep = 1; istep<=nstep; istep++ ) {
	
    // Accumulate total distance traveled
	for (i = 1; i <= npart; i++) {
	double speed_i = sqrt(
		  v(i,1)*v(i,1)
		+ v(i,2)*v(i,2)
		+ v(i,3)*v(i,3)
	);
	totalDist += speed_i * tau;
	}
	
	// Move all the particles ballistically
	for( i=1; i<=npart; i++ ) {
      x(i) += v(i,1)*tau;          // Update x position of particle
      x(i) = fmod(x(i)+L,L);       // Periodic boundary conditions
      
      y(i) += v(i,2)*tau;
	  y(i)  = fmod(y(i) + L, L);

	  z(i) += v(i,3)*tau;
	  z(i)  = fmod(z(i) + L, L);
      
      // "Unwrapped" coords (no BCs):
	  xU(i) += v(i,1)*tau;
	  yU(i) += v(i,2)*tau;
	  zU(i) += v(i,3)*tau;
	}
    // Sort the particles into cells
    sorter(x,L,sortData);
  
    // Evaluate collisions among the particles
    int col = colider(v,vrmax,tau,seed,selxtra,coeff,sortData);
    coltot += col;	// Increment collision count
    
    // Periodically display the current progress
    if( (istep%10) < 1 )
      cout << "Done " << istep << " of " << nstep << " steps; " << 
	         coltot << " collisions" << endl;
  }

  // Record final particle speeds
  Matrix vmagF(npart), vxF(npart);
  for (i = 1; i <= npart; i++) {
    vmagF(i) = sqrt( v(i,1)*v(i,1) + v(i,2)*v(i,2) + v(i,3)*v(i,3) );
    vxF(i)   = v(i,1);  // Record final x velocity
  }

  // Print out the plotting variables: vmagI, vmagF, & vxI, vxF
  ofstream vmagIOut("vmagI.txt"), vmagFOut("vmagF.txt");
  for( i=1; i<=npart; i++ ) {
    vmagIOut << vmagI(i) << endl;
    vmagFOut << vmagF(i) << endl;
  }
  /// My addition print out Vxs fro initial and final terms
  ofstream vxIOut("vxI.txt"), vxFOut("vxF.txt");
	for (i = 1; i <= npart; i++) {
    vxIOut << vxI(i) << endl;
    vxFOut << vxF(i) << endl;
  }
  
  double lambda_sim = totalDist / (2*coltot);
  cout << "Total collisions = " << coltot << endl;
  cout << "Total distance traveled (3D) = " << totalDist << endl;
  cout << "Mean free path from simulation = " << lambda_sim << " m" << endl;
  cout << "Theoretical mean free path: " << lambda_th << " m" << endl;

  return 0;
}
/***** To plot in MATLAB; use the script below ********************
load vmagI.txt; load vmagF.txt;
%* Plot the histogram of the initial speed distribution
vbin = 50:100:1050;    % Bins for histogram
hist(vmagI,vbin);  title('Initial speed distribution');
xlabel('Speed (m/s)');  ylabel('Number');
%* Plot the histogram of the final speed distribution
figure(2); clf;
hist(vmagF,vbin);  
title(sprintf('Final speed distribution'));
xlabel('Speed (m/s)');  ylabel('Number');
******************************************************************/
   

/*
  Josh Carter, 2013

  Demonstration of 'photodynamics' code.

  Call code as photodynam <input_file> <report_times> [> <output_file>].

  <input_file> is file of initial coordinates and properties in 
  following format:

  <N> <time0>
  <step_size> <orbit_error>

  <mass_1> <mass_2> ... <mass_N>
  <radius_1> <radius_2> ... <radius_N>
  <flux_1> <flux_2> ... <flux_N>
  <u1_1> <u1_2> ... <u1_N>
  <u2_1> <u2_2> ... <u2_N>

  <a_1> <e_1> <i_1> <o_1> <l_1> <m_1>
  ...
  <a_(N-1)> <e_(N-1)> <i_(N-1)> <o_(N-1)> <l_(N-1)> <m_(N-1)>

  where the Keplerian coordinates 
  (a = semimajor axis, e = eccentricity, i = inclination,
  o = argument periapse, l = nodal longitude, m = mean anomaly) are the 
  N-1 Jacobian coordinates associated with the masses as ordered above. 
  Angles are assumed to be in radians. The observer is along the positive z axis.
  Rotations are performed according to Murray and Dermott.

  For example, for Kepler-16, kepler16_pd_input.txt:
  
  3 212.12316
  0.01 1e-16

  0.00020335520 5.977884E-05    9.320397E-08
  0.00301596700 0.00104964500   0.00035941463
  0.98474961000	0.01525038700	0.00000000000
  0.65139908000	0.2		0.0
  0.00587581200	0.3		0.0	

  2.240546E-01 1.595442E-01 1.576745E+00 4.598385E+00 0.000000E+00 3.296652E+00
  7.040813E-01 7.893413E-03 1.571379E+00 -5.374484E-01 -8.486496E-06 2.393066E+00

  <report_times> is a list of times to report the outputs.
  First line is a space-separated list of single character-defined 
  output fields according to:
  
  t = time
  F = flux
  a = semi-major axes
  e = eccentricities
  i = sky-plane inclinations
  o = arguments of periapse
  l = nodal longitudes
  m = mean anomalies
  K = full keplerian osculating elements
  x = barycentric, light-time corrected coordinates
  v = barycentric, light-time corrected velocities
  M = masses
  E = fractional energy change from t0
  L = fraction Lz change from t0

  For example, the first line could be

  t F E

  and the output would have three columns of time flux and 
  fractional energy loss.
  
  Output is written to standard out.
  
*/

#include "n_body_state.h"
#include "n_body_lc.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include <stdlib.h>

using namespace std;

int main(int argc, char * argv[]) {

  ifstream in;
  string line;
  int N = 0;
  double t0,maxh,orbit_error;

  if (argc < 3) {
    cerr << "photodynam <input_file> <report_times>"<< endl;
  }

  in.open(argv[1],ios_base::in);

  char input[100];
  if (in.is_open()) {
    in >> input;
    N = atoi(input);
    in >> input;
    t0 = atof(input);
    in >> input;
    maxh = atof(input);
    in >> input;
    orbit_error = atof(input);
  } else { cerr << "File: " << argv[1] << " not found." << endl; exit(1); }
  
  double mass[N],radii[N],flux[N],u1[N],u2[N];

  for (int i = 0; i < N; i++) {
    in >> input;
    mass[i] = atof(input);
  }

  for (int i = 0; i < N; i++) {
    in >> input;
    
    radii[i] = atof(input);
  }

  for (int i = 0; i < N; i++) {
    in >> input;
    flux[i] = atof(input);
  }

  for (int i = 0; i < N; i++) {
    in >> input;
    u1[i] = atof(input);
  }

  for (int i = 0; i < N; i++) {
    in >> input;
    u2[i] = atof(input);
  }

  double mj[N-1],a[N-1],e[N-1],inc[N-1],om[N-1],ln[N-1],ma[N-1];

  for (int i = 0; i < N-1; i++) {
    in >> input;
    a[i] = atof(input);
    in >> input;
    e[i] = atof(input);
    in >> input;
    inc[i] = atof(input);
    in >> input;
    om[i] = atof(input);
    in >> input;
    ln[i] = atof(input);
    in >> input;
    ma[i] = atof(input);
  }

  in.close();
  
  NBodyState state(mass,a,e,inc,om,ln,ma,N,t0);

  state.kep_elements(mj,a,e,inc,om,ln,ma);

  cout.setf(ios_base::scientific);

  cout << "#Photodynamics, Josh Carter 2013\n#" << endl;
  cout << "#Input file is " << argv[1] << "." << endl;
  cout << "#Orbit Error tolerance is " << orbit_error << ".\n";
  cout << "#Step size is " << maxh << ".\n#\n";
  cout << "#" << N << " bodies with following properties: " << endl;
  cout <<"#Index\t|\tMass\t\tRadius\t\tFlux\t\tu1\t\tu2" << endl;
  for (int i = 0; i < N; i++) {
    cout << "# " << i << "\t|\t" << mass[i] << "\t" << radii[i] << "\t" << flux[i] << "\t"
	      << u1[i] << "\t" << u2[i] << endl;
  } 
  cout << "#\n#At Time = " << t0 << endl;
  cout << "#Jacobian coordinates:" << endl;
  cout << "#Index\t|\tMass\t\tx\t\ty\t\tz\t\tv_x\t\tv_y\t\tv_z" << endl;
  for (int i= 1; i < N; i++) {
    cout << "# " <<  i << "\t|\t" <<mj[i-1] <<"\t"<< state.X_J(i) << "\t" << state.Y_J(i) << "\t" << state.Z_J(i) << "\t" <<
      state.V_X_J(i) << "\t" << state.V_Y_J(i) << "\t" << state.V_Z_J(i) << endl;
  }
  cout << "#Keplerian Elements of Jacobian coordinates:" << endl;
  cout << "#Index\t|\tMass\t\ta\t\tPeriod\t\tEcc\t\tInc\t\tomega\t\tOmega\t\tMean Anom." << endl;
  for (int i= 0; i < N-1; i++) {
    cout << "# " << i << "\t|\t" << mj[i] << "\t" << a[i] << "\t" << 2*M_PI*sqrt(a[i]*a[i]*a[i]/mj[i]) << "\t" << e[i] << "\t" << inc[i]*180/M_PI << "\t" << om[i]*180/M_PI << "\t" << ln[i]*180/M_PI << "\t" << ma[i]*180/M_PI << endl;
  }

  cout << "#Barycentric coordinates:" << endl;
  cout << "#Index\t|\tMass\t\tx\t\ty\t\tz\t\tv_x\t\tv_y\t\tv_z" << endl;
  for (int i= 0; i < N; i++) {
    cout <<"# "<< i << "\t|\t" <<state.getMass(i) <<"\t"<< state.X_B(i) << "\t" << state.Y_B(i) << "\t" << state.Z_B(i) << "\t" <<
      state.V_X_B(i) << "\t" << state.V_Y_B(i) << "\t" << state.V_Z_B(i) << endl;
  }

  cout << "#Barycentric coordinates (LTE):" << endl;
  cout << "#Index\t|\tMass\t\tx\t\ty\t\tz" << endl;
  for (int i= 0; i < N; i++) {
    cout << "# "<< i << "\t|\t" <<state.getMass(i) <<"\t"<< state.X_LT(i) << "\t" << state.Y_LT(i) << "\t" << state.Z_LT(i) << endl;
  }

  double e0 = state.getE();
  double lx0,ly0,lz0,lx,ly,lz; state.getL(&lx0,&ly0,&lz0);

  in.open(argv[2],ios_base::in);

  double t;
  cout.setf(ios_base::scientific);

  if (in.is_open()) {
    // read first line and parse output format
    string str;
    getline(in,str);
    std::istringstream iss(str);
    std::vector<std::string> words;

    while (std::getline(iss, str, ' '))
      words.push_back(str);

    cout << "#\n#Begin output from times in " << argv[2] << "." << endl;
    cout << "#";
    
    for (int i = 0; i < words.size(); i++) {
      switch ((words[i].c_str())[0]) {
      case 't':cout << "Time" << "\t\t"; break;
      case 'F':cout << "Flux" << "\t\t"; break;
      case 'x':
	for (int j = 0; j < N; j++) cout << "x_" << j  << "\t\t" << "y_" << j << "\t\t" << "z_" << j << "\t\t";
	break;
      case 'v':
	for (int j = 0; j < N; j++) cout << "vx_" << j << "\t\t" << "vy_" << j << "\t\t" << "vz_" << j << "\t\t";
	break;
      case 'E': cout << "Delta_E/E0" << "\t"; break;
      case 'L': cout << "Delta_Lz/Lz0" << "\t"; break;
      case 'a':
	for (int j = 0; j < N-1; j++) cout << "a_" << j+1 << "\t\t";
	break;
      case 'e':
	state.kep_elements(mj,a,e,inc,om,ln,ma);
	for (int j = 0; j < N-1; j++) cout << "e_" << j+1 << "\t\t";
	break;
      case 'i':
	state.kep_elements(mj,a,e,inc,om,ln,ma);
	for (int j = 0; j < N-1; j++) cout << "inc_" << j+1 << " (deg)" << "\t";
	break;
      case 'o':
	state.kep_elements(mj,a,e,inc,om,ln,ma);
	for (int j = 0; j < N-1; j++) cout << "omega_" << j+1 << " (deg)" << "\t";
	break;
      case 'l':
	
	for (int j = 0; j < N-1; j++) cout << "long_" << j+1 << " (deg)" << "\t";
	break;
      case 'm':
	
	for (int j = 0; j < N-1; j++) cout << "meanan_" << j+1 << " (deg)" << "\t";
	break;
      case 'M':
	for (int j = 0; j < N; j++) cout << "mass_" << j+1 << "\t\t";
	break;
      case 'K':
	
	for (int j = 0; j < N-1; j++)
	  cout << "a_" << j+1 << "\t\t" << "e_" << j+1 << "\t\t" << "inc_" << j+1 << " (deg)" << "\t" 
	       << "omega_" << j+1 << " (deg)" << "\t" << "long_" << j+1 << " (deg)" << "\t" << "meanan_" << j+1 << " (deg)" << "\t";
	break;
      default:
	  cout << "????" << "\t\t" << endl;
      }
    }
    cout << endl;
    // output #column format
    while (in >> input) {
      t = atof(input);
      state(t,maxh,orbit_error,1e-16);
      for (int i = 0; i < words.size(); i++) {
	switch ((words[i].c_str())[0]) {
	case 't':cout << t << "\t"; break;
	case 'F':cout << occultn(state.getBaryLT(),radii,u1,u2,flux,N) << "\t"; break;
	case 'x':
	  for (int j = 0; j < N; j++) cout << state.X_LT(j) << "\t" << state.Y_LT(j) << "\t" << state.Z_LT(j) << "\t";
	  break;
	case 'v':
	  for (int j = 0; j < N; j++) cout << state.V_X_LT(j) << "\t" << state.V_Y_LT(j) << "\t" << state.V_Z_LT(j) << "\t";
	  break;
	case 'E': cout << (state.getE()-e0)/e0 << "\t"; break;
	case 'L': 
	  state.getL(&lx,&ly,&lz);
	  cout << (lz-lz0)/lz0 << "\t"; break;
	case 'a':
	  state.kep_elements(mj,a,e,inc,om,ln,ma);
	  for (int j = 0; j < N-1; j++) cout << a[j] << "\t";
	  break;
	case 'e':
	  state.kep_elements(mj,a,e,inc,om,ln,ma);
	  for (int j = 0; j < N-1; j++) cout << e[j] << "\t";
	  break;
	case 'i':
	  state.kep_elements(mj,a,e,inc,om,ln,ma);
	  for (int j = 0; j < N-1; j++) cout << inc[j]*180/M_PI << "\t";
	  break;
	case 'o':
	  state.kep_elements(mj,a,e,inc,om,ln,ma);
	  for (int j = 0; j < N-1; j++) cout << om[j]*180/M_PI << "\t";
	  break;
	case 'l':
	  state.kep_elements(mj,a,e,inc,om,ln,ma);
	  for (int j = 0; j < N-1; j++) cout << ln[j]*180/M_PI << "\t";
	  break;
	case 'm':
	  state.kep_elements(mj,a,e,inc,om,ln,ma);
	  for (int j = 0; j < N-1; j++) cout << ma[j]*180/M_PI << "\t";
	  break;
	case 'M':
	  for (int j = 0; j < N; j++) cout << state.getMass(j) << "\t";
	  break;
	case 'K':
	  state.kep_elements(mj,a,e,inc,om,ln,ma);
	  for (int j = 0; j < N-1; j++)
	    cout << a[j] << "\t" << e[j] << "\t" << inc[j]*180/M_PI << "\t" 
		 << om[j]*180/M_PI << "\t" << ln[j]*180/M_PI << "\t" << ma[j]*180/M_PI << "\t";
	  break;
	default:
	  cout << "????" << "\t\t";
	}
      }
	  
      cout << endl;
    } 
    
  } else { cerr << "File: " << argv[2] << " not found." << endl; exit(1); }

  state.getL(&lx,&ly,&lz);
  state.kep_elements(mj,a,e,inc,om,ln,ma);

  cout << "#End of output\n#\n";
  cout << "#Final Time = " << state.getTime()<< endl;
  cout << "#Fractional Energy Change: " << (e0-state.getE())/e0 << endl;
  cout << "#Fractional change in Lz: " << (lz0-lz)/lz0 << endl;
  cout << "#Keplerian Elements of Jacobian coordinates:" << endl;
  cout << "#Index\t|\tMass\t\ta\t\tPeriod\t\tEcc\t\tInc\t\tomega\t\tOmega\t\tMean Anom." << endl;
  for (int i= 0; i < N-1; i++) {
    cout << "# "<<i << "\t|\t" << mj[i] << "\t" << a[i] << "\t" << 2*M_PI*sqrt(a[i]*a[i]*a[i]/mj[i]) << "\t" << e[i] << "\t" << inc[i]*180/M_PI << "\t" << om[i]*180/M_PI << "\t" << ln[i]*180/M_PI << "\t" << ma[i]*180/M_PI << endl;
  }

  return 0;
}

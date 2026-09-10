/** 3D simulation of bubble migration in wedge channel. The geometry and fluid parameters 
are similar to the drop tower experiments, reported in Jenson et al. (2014 IJMF). */

#define DROP_TOWER 1

#include "ibm/src/ibm-gcm.h"
#include "ibm/src/my-centered.h"
#include "ibm/src/ibm-gcm-events.h"
#include "ibm/src/my-two-phase.h"
#include "ibm/src/my-tension.h"
#include "tag.h"
#include "view.h"
#include "iso3D.h"
#include "lambda2.h"
#include "navier-stokes/perfs.h"

#define LARGE 1e36
#define f_cut 1.e-6

int max_level = 6;
int min_level = 3;
double L = 1.6e-3;
double t_out = 1e-4;
double t_end = 1.001e-4;

/** Input physical parameters (all in SI unit) 
 */
double rhod   = 50.;
double rhof   = 1000.;
double mud    = 0.0005;
double muf    = 0.01;
double sigma  = 0.005;
double Uf     = 1.0;
double H      = 1.e-4;
double femax  = 0.001;
double uemax  = 0.001;
double maxruntime = 1442;
double z0     = 3e-3;
double x0     = 0.8e-3;
double R0     = 2.06e-3;
double slope  = 0.5; 
double init_film  = 0.1e-3; 
double vmax   = 1.0; 
double pmax   = 1.0;
double grav   = 0;
double init_grav = -9.81;
double grav_plus = -9.81;
double t_grav = 0.5; 

/** Dependent variables */
scalar omega_y[];

double DissRate2, DissWork2, DissWork1, DissRate1;

/** Boundary conditions */

/** top: embedded no-slip wall 
*/
u.n[immersed] = dirichlet(0.);
u.t[immersed] = dirichlet(0.);
u.r[immersed] = dirichlet(0.);

/* bottom: symmetric bc
*/
u.t[bottom] = neumann(0);
u.n[bottom] = dirichlet(0);
p[bottom]   = neumann(0);
pf[bottom]  = neumann(0);
f[bottom]   = neumann(0);

/* bottom: no-slip wall
*/
u.t[back] = dirichlet(0);
u.n[back] = dirichlet(0);
f[back]   = dirichlet(0);

/* left: symmetric bc
*/
u.t[left]   = neumann(0);
u.n[left]   = dirichlet(0);
p[left]     = neumann(0);
pf[left]    = neumann(0);
f[left]     = neumann(0); 

//right: outflow
u.t[right]   = neumann(0);
u.n[right]   = neumann(0);
p[right]     = dirichlet(0);
pf[right]    = dirichlet(0);
f[right]     = dirichlet(0);

int main(int argc, char * argv[])
{
  if (argc > 1)
    max_level = atoi (argv[1]);
  if (argc > 2)
    L         = atof (argv[2]);
  if (argc > 3)
    t_out     = atof (argv[3]);
  if (argc > 4)
    t_end     = atof (argv[4]);
  if (argc > 5)
    rhod      = atof (argv[5]);
  if (argc > 6)
    rhof      = atof (argv[6]);
  if (argc > 7)
    mud       = atof (argv[7]);
  if (argc > 8)
    muf       = atof (argv[8]);
  if (argc > 9)
    sigma     = atof (argv[9]);
  if (argc > 10)
    Uf        = atof (argv[10]);
  if (argc > 11)
    H         = atof (argv[11]);
  if (argc > 12)
    femax     = atof (argv[12]);
  if (argc > 13)
    uemax     = atof (argv[13]);
  if (argc > 14)
    R0        = atof (argv[14]);
  if (argc > 15)
    z0        = atof (argv[15]);
  if (argc > 16)
    slope     = atof (argv[16]);
  if (argc > 17)
    x0        = atof (argv[17]);
  if (argc > 18)
    init_film = atof (argv[18]);
  if (argc > 19)
    vmax      = atof (argv[19]);
  if (argc > 20)
    pmax      = atof (argv[20]);
  if (argc > 21)
    init_grav = atof (argv[21]);
  if (argc > 22)
    grav_plus = atof (argv[22]);
  if (argc > 23)
    t_grav    = atof (argv[23]);
  if (argc > 24)
    min_level = atoi (argv[24]);

  init_grid (1 << (max_level - 2));
  size (L);

  rho1 = rhod, rho2 = rhof;
  mu1  = mud,  mu2  = muf;
  f.sigma = sigma;
 
  TOLERANCE = 1.e-6;
  DT = 1.e-4;

  run();
}



/** Initial condition */
event init (i = 0) {

  if (!restore (file = "dump")) {
    scalar ibm1[], f1[];
    face vector ibmf1[];

    fprintf(stderr, "starting initial refinement\n");
    astats ad;
    int count = 0;
    do {
      solid (ibm1, ibmf1, H-slope*(z)-y);
      fraction(f1, -sq(x-x0)-sq(z-z0)-sq(y) +sq(R0));
      ad = adapt_wavelet ({ibm1,f1}, (double[]){1e-3,femax}, maxlevel=max_level, minlevel=min_level);
      count++;
    fprintf(stderr, "%d refine = %d coarsened = %d\n", count, ad.nf, ad.nc);
    } while ((ad.nc || ad.nf) && count < 50); // adapt until no more cells are
                                              // coarsened or refined
    fprintf(stderr, "initial refinement done\n");
    // Initial mesh refinement to specify liquid volume fraction for the bubble.
    solid (ibm, ibmf, H-slope*(z)-y);
    fraction(f, -sq(x-x0)-sq(z-z0)-sq(y) +sq(R0));

    boundary({ibm, f});
  }
  // restart from dump, but we need to reinitialize the solid volume fraction. 
  else{ 
    solid (ibm, ibmf, H-slope*(z)-y);
    event("update_metric");
    boundary(all);
  }
  // Set culumulative dissipation work for each phase to zero
  DissWork2 = 0.; 
  DissWork1 = 0.; 
}

/** Interfacial force. 
Two interfacial forces are added. The first one is the artificial van der Waals 
force, which will push the interface away from the embedeed wall. The parameters
to control the potential magnitude, the minimum distance between the interface 
and the wall, are specified in the curvature.h. The VDW force is normal to the 
wall and is inversely proportial to the distance between the interface and the 
wedge wall to the power 3. 

The second force is the gravity force in z direction (away from the wedge vertex).
We apply gravity only for a specified period of time (t<t_grav), and will turn 
it off for microgravity phase. 
*/
event acceleration (i++)
{
  // van der Waals force
  scalar phi = f.phi;

  coord G = {0.,0.,0.}, Z = {0.,0.,0.};
  G.y = 1./sqrt(1.+slope*slope), G.z = slope/sqrt(1.+slope*slope);
  Z.y = H;

  if (phi.i)
    potential (f, phi, G, Z, add = true);
  else {
    phi = new scalar;
    potential (f, phi, G, Z, add = false);
    f.phi = phi;
  }
  
  // gravity force
  coord G1 = {0.,0.,0.}, Z1 = {0.,0.,0.};
  
   if (t >= 0 && t <= 1) {
    grav = init_grav*erf(t);
   } else if (t > 1) {
    grav = grav_plus;
   } else {
    grav = init_grav;
   }

  G1.z = grav * (rho2 - rho1);

  if (phi.i)
    position (f, phi, G1, Z1, add = true);
  else {
    phi = new scalar;
    position (f, phi, G1, Z1, add = false);
    f.phi = phi;
  }
}

/** Calculate viscous dissipation rate in gas and liquid phases. 
*/
int dissipation_rate (double* rates)
{
  double DissRate2 = 0.0;
  double DissRate1 = 0.0;
  foreach (reduction (+:DissRate2) reduction (+:DissRate1)) {

    /**
    Gradient of the velocity*/

    double dudx = (u.x[1]     - u.x[-1]    )/(2.*Delta);
    double dudy = (u.x[0,1]   - u.x[0,-1]  )/(2.*Delta);
    double dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
    double dvdx = (u.y[1]     - u.y[-1]    )/(2.*Delta);
    double dvdy = (u.y[0,1]   - u.y[0,-1]  )/(2.*Delta);
    double dvdz = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
    double dwdx = (u.z[1]     - u.z[-1]    )/(2.*Delta);
    double dwdy = (u.z[0,1]   - u.z[0,-1]  )/(2.*Delta);
    double dwdz = (u.z[0,0,1] - u.z[0,0,-1])/(2.*Delta);
    
    double sqterm = 2.*sq(dudx) + 2.*sq(dvdy) + 2.*sq(dwdz) 
                  + sq(dvdx + dudy) + sq(dwdy+dvdz) + sq(dudz+dwdx);
     
    // We only calculate dissipation in non-solid cells
    if ( ibm[] >= 1. ) {
      DissRate2 += (1 - f[])*dv()*sqterm; //water
      DissRate1 +=      f[] *dv()*sqterm; //air
    }
  }
  rates[0] = mu2*DissRate2;
  rates[1] = mu1*DissRate1;
  return 0;
}


/** Calculate dissipation. */
event energy (i++) {
  double rates[2];
  dissipation_rate(rates);
  DissRate2 = rates[0];
  DissRate1 = rates[1];

  DissWork2 += DissRate2*dt;
  DissWork1 += DissRate1*dt;

}

/** Log for temporal outputs */
event logfile (i+=10)
{
  scalar posx[],posy[],posz[];
  position (f,posx,{1,0,0});
  position (f,posy,{0,1,0});
  position (f,posz,{0,0,1});

  double vol=0.,xd=0.,ud=0.,zd=0.,wd=0.,ke2=0.,ke1=0.;

  foreach(reduction(+:vol)
          reduction(+:ud ) reduction(+:xd )
          reduction(+:wd ) reduction(+:zd )
          reduction(+:ke1) reduction(+:ke2)
    ) {
    double dv1 = f[]*dv(), dv2 = dv()*(1. - f[]);

    /** Volume of bubble*/
    vol += dv1;

    double ke = 0.;
    /** kinetic energy of bubble and ambient fluid*/
    foreach_dimension() {
      ke += sq(u.x[]);
    }
    ke1 += dv1*ke;
    ke2 += dv2*ke;

    /** mean velocity of bubble*/
    ud += dv1*u.x[];
    wd += dv1*u.z[];

    /** centroid of bubble */
    xd += dv1*x;
    zd += dv1*z;
  }
  xd /= vol;
  zd /= vol;
  ud /= vol;
  wd /= vol;

  ke1 = ke1*rho1/2.;
  ke2 = ke2*rho2/2.;

  double area = interface_area(f);

  /**
  Measure characteristic length for the bubble. Small bubble fragments may have been prouced. 
  So we tag the gas cells connected and measure the volumes of all bubbles, and we will only 
  consider the largest bubble. The rest are likely fragements.  */

  //Create tags.
  scalar m[];
  foreach()
    m[] = f[] > f_cut;
  int n = tag (m);

  //Calculate the volume of each tagged structure.
  double vol_2[n];
  for (int j = 0; j < n; j++)
    vol_2[j] = 0.0;

  foreach()
    if (m[] > 0){
    int j = m[] - 1;
    vol_2[j] += dv() * f[];
  }
  
  //Find the tag ID for the largest bubble.
  int largest_index = 0;
  double max_volume = 0.0;
  for (int j = 0; j < n; j++) {
    if (vol_2[j] > max_volume) {
      max_volume = vol_2[j];
      largest_index = j;
    }
  }

  //Measure characteristic lengths
  double max_x = -LARGE, max_y = -LARGE, max_z= -LARGE, z_at_x_max = -LARGE;
  double min_x =  LARGE, min_z =  LARGE; 
  foreach(reduction(max:max_x) reduction(min:min_x) reduction(max:max_y) reduction(max:max_z) reduction(min:min_z) reduction(max:z_at_x_max))
    if (m[] - 1 == largest_index && f[] < 1. - f_cut){
      if (posx[] > max_x) {
        max_x = posx[];
        z_at_x_max = posz[]; // Store the z-coordinate corresponding to max_x
      }
      if (posx[] < min_x)
        min_x = posx[];
      if (posz[] > max_z)
        max_z = posz[];
      if (posz[] < min_z)
        min_z = posz[];
      if (posy[] > max_y)
        max_y = posy[];
    }

  /** 
  Calculate sum of energy budget terms. */
  double se  = sigma*(area-pi*sq(R0));
  double gpe = vol*(rho2-rho1)*grav*(zd-z0);
  double energy_sum = ke2 + ke1 + se + DissWork2 + DissWork1 + gpe;

  if ( i == 0 )
    fprintf(ferr,
      "#1:t; 2:dt; 3:xci; 4:uc; 5:zc; 6:wc; 7:y_max; 8:x_max; 9:x_min; 10:z_max; 11:z_min; 12:z_at_xmax; 13:AR_xz; 14:AR_yx 15:vol; 16:area; 17:SE; 18:KE1; 19:KE2; 20:DissRate1; 21:DissRate2; 22:DissWork1; 23:DissWork2; 24:GPE; 25:EnergySum; 26:n_grid; 27:cput; 28:speed; 29:grav; \n");

    fprintf (ferr, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %ld %g %g %g\n",
      t, dt, xd, ud, zd, wd, max_y, max_x, min_x, max_z, min_z, 
      z_at_x_max, (2*max_x)/(max_z-min_z), (2*max_y)/(2*max_x), vol, area,
	    se, ke1, ke2, DissRate1, DissRate2, DissWork1, DissWork2, gpe, energy_sum,
      grid->tn, perf.t, perf.speed, grav);
}


/** Remvoe small bubbles (phase 1) or droplets (phase 2) if they are generated due to VOF advection.
*/
event cleanup_bubbles (i++) {
  remove_droplets(f,5);
}

event cleanup_drops (i++) {
  // Remove small bubbles using remove_droplets with bubbles parameter set to true
  remove_droplets(f,5,true);
}
#endif 


/** 
Adapt mesh based on solid volume fraction, gas volume fraction, and velocity.
*/
event adapt (i++) {
  scalar ibm1[];
  foreach()
    ibm1[] = ibm[];
  adapt_wavelet ({ibm1,f,u}, (double[]){1e-3,femax,uemax,uemax,uemax}, maxlevel=max_level, minlevel=min_level);
}

event dumps (t += t_out; t <= t_end) {
    scalar * dumplist = {u, p, f, ibm};
    char name[80];
    sprintf(name, "../out/dump-%06.4f", t);
    dump(file = name, list = dumplist);
}


#if 1 
/** 
Output images.
*/
event output_images (t += t_out ) {
  char name[80];
  scalar omega[];
  scalar l2[];
 
  /** output figures */
  double image_width=1024;
  double fov1=20;
  //double fov1 = 40;
  double tx1=-0.5;
  double ty1=-0.5;
  double tz1=-0.5;

  foreach() {
    omega_y[] = (u.z[1,0,0] - u.z[-1,0,0] - u.x[0,0,1] + u.x[0,0,-1])/(2*Delta);
  }

  /** view 1: xz plan */
  view (fov = fov1*1, tx=tx1, ty=ty1, tz=0, quat = {0.707,0,0,0.707}, width = image_width, height = image_width);

  clear();
  draw_vof ("f"); mirror ( n= {0,1,0} ) { draw_vof ("f");} 
  squares("u.x", linear = true, n = {0,1,0}, max = vmax, min = -vmax, cbar = true, border = true, pos = {0.47, -0.9}, label = "u.x (m/s)", mid = true, levels = 10);
  cells(n={0,1,0});
  sprintf (name, "../out/vof-ux_view1_%06.4f.png",t);
  save(file=name);
  
  clear();
  draw_vof ("f");mirror ( n= {0,1,0} ) { draw_vof ("f");} 
  squares("u.z", linear = true, n = {0,1,0}, max = vmax, min = -vmax, cbar = true, border = true, pos = {0.47, -0.9}, label = "u.z (m/s)", mid = true, levels = 10); 
  sprintf (name, "../out/vof-uz_view1_%06.4f.png",t);
  save(file=name);

  clear();
  draw_vof ("f");mirror ( n= {0,1,0} ) { draw_vof ("f");} 
  squares("(1-f)*p", linear = true, n = {0,1,0}, max = pmax, min = -pmax, cbar = true, border = true, pos = {0.47, -0.9}, label = "p (Pa)", mid = true, levels = 20); 
  cells(n={0,1,0});
  sprintf (name, "../out/vof-p_view1_%06.4f.png",t);
  save(file=name);

  clear();
  draw_vof ("f");mirror ( n= {0,1,0} ) { draw_vof ("f");}
  squares("(1-f)*p", linear = true, n = {0,1,0}, max = pmax, min = -pmax, cbar = true, border = true, pos = {0.47, -0.9}, label = "p (Pa)", mid = true, levels = 20);
  sprintf (name, "../out/vof-p_view2_%06.4f.png",t);
  save(file=name);


  view (fov = fov1*1, tx=tx1, ty=ty1, tz=0, quat = {0.707, 0, 0, 0.707}, width = image_width, height = image_width);
  clear();
  vorticity (u, omega);
  draw_vof ("f");mirror ( n= {0,1,0} ) { draw_vof ("f");}
  squares ("omega", linear = true);
  lambda2 (u, l2);
  stats s = statsf(l2);
  isosurface ("l2", -0.01);
  sprintf(name,"../out/iso_no_vof-%06.4f.png",t);
  save(file=name);

  clear();
  draw_vof ("f");mirror ( n= {0,1,0} ) { draw_vof ("f");} 
  squares("omega_y", linear = true, n = {0,1,0}, cbar = true, max = 5, min = -5, border = true, pos = {0.47, -0.9}, label = "y-vorticity", mid = true, levels = 20);
  sprintf (name, "../out/vof-omega_y_view1_%06.4f.png",t);
  save(file=name);

  /** view4: yz plane */
  view (fov = fov1*0.8, tx=-0.65, ty=0.06, tz=tz1, quat = {0.707,0,0.707,0}, width = image_width, height = image_width);
  clear();
  squares("ibm+ 2*f", linear = true, n = {1,0,0}, alpha = x0, max = 2, min = 0);
  cells(n={1,0,0},alpha = x0);
  sprintf (name, "../out/vof_view4_%06.4f.png",t);
  save(file=name);

}

event interface (t += 5*t_out) {

   char names[80];
   sprintf(names, "../out/interface%d", pid());
   FILE * fp = fopen (names, "w");
   output_facets (f,fp);
   fclose(fp);
   char command[80];
   sprintf(command, "LC_ALL=C  cat interfa* > ../out/infc%06.4f.dat",t);
   system(command);

}

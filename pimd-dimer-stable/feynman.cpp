/* ----------------------------------------------------------------------   
   This file contains the functions for single-partition job
   -------------------------------------------------------------------
   * For other informations please rper to fix_feynman.h
------------------------------------------------------------------------- */

#include "fix_feynman.h"

#include "verlet.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "atom.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "output.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "timer.h"
#include "memory.h"
#include "random_park.h"
#include "error.h"
#include "universe.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- *
   Initilize the positions and velocities of the beads on the ring-polymer
/* ---------------------------------------------------------------------- */   

void FixFeynman::rp_init()
{
  // random number generator
  random = new RanPark(lmp, 920416+11*comm->me);
  // mass
  
  if(ms_restart==false)
  {
    if(comm->me==0 && screen) fprintf(screen,"FEYNMAN| Generate masses for ring-polymers.\n");

    for(int i=0; i<atom->nlocal; i++)
    {
      mass[i] = atom->mass[atom->type[i]];
    
      for(int j=0; j<input.nbeads; j++)
      {      
        if(input.nm) 
        {
          if(j==0) mass_rp[j][i] = mass[i];
          else mass_rp[j][i] = input.mscale * mass[i];
        }
        else mass_rp[j][i] = mass[i] / input.nbeads;
      }
    }
  }
  
  if(xv_restart) return;
   
  // Init x/v by external file
  
  FILE *fp = fopen("xv.init","r");
    
  if(fp)
  {  
    if(comm->me==0 && screen) fprintf(screen,"FEYNMAN| Reading ring-polymers from file [xv.init].\n");
    
    char line[1001];
    
    while(fgets(line,1000,fp))
    {
      int _iatom, _ibead;
      double _x[3], _v[3];
      
      _iatom = atoi(strtok(line," \t\n"));
      _ibead = atoi(strtok(NULL," \t\n"));
      for(int _i=0; _i<3; _i++) _x[_i] = atof(strtok(NULL," \t\n"));
      for(int _i=0; _i<3; _i++) _v[_i] = atof(strtok(NULL," \t\n"));
      
      int index = 0;
      for(index=0; index<atom->nlocal; index++) if(atom->tag[index]==_iatom) break;
      
      if(index < atom->nlocal) for(int d=0; d<3; d++)
      {
        x_rp[_ibead-1][index][d] = _x[d];
        v_rp[_ibead-1][index][d] = _v[d];
      }
    }
    
    if(input.nm)
    {
      nm_transform(x_rp, M_x2xp, atom->nlocal);
      nm_transform(v_rp, M_x2xp, atom->nlocal);
    }
    
    compute_centroid(x_rp, atom->x, atom->nlocal);
    compute_centroid(v_rp, atom->v, atom->nlocal);
    
    if(input.nm) nm_transform(x_rp, M_xp2x, atom->nlocal);
    
    fclose(fp);
    return;
  }
  if(fp) fclose(fp);
  
  //printf("initial v[110] = %12.8lf v_rp = %12.8lf\n", atom->v[110][0], v_rp[0][110][0]);
  // Init x/v by random assignment
  
  if(comm->me==0 && screen) fprintf(screen,"FEYNMAN| Initialize ring-polymers randomly.\n");
  //srand(static_cast<int>(MPI_Wtime()));
  srand(atom->nlocal);

  for(int i=0; i<atom->nlocal; i++) if(atom->mask[i] & groupbit)
  {
    double r     = 0.01 * rand() / RAND_MAX;
    double phi   = M_PI * rand() / RAND_MAX;
    double theta = M_PI * rand() / RAND_MAX;
    
    for(int j=0; j<input.nbeads; j++)
    {            
      // Position
      
      double rx = r * cos(M_PI * 2.0 * j / input.nbeads); 
      double ry = r * sin(M_PI * 2.0 * j / input.nbeads); 
      
      x_rp[j][i][0] =  rx * cos (phi) + atom->x[i][0];
      x_rp[j][i][1] =  rx * sin(theta)*sin(phi) + ry * cos(theta) + atom->x[i][1];
      x_rp[j][i][2] = -rx * cos(theta)*sin(phi) + ry * sin(theta) + atom->x[i][2];
      
      // Velocity
      
      if(input.nm && j==0) 
      {
        v_rp[j][i][0] = atom->v[i][0];
        v_rp[j][i][1] = atom->v[i][1];
        v_rp[j][i][2] = atom->v[i][2];
      }
      else
      {
        double v_phi   = M_PI * 1.0 * rand() / RAND_MAX;
        double v_theta = M_PI * 2.0 * rand() / RAND_MAX;
    
        double v = sqrt ( atom->v[i][0] * atom->v[i][0]
                          + atom->v[i][1] * atom->v[i][1]
                          + atom->v[i][2] * atom->v[i][2] );
                      
        v *= sqrt(mass[i]/mass_rp[j][i]);
        v_rp[j][i][0] = v * cos(v_theta) * sin(v_phi);
        v_rp[j][i][1] = v * sin(v_theta) * sin(v_phi); 
        v_rp[j][i][2] = v * cos(v_phi);        
        
      }
    }
  }
  else
  {
    for(int j=0; j<input.nbeads; j++)
    {                
      x_rp[j][i][0] = atom->x[i][0];
      x_rp[j][i][1] = atom->x[i][1];
      x_rp[j][i][2] = atom->x[i][2];
      v_rp[j][i][0] = atom->v[i][0];
      v_rp[j][i][1] = atom->v[i][1];
      v_rp[j][i][2] = atom->v[i][2];
    }
  }
  //printf("_initial v[110] = %12.8lf v_rp = %12.8lf\n", atom->v[110][0], v_rp[0][110][0]);
}

/* ---------------------------------------------------------------------- *
   Compute centroid position, force and velocity
/* ---------------------------------------------------------------------- */   

void FixFeynman::compute_centroid(double*** rp, double** c, int n)
{
  if(input.nm) 
  {
    double scale = sqrt(input.nbeads);
    for(int i=0; i<n; i++) if(atom->mask[i] & groupbit)
    {
      c[i][0] = rp[0][i][0] / scale;
      c[i][1] = rp[0][i][1] / scale;
      c[i][2] = rp[0][i][2] / scale;
    }
  }
  else
  {
    for(int i=0; i<n; i++) if(atom->mask[i] & groupbit) for(int d=0; d<3; d++)
    {
      c[i][d] = 0.0;
      for(int ibead=0; ibead<input.nbeads; ibead++) c[i][d] += rp[ibead][i][d];
      c[i][d] /= input.nbeads;
    }
  } 
}

void FixFeynman::compute_centroid_force(double*** rp, double** c, int n)
{
   double scale = (double)input.nbeads;
   for (int i=0; i<n; i++) for (int d=0; d<3; d++)
   {
      c[i][d] = 0.0;
      for (int ibead=0; ibead<input.nbeads; ibead++) c[i][d] += rp[ibead][i][d];
      c[i][d] /= scale;
   }
}

void FixFeynman::reduce_bead(double*** rp, double** c, int n)
{
   for (int i=0; i<n; i++)
   {
      //if (atom->tag[i]==3540) 
      if (atom->mask[i] & groupbit) 
      {
	c[i][2]=rp[0][i][2];
	c[i][1]=rp[0][i][1];
	c[i][0]=rp[0][i][0];
	//break;
      }
   }
}

/* ---------------------------------------------------------------------- */

void FixFeynman::normal_modes_init()
{
  memory->create(M_x2xp, input.nbeads, input.nbeads, "fix_feynman:M_x2xp");
  memory->create(M_xp2x, input.nbeads, input.nbeads, "fix_feynman:M_xp2x");
  memory->create(M_v2vp, input.nbeads, input.nbeads, "fix_feynman:M_v2vp");
  memory->create(M_vp2v, input.nbeads, input.nbeads, "fix_feynman:M_vp2v");
  memory->create(M_xp2x_C, input.nprime, input.nprime, "fix_feynman:M_xp2x_C");
  memory->create(T_x2p, input.nprime, input.nbeads, "fix_feynman:T_x2p");

  lam = (double*) memory->smalloc(sizeof(double)*input.nbeads, "fix_feynman:lam");
  cosomg = (double*) memory->smalloc(sizeof(double)*input.nbeads, "fix_feynman:cosomg");
  sinomg = (double*) memory->smalloc(sizeof(double)*input.nbeads, "fix_feynman:sinomg");
  trans_buf = (double*) memory->smalloc(sizeof(double)*input.nbeads, "fix_feynman:trans_buf");
  
  // Set up  eigenvalues
  const double Boltzmann = force->boltz;
  const double Plank	 = force->hplanck;
  double hbar = Plank / (2.0*M_PI);
  double hbar2 = hbar * input.pscale;
  double beta = 1.0 / (Boltzmann * input.temp);
  beta_n = beta / (double)input.nbeads;
  omega_n = 1.0 / (beta_n*hbar2);  

  lam[0] = 0.0;
  for (int k=1; k<input.nbeads; k++)
  {
     lam[k] = 2.0 * omega_n * sin ( M_PI * k / (double)input.nbeads);
     cosomg[k] = cos( lam[k] * update->dt);
     sinomg[k] = sin( lam[k] * update->dt);
  }
  
  // Set up eigenvectors for non-degenerated modes
  
  for (int j=0; j<input.nbeads; j++)
  {
     int beads_half = input.nbeads/2;
     M_xp2x[j][0] = sqrt(1.0/input.nbeads);
     for (int k=1; k < beads_half; k++)
     {
	M_xp2x[j][k] = sqrt(2.0/input.nbeads) * cos ( 2.0 * M_PI * (j+1) * k / input.nbeads);
     }
     if ( j%2 == 0 )  M_xp2x[j][beads_half] = - sqrt(1.0/input.nbeads);
     else M_xp2x[j][beads_half] = sqrt(1.0/input.nbeads);
     for (int k=beads_half+1; k<input.nbeads; k++)
     {
	M_xp2x[j][k] = sqrt(2.0/input.nbeads) * sin ( 2.0 * M_PI * (j+1) * k / input.nbeads);
     }
  }
  
  
  // Set up Ut
  
  for(int i=0; i<input.nbeads; i++) 
    for(int j=0; j<input.nbeads; j++)
    {
      M_x2xp[i][j] = M_xp2x[j][i]; 
    }
  
  for(int i=0; i<input.nbeads; i++) 
    for(int j=0; j<input.nbeads; j++)
    {
      M_v2vp[i][j] = M_x2xp[i][j];
      M_vp2v[i][j] = M_xp2x[i][j];
    }

  // set up M_PRIME for ring-polymer contraction
  for (int j=0; j<input.nprime; j++)
  {
     M_xp2x_C[j][0] = sqrt(1.0/input.nprime);
     for (int k=1; k<=(input.nprime-1)/2; k++)
     {
	M_xp2x_C[j][k] = sqrt(2.0/input.nprime) * cos ( 2.0 * M_PI * (j+1) * k / input.nprime);
     }
     for (int k=(input.nprime+1)/2; k<input.nprime; k++)
     {
	M_xp2x_C[j][k] = sqrt(2.0/input.nprime) * sin ( 2.0 * M_PI * (j+1) * k / input.nprime);
     }
  }

  double fac_t = sqrt( double(input.nprime) / double(input.nbeads) );
  // set up T_J'J for ring polymer contraction
  for (int j_prime=0; j_prime<input.nprime; j_prime++)
  {
     for (int j=0; j<input.nbeads; j++)
     {
	T_x2p[j_prime][j] = 0;
	for (int k=0; k<=(input.nprime-1)/2; k++)
	{
	   T_x2p[j_prime][j] += fac_t * M_xp2x[j][k] * M_xp2x_C[j_prime][k];
	}
	for (int k=1; k<=(input.nprime-1)/2; k++)
	{
	   T_x2p[j_prime][j] += fac_t * M_xp2x[j][input.nbeads-k] * M_xp2x_C[j_prime][input.nprime-k];
	}
     }
  }

  /**** debug ****/
  
  if(comm->me==0)
  {
    FILE* fp = fopen("NM.info","w");
    fprintf(fp, "EIGENVALUES\n\n");
    for(int i=0; i<input.nbeads; i++) fprintf(fp,"%3d %12.8lf\n", i+1, lam[i]);
    fprintf(fp,"\n");
    
    fprintf(fp, "EIGENVECTORS\n\n");
    for(int i=0; i<input.nbeads; i++) 
    {
      for(int j=0; j<input.nbeads; j++) fprintf(fp, " %12.8lf", M_x2xp[j][i]);
      fprintf(fp,"\n");
    }
   
    fprintf(fp, "EIGENVECTORS->C'\n\n");
    for (int i=0; i<input.nprime; i++)
    {
      for(int j=0; j<input.nprime; j++)  fprintf(fp, " %12.8lf", M_xp2x_C[i][j]);
      fprintf(fp, "\n");
    } 

    fprintf(fp, "TransformationMatrix->T'\n\n");
    for (int i=0; i<input.nprime; i++)
    {
      for(int j=0; j<input.nbeads; j++)  fprintf(fp, " %12.8lf", T_x2p[i][j]);
      fprintf(fp, "\n");
    }
 
    fprintf(fp, "OMEGA_N\n\n");
    fprintf(fp, "%d %12.8lf %12.8lf %12.8f\n", input.nbeads,omega_n, beta_n, omega_nbeads); 

    fclose(fp);
  }
  
  /**************/
}

/* ---------------------------------------------------------------------- *
   Evaluate the forces by PIMD algorithm
/* ---------------------------------------------------------------------- */   

void FixFeynman::effective_force()
{
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int eflag = true, vflag = false;
  
  int world_size;
  MPI_Comm_size(world, &world_size);

  double **x_save = atom->x;
  double **f_save = atom->f;

  int np = input.nbeads;
  if(input.contract>1) np /= input.contract;
  
  double npinv = 1.0 / np;
  double nbinv = 1.0 / input.nbeads;
  double ncinv = 1.0 / input.contract;
  
  double rmix = input.rmix;
  
  // V(p_i)

  evdwl = ecoul = ebond = eangle = edihed = eimp = ekspace = ekinetic = 0.0;
  
  for(int i=0; i<input.nbeads; i++) for(int j=0; j<nall; j++)
  {
    f_rp[i][j][0] = 0.0;
    f_rp[i][j][1] = 0.0;
    f_rp[i][j][2] = 0.0;
  }

  if (lmp_kspace)
  {
    for (int i=0; i<nall; i++)
    {
	x_save[i][0] = x_rp[0][i][0];
	x_save[i][1] = x_rp[0][i][1];
	x_save[i][2] = x_rp[0][i][2];
    }

    memset(&(atom->f[0][0]), 0, nall*sizeof(double)*3);
   
    lmp_kspace->compute(eflag, vflag);
    ekspace += lmp_kspace->energy;
   
    for (int ip=0; ip<input.contract; ip++)
    {
      double fscale = contract_force[ip];
      if (ip==1) ip = input.nbeads/2 - 1; 
      for (int i=0; i<nall; i++) {
	f_rp[ip][i][0] += f_save[i][0] * fscale;
	f_rp[ip][i][1] += f_save[i][1] * fscale;
	f_rp[ip][i][2] += f_save[i][2] * fscale;
      }
    }
  }

  for (int ip=0; ip<input.contract; ip++)
  {
    double fscale = contract_force[ip];
    double escale = fscale;
    if (ip==1) ip= input.nbeads/2 - 1;
    for (int i=0; i<nall; i++)
    {
      x_save[i][0] = x_rp[ip][i][0];
      x_save[i][1] = x_rp[ip][i][1];
      x_save[i][2] = x_rp[ip][i][2];
    }
    
    memset(&(atom->f[0][0]), 0, nall*sizeof(double)*3);

    if (lmp_pair) {
	lmp_pair->compute(eflag, vflag);
    }

    if (atom->molecular) {
	if (lmp_bond) lmp_bond->compute(eflag, vflag);
	if (lmp_angle) lmp_angle->compute(eflag, vflag);
	if (lmp_dihedral) lmp_dihedral->compute(eflag, vflag);
	if (lmp_improper) lmp_improper->compute(eflag, vflag);
    }

    if(lmp_bond)	ebond	+= lmp_bond->energy * escale;
    if(lmp_angle)	eangle	+= lmp_angle->energy * escale;
    if(lmp_dihedral)	edihed	+= lmp_dihedral->energy * escale;
    if(lmp_improper)	eimp	+= lmp_improper->energy * escale;
    if(lmp_pair)	evdwl	+= lmp_pair->eng_vdwl * escale;
    if(lmp_pair)	ecoul	+= lmp_pair->eng_coul * escale;

    //assigning the force values

    for (int i=0; i<nall; i++)
    {
	f_rp[ip][i][0] += f_save[i][0] * fscale;
	f_rp[ip][i][1] += f_save[i][1] * fscale;
	f_rp[ip][i][2] += f_save[i][2] * fscale;
    }

  }

  /* Other forces equals to zero */
  if (input.contract<input.nbeads)
  {
    //for (int ip=input.contract; ip<input.nbeads; ip++)
    for (int ip=0; ip<input.nbeads; ip++)
    {
      if (ip==0)
	continue;
      else if (input.contract==2 && ip== input.nbeads/2 - 1)
	continue;
      else
      {
        for (int i=0; i<nall; i++)
        {
	  f_rp[ip][i][0] = 0.0;
	  f_rp[ip][i][1] = 0.0;
	  f_rp[ip][i][2] = 0.0;
        }
      }
    }
  }

  if(lmp_pair)		lmp_pair->eng_vdwl = evdwl;
  if(lmp_pair)		lmp_pair->eng_coul = ecoul;
  if(lmp_bond)		lmp_bond->energy = ebond;
  if(lmp_angle)		lmp_angle->energy = eangle;
  if(lmp_dihedral)	lmp_dihedral->energy = edihed;
  if(lmp_improper)	lmp_improper->energy = eimp;
  //if(lmp_kspace)	lmp_kspace->energy = ekspace *= 1.0/(double)world_size;
  if(lmp_kspace)	lmp_kspace->energy = ekspace *= (double)input.nbeads/(double)world_size;
  // Conctraction calculation of KSPACE

  double espring = 0.0;
  double x_c = input.center; 
  double k = input.k_spring;

  for(int ip=0; ip<input.contract; ip++) for(int i=0; i<atom->nlocal; i++)
  {
    double fscale = contract_force[ip];
    if(ip==1) ip = input.nbeads/2 - 1;
    if (atom->tag[i] == 3540) 
    {
       f_rp[ip][i][2] -= fscale * k * (x_rp[ip][i][2] - x_c);
       espring += fscale * 0.5 * k * (x_rp[ip][i][2] - x_c) * (x_rp[ip][i][2] - x_c);
    }
  }

  atom->x = x_save;
  atom->f = f_save;
  
  // Centroid 
  compute_centroid_force(f_rp, atom->f, atom->nlocal);

  last_compute = update->ntimestep;
  double pe = evdwl + ecoul + ebond + eangle + edihed + eimp + ekspace + espring;
  pe = pe / (double)input.nbeads;
  MPI_Allreduce(&pe, &pe_total, 1, MPI_DOUBLE, MPI_SUM, world);  
}

void FixFeynman::apply_wall()
{
  double epsilon = 4.0;
  double rc = 2.0;
  double left = -13.5;
  double right = 13.5;
  for (int i=0; i<atom->nlocal; i++)
  {
    //if (mask[i] & groupbit)
    if (atom->tag[i] == 3540)
    {
       for (int ii=0; ii<input.nbeads; ii++)
       {
	  if ((x_rp[ii][i][2] - left < rc))
	  {
	     f_rp[ii][i][2] += -epsilon * 2.0 * (x_rp[ii][i][2] - left - rc);
	  }
	  else if ((right - x_rp[ii][i][2]) < rc)
          {
	     f_rp[ii][i][2] += epsilon * 2.0 * (right - x_rp[ii][i][2] - rc);
	  }
       }	
       break;
    }
  }
}

double FixFeynman::apply_spring()
{
    double x_c = input.center; 
    double k = input.k_spring;
    double espring = 0.0;

    for(int ip=0; ip<input.contract; ip++) for(int i=0; i<atom->nlocal; i++)
    {
      double fscale = contract_force[ip];
      if(ip==1) ip = input.nbeads/2 - 1;
      if (atom->tag[i] == 3540) 
      {
	 f_rp[ip][i][2] -= fscale * k * (x_rp[ip][i][2] - x_c);
	 espring += fscale * 0.5 * k * (x_rp[ip][i][2] - x_c) * (x_rp[ip][i][2] - x_c);
      }
    }

    return espring;
}

/* ---------------------------------------------------------------------- */

void FixFeynman::pbc_ring()
{
  for(int i=0; i<atom->nlocal; i++)
  {
    double dx = 0.0, dy = 0.0, dz = 0.0;
    
    if(x_rp[0][i][0] - atom->x[i][0] > domain->xprd_half) dx = - domain->xprd;
    else if(atom->x[i][0] - x_rp[0][i][0] > domain->xprd_half) dx = domain->xprd;
  
    if(x_rp[0][i][1] - atom->x[i][1] > domain->yprd_half) dy = - domain->yprd;
    else if(atom->x[i][1] - x_rp[0][i][1] > domain->yprd_half) dy = domain->yprd;
  
    if(x_rp[0][i][2] - atom->x[i][2] > domain->zprd_half) dz = - domain->zprd;
    else if(atom->x[i][2] - x_rp[0][i][2] > domain->zprd_half) dz = domain->zprd;
    
    for(int j=0; j<input.nbeads; j++)
    {
      x_rp[j][i][0] += dx;
      x_rp[j][i][1] += dy;
      x_rp[j][i][2] += dz;
    }
  }
}

/* ---------------------------------------------------------------------- *
   Integrator
/* ---------------------------------------------------------------------- */   
void FixFeynman::langevin_update_v()
{
    double c1_0 = exp (-dthalf * input.friction);
    double c2_0 = sqrt (1.0 - c1_0*c1_0);
    double factor = 1.0 / force->mvv2e;
	   //printf("before v = %12.8lf\n", v_rp[1][110][0]);
    for (int ii=0; ii<input.nbeads; ii++)
    for (int i=0; i<atom->nlocal; i++)
    if (atom->mask[i] & groupbit)
    {
	if (ii == 0)
	{
	  v_rp[ii][i][0] = c1_0 * v_rp[ii][i][0] + c2_0 * sqrt(factor/(beta_n * mass_rp[0][i])) * random->gaussian();
	  v_rp[ii][i][1] = c1_0 * v_rp[ii][i][1] + c2_0 * sqrt(factor/(beta_n * mass_rp[0][i])) * random->gaussian();
	  v_rp[ii][i][2] = c1_0 * v_rp[ii][i][2] + c2_0 * sqrt(factor/(beta_n * mass_rp[0][i])) * random->gaussian();
	}
	else
	{
	   double c1 = exp (-dthalf * 2.0 * lam[ii]);
	   double c2 = sqrt (1.0 - c1*c1);
	   v_rp[ii][i][0] = c1 * v_rp[ii][i][0] + c2 * sqrt(factor/(beta_n * mass_rp[ii][i])) * random->gaussian();
	   v_rp[ii][i][1] = c1 * v_rp[ii][i][1] + c2 * sqrt(factor/(beta_n * mass_rp[ii][i])) * random->gaussian();
	   v_rp[ii][i][2] = c1 * v_rp[ii][i][2] + c2 * sqrt(factor/(beta_n * mass_rp[ii][i])) * random->gaussian();
	}	
    }
}

void FixFeynman::update_v_vv()
{
  for(int ii=0; ii<input.nbeads; ii++)
    for(int i=0; i<atom->nlocal; i++)
    if(atom->mask[i] & groupbit)
    {
      double dtfm = dtf / mass_rp[ii][i];
      v_rp[ii][i][0] += dtfm * f_rp[ii][i][0];
      v_rp[ii][i][1] += dtfm * f_rp[ii][i][1];
      v_rp[ii][i][2] += dtfm * f_rp[ii][i][2];
    }
    else
    {
      if(!input.nm || ii==0)
	{
          v_rp[ii][i][0] = atom->v[i][0];
          v_rp[ii][i][1] = atom->v[i][1];
          v_rp[ii][i][2] = atom->v[i][2];
	}
	else
	{
          v_rp[ii][i][0] = 0.0;
          v_rp[ii][i][1] = 0.0;
          v_rp[ii][i][2] = 0.0;
	}
    }
}


void FixFeynman::update_x()
{
  // update quasi-beads -> exact evolution of ring polymer beads!
  for(int ii=0; ii<input.nbeads; ii++) 
  for(int i=0; i<atom->nlocal; i++) 
    if(atom->mask[i] & groupbit)
    {
     if (ii==0)
     {
	x_rp[ii][i][0] += v_rp[ii][i][0] * dtv;
	x_rp[ii][i][1] += v_rp[ii][i][1] * dtv;
	x_rp[ii][i][2] += v_rp[ii][i][2] * dtv;
     }
     else
     {
	//double wk = lam[ii] / sqrt(mass_rp[ii][i]);	
	double v_new = cosomg[ii] * v_rp[ii][i][0] - lam[ii] * sinomg[ii] * x_rp[ii][i][0];
	x_rp[ii][i][0] = (sinomg[ii] / lam[ii]) * v_rp[ii][i][0] + cosomg[ii] * x_rp[ii][i][0];
	v_rp[ii][i][0] = v_new;

	v_new = cosomg[ii] * v_rp[ii][i][1] - lam[ii] * sinomg[ii] * x_rp[ii][i][1];
	x_rp[ii][i][1] = (sinomg[ii] / lam[ii]) * v_rp[ii][i][1] + cosomg[ii] * x_rp[ii][i][1];
	v_rp[ii][i][1] = v_new;

	v_new = cosomg[ii] * v_rp[ii][i][2] - lam[ii] * sinomg[ii] * x_rp[ii][i][2];
	x_rp[ii][i][2] = (sinomg[ii] / lam[ii]) * v_rp[ii][i][2] + cosomg[ii] * x_rp[ii][i][2];
	v_rp[ii][i][2] = v_new;
     }
    }
    else
    {
	if (!input.nm || ii==0)
	{
	   x_rp[ii][i][0] = atom->x[i][0] * sqrt(input.nbeads);
	   x_rp[ii][i][1] = atom->x[i][1] * sqrt(input.nbeads);
	   x_rp[ii][i][2] = atom->x[i][2] * sqrt(input.nbeads);
	}
	else
	{
	   x_rp[ii][i][0] = 0.0;
	   x_rp[ii][i][1] = 0.0;
	   x_rp[ii][i][2] = 0.0;
	}
    }
}

/* ---------------------------------------------------------------------- *
   Normal-modes transform
/* ---------------------------------------------------------------------- */   


void FixFeynman::nm_transform(double*** p, double** T, int n)
{
  for(int i=0; i<n; i++) 
    for(int d=0; d<3; d++)
    {
      // transform p[i][d]
      
      for(int ib=0; ib<input.nbeads; ib++)
      { 
        trans_buf[ib] = p[ib][i][d];
        p[ib][i][d] = 0.0;
      }
      
      for(int ib=0; ib<input.nbeads; ib++)  // mode [ib]
        for(int j=0; j<input.nbeads; j++)   // component of particle [j]
          p[ib][i][d] += (trans_buf[j] * T[ib][j]);
    }
}


/* ---------------------------------------------------------------------- *
   Compute vectors:
   
     (1) ring-potential
     (2) ring-temperature
     (3) ring-scattering
     
/* ---------------------------------------------------------------------- */   

double FixFeynman::compute_ring_temp()
{
  double ke_total = 0.0;

  int ibead_start;
  if(input.method==CMD) ibead_start=1;
  else ibead_start=0;

  for(int ibead=ibead_start; ibead<input.nbeads; ibead++)
  {
    for (int iatom=0; iatom<atom->nlocal; iatom++)
    {
      double *v1 = v_rp[ibead][iatom];
      double kecurrent = mass_rp[ibead][iatom] * (v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]) * force->mvv2e;
      ke_total += kecurrent;
    }
  }

  double scalar[2], scalar_all[2];
  scalar[0] = atom->nlocal;
  scalar[1] = ke_total;
  MPI_Allreduce(scalar,scalar_all,2,MPI_DOUBLE,MPI_SUM,world);
  double ring_temp = scalar_all[1]/(3.0 * force->boltz)/(scalar_all[0]*(input.nbeads-ibead_start));
  
  return ring_temp;
}  

/* ---------------------------------------------------------------------- */

double FixFeynman::compute_ring_dev()
{
  double max_dev2 = 0.0;
  
  for(int i=0; i<atom->nlocal; i++)
    for(int j=0; j<input.nbeads; j++)
    {
      double r2 = 0.0;
      
      for(int d=0; d<3; d++) 
      {
        double dr = x_rp[j][i][d] - atom->x[i][d];
        r2 += dr*dr;
      }
      
      if(r2 > max_dev2) max_dev2 = r2;

    }
  
  double max_dev2_all;
  MPI_Allreduce(&max_dev2, &max_dev2_all, 1, MPI_DOUBLE, MPI_MAX, world);
  return sqrt(max_dev2_all);
}

/* ---------------------------------------------------------------------- */
//FOR DEBUGGING ONLY
void FixFeynman::print_x()
{
   for (int i=0; i<atom->nlocal; i++)
	printf("%d\t%12.8lf\t%12.8lf\t%12.8lf\t%d\t%d\n", atom->type[i], atom->x[i][0], atom->x[i][1], atom->x[i][0], universe->me, atom->nlocal);
}

void FixFeynman::print_coordinates()
{
   for(int i=0; i<atom->nlocal; i++)
     for(int j=0; j<input.nbeads; j++)
     {
	printf("Atom Type %d Bead %d, x= %12.8lf y= %12.8lf z=%12.8lf proc %d\n",atom->type[i],j,x_rp[j][i][0],x_rp[j][i][1],x_rp[j][i][2], universe->me);
     }
}
void FixFeynman::print_velocities()
{
   for(int i=0; i<atom->nlocal; i++)
     for(int j=0; j<input.nbeads; j++)
     {
	printf("Atom %d Bead %d, v_x= %12.8lf v_y= %12.8lf v_z=%12.8lf\n",i,j,v_rp[j][i][0],v_rp[j][i][1],v_rp[j][i][2]);
     }
}
// For HREMD
void FixFeynman::effective_potential()
{
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int eflag = true, vflag = false;
  
  int world_size;
  MPI_Comm_size(world, &world_size);

  double **x_save = atom->x;
  double **f_save = atom->f;

  int np = input.nbeads;
  if(input.contract>1) np /= input.contract;
  
  double npinv = 1.0 / np;
  double nbinv = 1.0 / input.nbeads;
  double ncinv = 1.0 / input.contract;
  
  double rmix = input.rmix;
  
  // V(p_i)

  evdwl = ecoul = ebond = eangle = edihed = eimp = ekspace = ekinetic = 0.0;
  
  if (lmp_kspace)
  {
    for (int i=0; i<nall; i++)
    {
	x_save[i][0] = x_rp[0][i][0];
	x_save[i][1] = x_rp[0][i][1];
	x_save[i][2] = x_rp[0][i][2];
    }

    memset(&(atom->f[0][0]), 0, nall*sizeof(double)*3);
   
    lmp_kspace->compute(eflag, vflag);
    ekspace += lmp_kspace->energy;
  }

  for (int ip=0; ip<input.contract; ip++)
  {
    double fscale = contract_force[ip];
    double escale = fscale;
    if(ip==1) ip= input.nbeads/2 - 1;
    for (int i=0; i<nall; i++)
    {
      x_save[i][0] = x_rp[ip][i][0];
      x_save[i][1] = x_rp[ip][i][1];
      x_save[i][2] = x_rp[ip][i][2];
    }
    
    memset(&(atom->f[0][0]), 0, nall*sizeof(double)*3);

    if (lmp_pair) {
	lmp_pair->compute(eflag, vflag);
    }

    if (atom->molecular) {
	if (lmp_bond) lmp_bond->compute(eflag, vflag);
	if (lmp_angle) lmp_angle->compute(eflag, vflag);
	if (lmp_dihedral) lmp_dihedral->compute(eflag, vflag);
	if (lmp_improper) lmp_improper->compute(eflag, vflag);
    }

    if(lmp_bond)	ebond	+= lmp_bond->energy * escale;
    if(lmp_angle)	eangle	+= lmp_angle->energy * escale;
    if(lmp_dihedral)	edihed	+= lmp_dihedral->energy * escale;
    if(lmp_improper)	eimp	+= lmp_improper->energy * escale;
    if(lmp_pair)	evdwl	+= lmp_pair->eng_vdwl * escale;
    if(lmp_pair)	ecoul	+= lmp_pair->eng_coul * escale;

    //assigning the force values
  }

  if(lmp_pair)		lmp_pair->eng_vdwl = evdwl;
  if(lmp_pair)		lmp_pair->eng_coul = ecoul;
  if(lmp_bond)		lmp_bond->energy = ebond;
  if(lmp_angle)		lmp_angle->energy = eangle;
  if(lmp_dihedral)	lmp_dihedral->energy = edihed;
  if(lmp_improper)	lmp_improper->energy = eimp;
  if(lmp_kspace)	lmp_kspace->energy = ekspace *= (double)input.nbeads/(double)world_size;
  // Conctraction calculation of KSPACE
 
    double x_c = input.center; 
    double k = input.k_spring;
    double espring = 0.0;

    for(int ip=0; ip<input.contract; ip++) for(int i=0; i<atom->nlocal; i++)
    {
      double fscale = contract_force[ip];
      if(ip==1) ip = input.nbeads/2 - 1;
      if (atom->tag[i] == 3540) 
      {
	 espring += fscale * 0.5 * k * (x_rp[ip][i][2] - x_c) * (x_rp[ip][i][2] - x_c);
      }
    }

  atom->x = x_save;
  atom->f = f_save;
  
  last_compute = update->ntimestep;
  double pe = evdwl + ecoul + ebond + eangle + edihed + eimp + ekspace + espring;
  pe = pe / (double)input.nbeads;
  MPI_Allreduce(&pe, &pe_total, 1, MPI_DOUBLE, MPI_SUM, world);  
}


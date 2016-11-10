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
  
  double fscale = 1.0;
  double escale = input.contract;

  // V(p_i)

  evdwl = ecoul = ebond = eangle = edihed = eimp = ekspace = ekinetic = 0.0;
  
  for(int i=0; i<input.nbeads; i++) for(int j=0; j<nall; j++)
  {
    f_rp[i][j][0] = 0.0;
    f_rp[i][j][1] = 0.0;
    f_rp[i][j][2] = 0.0;
  }

  // Conctraction calculation of KSPACE
  
  if (lmp_kspace)
  {
    for (int i=0; i<nall; i++)
    {
       f_save[i][0] = 0.0;
       f_save[i][1] = 0.0;
       f_save[i][2] = 0.0;
       x_save[i][0] = 0.0;
       x_save[i][1] = 0.0;
       x_save[i][2] = 0.0;
    }

    for(int j=0; j<input.nbeads; j++) for(int i=0; i<nall; i++)
    {
       x_save[i][0] += x_rp[j][i][0] * nbinv;
       x_save[i][1] += x_rp[j][i][1] * nbinv;
       x_save[i][2] += x_rp[j][i][2] * nbinv;
    }
   
    //timer->stamp();
    lmp_kspace->compute(eflag, vflag);
    //timer->stamp(TIME_KSPACE);

    ekspace += lmp_kspace->energy * input.nbeads;
    
    for(int i=0; i<input.nbeads; i++) for (int j=0; j<nall; j++)
    {
       f_rp[i][j][0] += f_save[j][0] * input.nbeads;
       f_rp[i][j][1] += f_save[j][1] * input.nbeads;
       f_rp[i][j][2] += f_save[j][2] * input.nbeads;
    }
   
  }
  

  if (lmp_pair)
  {
    double scale = (double)input.nbeads / (double)input.nprime;

    for(int ip=0; ip<input.nprime; ip++)
    {
       for (int i=0; i<nall; i++)
       {
         x_save[i][0] = 0.0;
         x_save[i][1] = 0.0;
         x_save[i][2] = 0.0;
       }
  
       for (int j=0; j<input.nbeads; j++)
       {
	 for(int i=0; i<nall; i++)
	 {
	    x_save[i][0] += T_x2p[ip][j] * x_rp[j][i][0];
	    x_save[i][1] += T_x2p[ip][j] * x_rp[j][i][1];
	    x_save[i][2] += T_x2p[ip][j] * x_rp[j][i][2];
	 }
       }
 
       memset(&(atom->f[0][0]), 0, nall*sizeof(double)*3);
       //timer->stamp();
       lmp_pair->compute(eflag, vflag);
       //timer->stamp(TIME_PAIR);

       evdwl += lmp_pair->eng_vdwl * scale;
       ecoul += lmp_pair->eng_coul * scale;
	
       // assign the forces
       for (int j=0; j<input.nbeads; j++)
       {
       	 for(int i=0; i<nall; i++)
         {
	   f_rp[j][i][0] += f_save[i][0] * T_x2p[ip][j] * scale;
	   f_rp[j][i][1] += f_save[i][1] * T_x2p[ip][j] * scale;
	   f_rp[j][i][2] += f_save[i][2] * T_x2p[ip][j] * scale;
         }
       }
     }
  }
 
  /*for (int ip=0; ip<input.nbeads; ip++)
  {
     for (int i=0; i<nall; i++)
     {
	x_save[i][0] = x_rp[ip][i][0];
	x_save[i][1] = x_rp[ip][i][1];
	x_save[i][2] = x_rp[ip][i][2];
     }
     memset(&(atom->f[0][0]), 0, nall*sizeof(double)*3);
     timer->stamp();
     
     if (atom->molecular)  {
        if (lmp_bond) lmp_bond->compute(eflag, vflag);
	if (lmp_angle) lmp_angle->compute(eflag, vflag);
	if (lmp_dihedral) lmp_dihedral->compute(eflag, vflag);
        if (lmp_improper) lmp_improper->compute(eflag, vflag);
	timer->stamp(TIME_BOND);
     }

     if (lmp_bond) ebond+= lmp_bond->energy;
     if (lmp_angle) eangle+= lmp_angle->energy;
     if (lmp_dihedral) edihed+= lmp_dihedral->energy;
     if (lmp_improper) eimp+= lmp_improper->energy;

     for (int i=0; i<nall; i++)
     {
	f_rp[ip][i][0] += f_save[i][0];
	f_rp[ip][i][1] += f_save[i][1];
	f_rp[ip][i][2] += f_save[i][2];
     }
  }*/ 
  //Applying Contracted Version
  for(int ip=0; ip<np; ip++)
  {
    // assign the coordinates
    
    if(input.contract==1 && rmix<1E-3)
    {
      for (int i=0; i<nall; i++)
      {
         x_save[i][0] = x_rp[ip][i][0];
         x_save[i][1] = x_rp[ip][i][1];
         x_save[i][2] = x_rp[ip][i][2];
      }
    }
    else
    {
      for(int i=0; i<nall; i++) 
      {
        x_save[i][0] = 0.0;
	x_save[i][1] = 0.0;
	x_save[i][2] = 0.0;
      }
      
      int istart = ip * input.contract;
      
      for(int ik=0; ik<input.contract; ik++)
      {
        int j = istart + ik;
	
        for(int i=0; i<nall; i++) 
        {
          x_save[i][0] += x_rp[j][i][0];
	  x_save[i][1] += x_rp[j][i][1];
	  x_save[i][2] += x_rp[j][i][2];
        }
      }
     
      double fc = ncinv * (1-rmix);
      
      for(int i=0; i<nall; i++) 
      {
        x_save[i][0] *= fc;
	x_save[i][1] *= fc;
	x_save[i][2] *= fc;
      }
       
      istart = ip / 2 * input.contract * 2;
      fc = ncinv * rmix * 0.5;
      
      if(np>1) for(int ik=0; ik<input.contract * 2; ik++)
      {
        int j = istart + ik;
	
        for(int i=0; i<nall; i++) 
        {
          x_save[i][0] += x_rp[j][i][0] * fc;
	  x_save[i][1] += x_rp[j][i][1] * fc;
	  x_save[i][2] += x_rp[j][i][2] * fc;
        }
      }
      
    }
   
    // zero all force values
    
    memset(&(atom->f[0][0]), 0, nall*sizeof(double)*3);
    
    // initiate timer
    
    //timer->stamp();
    
    // call force calculators
    // <!---------
 
    if (atom->molecular) {
        if (lmp_bond) lmp_bond->compute(eflag,vflag);
        if (lmp_angle) lmp_angle->compute(eflag,vflag);
        if (lmp_dihedral) lmp_dihedral->compute(eflag,vflag);
        if (lmp_improper) lmp_improper->compute(eflag,vflag);
        //timer->stamp(TIME_BOND);
    }

    // total the energy
    
    if(lmp_bond)     ebond   += lmp_bond->energy;
    if(lmp_angle)    eangle  += lmp_angle->energy;
    if(lmp_dihedral) edihed  += lmp_dihedral->energy;
    if(lmp_improper) eimp    += lmp_improper->energy;
    // -------------!>
    // final force
    
    for(int i=0; i<nall; i++)
    { 
        atom->f[i][0] *= fscale;
        atom->f[i][1] *= fscale;
        atom->f[i][2] *= fscale;
    } 
      
    if(input.contract!=1 || rmix>1E-3)
    {
      int istart = ip * input.contract;
      double fc = 1.0 - rmix;
      
      for(int ik=0; ik<input.contract; ik++)
      {
          int j = istart + ik;
	
          for(int i=0; i<nall; i++) 
          {
            f_rp[j][i][0] += f_save[i][0] * fc;
	    f_rp[j][i][1] += f_save[i][1] * fc;
	    f_rp[j][i][2] += f_save[i][2] * fc;
          }
      }
     
      istart = ip / 2 * input.contract * 2;
      fc = rmix * 0.5;
      
      if(np>1) for(int ik=0; ik<input.contract*2; ik++)
      {
          int j = istart + ik;
	
          for(int i=0; i<nall; i++) 
          {
            f_rp[j][i][0] += f_save[i][0] * fc;
	    f_rp[j][i][1] += f_save[i][1] * fc;
	    f_rp[j][i][2] += f_save[i][2] * fc;
          }
      }
    
    }
    else
    {
	for (int i=0; i<nall; i++)
	{
	  f_rp[ip][i][0] += f_save[i][0];
	  f_rp[ip][i][1] += f_save[i][1];
	  f_rp[ip][i][2] += f_save[i][2];
	}
    }
  }
  
  // Average the energy

  if(lmp_pair)     lmp_pair->eng_vdwl = evdwl;
  if(lmp_pair)     lmp_pair->eng_coul = ecoul;
  if(lmp_bond)     lmp_bond->energy = ebond *= escale;
  if(lmp_angle)    lmp_angle->energy = eangle *= escale;
  if(lmp_dihedral) lmp_dihedral->energy = edihed *= escale;
  if(lmp_improper) lmp_improper->energy = eimp *= escale;
  if(lmp_kspace)   lmp_kspace->energy = ekspace;
  
  // Harmonic Spring for Hp(p_i)
  
  spring_energy = 0.0;

  for(int i=0; i<input.nbeads; i++)
  {
    int j = (i+1) % input.nbeads;
    
    for(int k=0; k<nlocal; k++)
    {
      double delx = x_rp[i][k][0] - x_rp[j][k][0];
      double dely = x_rp[i][k][1] - x_rp[j][k][1];
      double delz = x_rp[i][k][2] - x_rp[j][k][2];
      domain->minimum_image(delx,dely,delz);
      
      double r2 = delx*delx + dely*dely + delz*delz;
      //double fmass = fbond; //I think this may be wrong because the mass was dropped here --Yining 07/06/15
      double fmass = fbond * mass[k] * input.nbeads;
      spring_energy += -0.5 * fmass * r2;
      
    }
  }

  atom->x = x_save;
  atom->f = f_save;

 
  last_compute = update->ntimestep;
  double pe = evdwl + ecoul + ebond + eangle + edihed + eimp + ekspace + spring_energy;
  MPI_Allreduce(&pe, &pe_total, 1, MPI_DOUBLE, MPI_SUM, world);  
}

/* ---------------------------------------------------------------------- */

void FixFeynman::effective_force_respa(int ilevel, int iloop)
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
  
  double fscale = 1.0;
  double escale = input.contract;

  // V(p_i)

  evdwl = ecoul = ebond = eangle = edihed = eimp = ekspace = ekinetic = 0.0;
  
  for(int i=0; i<input.nbeads; i++) for(int j=0; j<nall; j++)
  {
    f_rp[i][j][0] = 0.0;
    f_rp[i][j][1] = 0.0;
    f_rp[i][j][2] = 0.0;
  }

  if (ilevel == level_pair)
  {
    //printf("pair force calculation level %d\n", ilevel);
    double scale = (double)input.nbeads / (double)input.nprime;

    for(int ip=0; ip<input.nprime; ip++)
    {
       for (int i=0; i<nall; i++)
       {
         x_save[i][0] = 0.0;
         x_save[i][1] = 0.0;
         x_save[i][2] = 0.0;
       }
  
       for (int j=0; j<input.nbeads; j++)
       {
	 for(int i=0; i<nall; i++)
	 {
	    x_save[i][0] += T_x2p[ip][j] * x_rp[j][i][0];
	    x_save[i][1] += T_x2p[ip][j] * x_rp[j][i][1];
	    x_save[i][2] += T_x2p[ip][j] * x_rp[j][i][2];
	 }
       }
 
       memset(&(atom->f[0][0]), 0, nall*sizeof(double)*3);
       //timer->stamp();
       lmp_pair->compute(eflag, vflag);
       //timer->stamp(TIME_PAIR);

       evdwl += lmp_pair->eng_vdwl * scale;
       ecoul += lmp_pair->eng_coul * scale;
	
       // assign the forces
       for (int j=0; j<input.nbeads; j++)
       {
       	 for(int i=0; i<nall; i++)
         {
	   f_rp[j][i][0] += f_save[i][0] * T_x2p[ip][j] * scale;
	   f_rp[j][i][1] += f_save[i][1] * T_x2p[ip][j] * scale;
	   f_rp[j][i][2] += f_save[i][2] * T_x2p[ip][j] * scale;
         }
       }
     }   
   }  // end of pair
 
  if (ilevel == level_bond)
  {
    //printf("bond force calculation level %d\n", ilevel);
    for(int ip=0; ip<np; ip++)
    {
      // assign the coordinates
    
      if(input.contract==1 && rmix<1E-3)
      {
        for (int i=0; i<nall; i++)
        {
           x_save[i][0] = x_rp[ip][i][0];
           x_save[i][1] = x_rp[ip][i][1];
           x_save[i][2] = x_rp[ip][i][2];
        }
      }
      else
      {
        for(int i=0; i<nall; i++) 
        {
          x_save[i][0] = 0.0;
	  x_save[i][1] = 0.0;
	  x_save[i][2] = 0.0;
        }
      
        int istart = ip * input.contract;
      
        for(int ik=0; ik<input.contract; ik++)
        {
          int j = istart + ik;	
          for(int i=0; i<nall; i++) 
          {
            x_save[i][0] += x_rp[j][i][0];
	    x_save[i][1] += x_rp[j][i][1];
	    x_save[i][2] += x_rp[j][i][2];
          }
        }
     
        double fc = ncinv * (1-rmix);
      
        for(int i=0; i<nall; i++) 
        {
           x_save[i][0] *= fc;
	   x_save[i][1] *= fc;
	   x_save[i][2] *= fc;
        }
       
        istart = ip / 2 * input.contract * 2;
        fc = ncinv * rmix * 0.5;
      
        if(np>1) for(int ik=0; ik<input.contract * 2; ik++)
        {
          int j = istart + ik;
	
          for(int i=0; i<nall; i++) 
          {
            x_save[i][0] += x_rp[j][i][0] * fc;
	    x_save[i][1] += x_rp[j][i][1] * fc;
	    x_save[i][2] += x_rp[j][i][2] * fc;
          }
        }
      
      }
   
      // zero all force values
    
      memset(&(atom->f[0][0]), 0, nall*sizeof(double)*3);
    
      // initiate timer
    
      //timer->stamp();
    
      // call force calculators
      // <!---------
 
      if (atom->molecular) {
          if (lmp_bond) lmp_bond->compute(eflag,vflag);
          if (lmp_angle) lmp_angle->compute(eflag,vflag);
          if (lmp_dihedral) lmp_dihedral->compute(eflag,vflag);
          if (lmp_improper) lmp_improper->compute(eflag,vflag);
          //timer->stamp(TIME_BOND);
      }

      // total the energy
    
      if(lmp_bond)     ebond   += lmp_bond->energy;
      if(lmp_angle)    eangle  += lmp_angle->energy;
      if(lmp_dihedral) edihed  += lmp_dihedral->energy;
      if(lmp_improper) eimp    += lmp_improper->energy;
      // -------------!>
      // final force
    
      for(int i=0; i<nall; i++)
      { 
         atom->f[i][0] *= fscale;
         atom->f[i][1] *= fscale;
         atom->f[i][2] *= fscale;
      } 
      
      if(input.contract!=1 || rmix>1E-3)
      {
         int istart = ip * input.contract;
         double fc = 1.0 - rmix;
      
         for(int ik=0; ik<input.contract; ik++)
         {
            int j = istart + ik;
	
            for(int i=0; i<nall; i++) 
            {
              f_rp[j][i][0] += f_save[i][0] * fc;
     	      f_rp[j][i][1] += f_save[i][1] * fc;
  	      f_rp[j][i][2] += f_save[i][2] * fc;
            }
         }
     
         istart = ip / 2 * input.contract * 2;
         fc = rmix * 0.5;
        
         if(np>1) for(int ik=0; ik<input.contract*2; ik++)
         {
           int j = istart + ik;
	
           for(int i=0; i<nall; i++) 
           {
             f_rp[j][i][0] += f_save[i][0] * fc;
	     f_rp[j][i][1] += f_save[i][1] * fc;
	     f_rp[j][i][2] += f_save[i][2] * fc;
           }
         }
    
      }
      else
      { 
 	 for (int i=0; i<nall; i++)
	 {
	   f_rp[ip][i][0] += f_save[i][0];
	   f_rp[ip][i][1] += f_save[i][1];
	   f_rp[ip][i][2] += f_save[i][2];
	 }
      }
    }
    
  }

  if (ilevel == level_kspace)
  {
    //printf("kspace force calculation level %d\n", ilevel);
    for (int i=0; i<nall; i++)
    {
       f_save[i][0] = 0.0;
       f_save[i][1] = 0.0;
       f_save[i][2] = 0.0;
       x_save[i][0] = 0.0;
       x_save[i][1] = 0.0;
       x_save[i][2] = 0.0;
    }

    for(int j=0; j<input.nbeads; j++) for(int i=0; i<nall; i++)
    {
       x_save[i][0] += x_rp[j][i][0] * nbinv;
       x_save[i][1] += x_rp[j][i][1] * nbinv;
       x_save[i][2] += x_rp[j][i][2] * nbinv;
    }
   
    //timer->stamp();
    lmp_kspace->compute(eflag, vflag);
    //timer->stamp(TIME_KSPACE);

    ekspace += lmp_kspace->energy * input.nbeads;
    
    for(int i=0; i<input.nbeads; i++) for (int j=0; j<nall; j++)
    {
       f_rp[i][j][0] = f_save[j][0] * input.nbeads;
       f_rp[i][j][1] = f_save[j][1] * input.nbeads;
       f_rp[i][j][2] = f_save[j][2] * input.nbeads;
    }
    
  }
}

/* ---------------------------------------------------------------------- */

void FixFeynman::pbc_ring()
{
  //printf("pbc_ring called!\n");
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
    for (int ii=0; ii<input.nbeads; ii++)
    for (int i=0; i<atom->nlocal; i++)
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
      //printf("update v w/ dtf = %12.8lf f[4][55][0] = %12.8lf v = %12.8lf\n", dtf, f_rp[4][55][0], v_rp[4][55][0]);
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


void FixFeynman::update_v_vv(int ilevel)
{
      //printf("dtf = %12.8lf \n", dtf);
  for(int ii=0; ii<input.nbeads; ii++)
    for(int i=0; i<atom->nlocal; i++)
    if(atom->mask[i] & groupbit)
    {
      double dtfm = dtf / mass_rp[ii][i];
      v_rp[ii][i][0] += dtfm * f_rp_level[ilevel][ii][i][0];
      v_rp[ii][i][1] += dtfm * f_rp_level[ilevel][ii][i][1];
      v_rp[ii][i][2] += dtfm * f_rp_level[ilevel][ii][i][2];
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
      //printf("update x w/ dtv = %12.8lf \n", dtv);
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
	   x_rp[ii][i][0] = atom->x[i][0];
	   x_rp[ii][i][1] = atom->x[i][1];
	   x_rp[ii][i][2] = atom->x[i][2];
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

void FixFeynman::print_coordinates()
{
   for(int i=0; i<atom->nlocal; i++)
     for(int j=0; j<input.nbeads; j++)
     {
	printf("Atom %d Bead %d, x= %12.8lf y= %12.8lf z=%12.8lf\n",i,j,x_rp[j][i][0],x_rp[j][i][1],x_rp[j][i][2]);
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


/* ----------------------------------------------------------------------   
   Package      FixFeynman
   Purpose      Feynman Path Integral Algorithm for Quantum Chemistry
   -------------------------------------------------------------------
   * For other informations please rper to fix_feynman.h
------------------------------------------------------------------------- */

#include "fix_feynman.h"

#define ERR_HEAD "FEYNMAN|"
#define ERR_LAST "Aborted in \"fix/feynman\""

#define STR2UPPER(str) { char *src=str; \
  while(*src) { if(*src>='a' && *src<='z') *src -= 32; src++; } \
}

#define _ECHO if(universe->me==0 && screen) fprintf(screen,
#define _END );

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "output.h"
#include "universe.h"
#include "update.h"

#include "timer.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "neighbor.h"
#include "integrate.h"
#include "respa.h"
#include "fix.h"

#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "dump_feynman.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- *
   The paramters [arg] when constructing a [Fix] object:
   arg[0] = Fix-ID   arg[1] = Group-ID   arg[2] = Fix-Style 
/* ---------------------------------------------------------------------- */   

FixFeynman::FixFeynman(LAMMPS *lmp, int narg, char **arg) : 
Fix(lmp, narg, arg)
{ 
  /* Mission 1: Open input file and read */
  
  _ECHO "\n" _END    
  
  // Step 1: Initialize all input parameters
  
  input.method    = PIMD;
  
  input.nm        = false;
  input.nbeads    = 8;
  input.mscale    = 1.0;
  input.pscale    = 1.0;
  
  input.temp   = 298.15;
  
  input.dump_rp   = 0;
  input.contract  = 1;
  input.rmix      = 0.0;
  
  label = 0;
  
  // Step 2: Check the number of parameters in command line
  
  if(narg!=4) 
  {
    if(comm->me==0)
    {
      _ECHO "%s Wrong number of parameters for [fix/feynman].\n",ERR_HEAD _END
      _ECHO "%s Example: fix fix-ID group-ID feynman input-file\n",ERR_HEAD _END
    }
    error->all(FLERR,ERR_LAST);
  }
  
  // Step 3: Open file at proc 0 and read
  
  if(comm->me==0)
  {
    FILE *fp = fopen(arg[3],"r+");
    
    // If the file cannot be open, STOP the program
    
    if(!fp)
    {
      _ECHO "%s Error: Failed to open input file [%s].\n",ERR_HEAD,arg[3] _END
      error->one(FLERR,ERR_LAST);
    }
    
    char line[1001];
    int iline = 0;
    
    while(fgets(line,1000,fp))
    {
      iline++;
      char *key   = strtok(line, " \t\n");
      char *value = strtok(NULL, " \t\n");
      
      // Skip all blank or comment lines
      
      if(!key || key[0]=='#') continue; 
      
      // If it is an incompleted line, STOP the program
      
      if(!value || value[0]=='#')
      {
        _ECHO "%s Error: Incompleted line at [%s:%d].\n",
          ERR_HEAD,arg[3],iline _END
        error->one(FLERR,ERR_LAST);
      }
      
      // Process the keyword: 
      // 1) Firstly turn all charactors into upper cases
      // 2) Some macros defined here just for beacuty and convinience
      
      STR2UPPER(key);
      STR2UPPER(value);
      
      #define CHECK_START if(0);
      #define CHECK_KEY(cmd) else if(strcmp(key,#cmd)==0)
      #define CHECK_END else
      #define ERROR_LINE _ECHO "%s Error: incorrect value for the keyword \"%s\" at [%s:%d].\n", \
         ERR_HEAD,key, arg[3],iline _END
      
      /* ======================================================================== */
      
      CHECK_START
        
        CHECK_KEY(METHOD) { 
          if     (strcmp(value,"PIMD")==0)   { input.method = PIMD; }          
          else if(strcmp(value,"NMPIMD")==0) { input.method = NMPIMD; }
	  else if(strcmp(value,"CMD")==0)    { input.method = CMD; }
	  else if(strcmp(value,"RPMD")==0)   { input.method = RPMD; }
          else
          {
            ERROR_LINE;
            _ECHO "%s Options: PIMD NMPIMD CMD RPMD.\n", ERR_HEAD _END
            error->one(FLERR,ERR_LAST);
          }
        }

        CHECK_KEY(NBEADS) { 
          input.nbeads = atoi(value);
          if (input.nbeads<3)
          {
            ERROR_LINE;
            _ECHO "%s It should be an integer larger than 2.\n", ERR_HEAD _END
            error->one(FLERR,ERR_LAST);
          }
        }
        
        CHECK_KEY(MSCALE) {           
          input.mscale = atof(value);
	  if (input.mscale<0.0)
          {
            ERROR_LINE;
            _ECHO "%s It should be a positive value.\n", ERR_HEAD _END
            error->one(FLERR,ERR_LAST);
          }
        }
        
	CHECK_KEY(PSCALE) {           
          input.pscale = atof(value);
	  if (input.pscale<0.0)
          {
            ERROR_LINE;
            _ECHO "%s It should be a positive value.\n", ERR_HEAD _END
            error->one(FLERR,ERR_LAST);
          }
        }
	
        CHECK_KEY(TEMP) { 
          input.temp = atof(value);
          if (input.temp<0.0)
          {
            ERROR_LINE;
            _ECHO "%s It should be a positive value.\n", ERR_HEAD _END
            error->one(FLERR,ERR_LAST);
          }
        }
             
        CHECK_KEY(FRICTION) { 
          input.friction = atof(value);
          if (input.friction<0.0)
          {
            ERROR_LINE;
            _ECHO "%s It should be a positive value.\n", ERR_HEAD _END
            error->one(FLERR,ERR_LAST);
          }
        }
             
        CHECK_KEY(DUMP_RP) { 
          input.dump_rp = atoi(value);
          if (input.dump_rp<1)
          {
            ERROR_LINE;
            _ECHO "%s It should be an integer larger than 1.\n", ERR_HEAD _END
            error->one(FLERR,ERR_LAST);
          }
        }
        
        CHECK_KEY(CONTRACT) { 
          input.contract = atoi(value);
          if (input.contract<1 || input.nbeads%input.contract!=0)
          {
            ERROR_LINE;
            _ECHO "%s It should be an integer larger than 1 and (nbeads MOD contract == 0).\n", ERR_HEAD _END
            error->one(FLERR,ERR_LAST);
          }
        }
        
	CHECK_KEY(RMIX) { 
          input.rmix = atof(value);
          if (input.rmix<0)
          {
            ERROR_LINE;
            _ECHO "%s It should be larger than 0.\n", ERR_HEAD _END
            error->one(FLERR,ERR_LAST);
          }
        }
	
	CHECK_KEY(NPRIME) { 
          input.nprime = atoi(value);
          if (input.nprime<1)
          {
            ERROR_LINE;
            _ECHO "%s It should be larger than 0.\n", ERR_HEAD _END
            error->one(FLERR,ERR_LAST);
          }
        }

      CHECK_END {
        _ECHO "%s Error: Unkown keyword \"%s\" at [%s:%d].\n",
          ERR_HEAD,key, arg[3],iline _END
        error->one(FLERR,ERR_LAST);
      }
      
      /* ======================================================================== */
      
      // Echo the setting to the screen
      
      _ECHO "%s %-10s = %s\n", ERR_HEAD, key, value _END
    }
    
    _ECHO "\n" _END 
    fclose(fp);
  }
  
  if(input.method==NMPIMD || input.method==CMD) input.nm = true;
    
  /* Mission 2: Broadcast the paramters */
  
  MPI_Bcast(&input, sizeof(input), MPI_CHAR, 0, world);
  
  /* Mission 3: Set up Fix internal variables and interfaces
  /* Some internal variables are listed following, rpered from [fix.h]
  
     comm_forward           size of forward communication
     peratom_flag           0/1 if per-atom data is stored
     peratom_freq           frequency per-atom data is available at
     size_peratom_cols      0 = vector, N = columns in peratom array
     restart_peratom        1 if Fix saves peratom state, 0 if not
  */
  
  thermo_energy = 1;
  vector_flag = 1;
  extscalar = 1;
  global_freq = 1;
  scalar_flag = 1;
  size_vector = 3;
  extvector   = 0;
  global_freq = 1;
  
  restart_global  = 1;
  restart_peratom = 1;
  peratom_flag    = 1;
  peratom_freq    = 1;
  
  array_atom      = NULL;
  
  //  double mass           +1
  //  
  //  double mass_rp        +1 x nbeads
  //  double x_rp           +3 x nbeads
  //  double v_rp           +3 x nbeads
  //  double f_rp           +3 x nbeads
  //  double f_rp_level     +3 x nbeads x nlevels_respa
    
  size_peratom_cols = input.nbeads * (10 + 3 + 3*nlevels_respa) + 1;
  
  // In leap_frog we have extra v_rp eta, and nh_eta_dot
  
  comm_forward = comm_reverse = 3 * input.nbeads;
  
  /* Mission 4: Set up other variables */

  //printf("FEYNMAN|checking....  %d\n", ((Respa *) update->integrate)->nlevels);
  if (strstr(update->integrate_style, "respa")) {
     nlevels_respa  = ((Respa *) update->integrate)->nlevels;
     step_respa     = ((Respa *) update->integrate)->step;
     level_bond     = ((Respa *) update->integrate)->level_bond;
     level_angle    = ((Respa *) update->integrate)->level_angle;
     level_dihedral = ((Respa *) update->integrate)->level_dihedral;
     level_improper = ((Respa *) update->integrate)->level_improper;
     level_pair     = ((Respa *) update->integrate)->level_pair;
     level_kspace   = ((Respa *) update->integrate)->level_kspace;
  }

  mass = NULL;
  x_rp = f_rp = v_rp = NULL;
  f_rp_level = NULL;
  x_rp = (double***) malloc(sizeof(double**)*input.nbeads);
  f_rp = (double***) malloc(sizeof(double**)*input.nbeads);
  v_rp = (double***) malloc(sizeof(double**)*input.nbeads);
  mass_rp = (double**) malloc(sizeof(double*)*input.nbeads);
  f_rp_level = (double****) malloc(sizeof(double***)*nlevels_respa);
    
  memset(x_rp,0,sizeof(double***)*input.nbeads);
  memset(f_rp,0,sizeof(double***)*input.nbeads);
  memset(v_rp,0,sizeof(double***)*input.nbeads);
  memset(mass_rp,0,sizeof(double**)*input.nbeads);
  memset(f_rp_level,0,sizeof(double****)*nlevels_respa);
    
  /* Mission 5: Initialize the memory for this fix object */
  
  atom->add_callback(0); // Call LAMMPS to allocate memory for per-atom array
  atom->add_callback(1); // Call LAMMPS to re-assign restart-data for per-atom array

  // setup the dump for the ring-polymer
  
  if(input.dump_rp && universe->iworld == 0)
  {
    char fr[10]; sprintf(fr,"%-d", input.dump_rp);
    char nb[10]; sprintf(nb,"%-d", input.nbeads);
    
    char *arg[] = {(char*)"rp",(char*)"all",(char*)"feynman",fr,(char*)"rp.xyz", nb};
    output->add_dump(6,arg);
    ((DumpFeynman*)(output->dump[output->ndump-1]))->fix_feynman = this;
  }
  
  // allocate per-atom array
  
  grow_arrays(atom->nmax);  
  //printf("FEYNMAN|NMAX = %d\n", atom->nmax);
  
  /* Mission 6: other initials */
  
  // setup normal-mode transform
  
  if(input.nm) normal_modes_init();
  
  // restart information
  
  ms_restart = xv_restart = allset = false;
}

FixFeynman::~FixFeynman()
{
  memory->sfree(array_atom);
  memory->sfree(mass);
  
  for (int i=0; i<nlevels_respa; i++)
  {
     for (int j=0; j<input.nbeads; j++)
     {
	memory->destroy(f_rp_level[i][j]);
     }
  }

  for(int i=0; i<input.nbeads; i++)
  {
    memory->destroy(x_rp[i]);
    memory->destroy(f_rp[i]);
    memory->destroy(v_rp[i]);
    memory->destroy(mass_rp[i]);
  }
  
  if(x_rp) free(x_rp);
  if(f_rp) free(f_rp);
  if(v_rp) free(v_rp);
  if(mass_rp) free(mass_rp);
  if(f_rp_level) free(f_rp_level);
}

/* ---------------------------------------------------------------------- *
   The following functions are all overloading from [Fix] class. They 
   are used for setting up the interface to other LAMMPS module, i.e.,
   [atom], [update], [dump] ... For explicit explanation people can 
   rper to their declarations and definitions in [fix.h] and [fix.cpp]
/* ---------------------------------------------------------------------- */

int FixFeynman::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE;
  mask |= FINAL_INTEGRATE_RESPA;
  mask |= PRE_FORCE;
  mask |= PRE_FORCE_RESPA;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixFeynman::init()
{ 
  if(universe->me==0 && screen) fprintf(screen,"%s Initiating Path-Integral ...\n",ERR_HEAD);
  
  last_compute = 0;
  
  // prepare the constants
  
  dtv    = update->dt;
  dtf    = 0.5   * update->dt * force->ftm2v;
  dthalf = 0.5   * update->dt;
  dt4    = 0.25  * update->dt;
  dt8    = 0.125 * update->dt;
  
  inverse_np = 1.0 / input.nbeads;
  
  /* The first solution for the force constant, using SI units 
  
  const double Boltzmann = 1.3806488E-23;    // SI unit: J/K
  const double Plank     = 6.6260755E-34;    // SI unit: m^2 kg / s
  
  double hbar = Plank / ( 2.0 * M_PI );
  double beta = 1.0 / ( Boltzmann * input.nh_temp);
  
  // - P / ( beta^2 * hbar^2)   SI unit: s^-2
  double _fbond = -1.0 / (beta*beta*hbar*hbar) * input.nbeads;
  
  // convert the units: s^-2 -> (kcal/mol) / (g/mol) / (A^2)
  fbond = _fbond * 4.184E+26;
  
  */
  
  /* The current solution, using LAMMPS internal real units */
  
  const double Boltzmann = force->boltz; 
  const double Plank     = force->hplanck; 
  
  double hbar   = Plank / ( 2.0 * M_PI );
  double hbar2  = hbar * input.pscale;
  double beta   = 1.0 / ( Boltzmann * input.temp);  
  beta_n  = beta / (double)input.nbeads;
  double _fbond = 1.0 / (beta*beta*hbar2*hbar2) * input.nbeads;
  
  omega_nbeads = sqrt(input.nbeads) / (hbar * beta) * sqrt(force->mvv2e);
  omega_n = 1.0 / (beta_n*hbar2); // Check This! /Yining
  fbond = - _fbond * force->mvv2e;
  
  //printf("BETA_N == %12.8lf, BETA == %12.8lf OMEG_NB == %12.8lf OMEG_N == %12.8lf\n", beta_n, beta, omega_nbeads, omega_n);
  if(universe->me==0)
    printf("%s -P/(beta^2 * hbar^2) = %20.7lE (kcal/mol/A^2)\n\n", ERR_HEAD, fbond);

  rp_init();
 
  // save force pointers

  save_force_pointers();
}

/* ---------------------------------------------------------------------- */

void FixFeynman::setup(int vflag)
{
  if(universe->me==0 && screen) fprintf(screen,"Setting up Path-Integral ...\n");

  // clear the force, energy and virial
  
  int n = atom->nlocal + atom->nghost;
  
  for(int i=0; i<n; i++)
    atom->f[i][0] = atom->f[i][1] = atom->f[i][2] = 0.0;
  
  if(force->pair) force->pair->eng_vdwl = 0.0;
  if(force->pair) force->pair->eng_coul = 0.0;
  if(force->kspace) force->kspace->energy = 0.0;
  if(force->bond) force->bond->energy = 0.0;
  if(force->angle) force->angle->energy = 0.0;
  if(force->dihedral) force->dihedral->energy = 0.0;
  if(force->improper) force->improper->energy = 0.0;
  
  if (force->pair) memset(force->pair->virial,0,sizeof(double)*6);
  if (force->bond) memset(force->bond->virial,0,sizeof(double)*6);
  if (force->angle) memset(force->angle->virial,0,sizeof(double)*6);
  if (force->dihedral) memset(force->dihedral->virial,0,sizeof(double)*6);
  if (force->improper) memset(force->improper->virial,0,sizeof(double)*6);
  if (force->kspace) memset(force->kspace->virial,0,sizeof(double)*6);

  // call PI calculator
  
  pbc_ring();
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixFeynman::pre_force(int vflag)
{
  // clear force pointers
  
  clear_force_pointers();
}

/* ---------------------------------------------------------------------- */

void FixFeynman::post_force(int vflag)
{
  // recover force pointers
  if (force->newton) load_force_pointers();

  // ring pbc
  
  if(neighbor->ago==0) pbc_ring();

  // update coordinates for ghost atoms
  
  timer->stamp();
  comm->forward_comm_fix(this);
  timer->stamp(TIME_COMM);
  // call path-integral force calculator
  
  timer->stamp();

  effective_force();
  timer->stamp(TIME_BOND);

  
  // update force for ghost atoms
    
  timer->stamp();
  comm->reverse_comm_fix(this);
  timer->stamp(TIME_COMM);
}

/* ----------------------------------------------------------------------
   For RESPA force evaluation
------------------------------------------------------------------------- */

void FixFeynman::pre_force_respa(int vflag, int ilevel, int iloop)
{
   clear_force_pointers();  //unload the force pointers, avoiding automatic force calculation
}

void FixFeynman::post_neighbor()
{
   pbc_ring();
}

void FixFeynman::post_force_respa(int vflag, int ilevel, int iloop)
{
   if (force->newton) load_force_pointers();  // recover force pointers
   
   /*if (ilevel == nlevels_respa-1 && neighbor->ago==0) 
   {
      pbc_ring();       
   }*/
   
   if (ilevel == 0)
   {
      comm->forward_comm_fix(this);
   }
  
   effective_force_respa(ilevel, iloop);
  
   comm->reverse_comm_fix(this);
     
}
/* ----------------------------------------------------------------------
   For the purpose of skipping orginal force->compute
------------------------------------------------------------------------- */

void FixFeynman::save_force_pointers()
{
  lmp_pair     = force->pair;
  lmp_bond     = force->bond;
  lmp_angle    = force->angle;
  lmp_dihedral = force->dihedral;
  lmp_improper = force->improper;
  lmp_kspace   = force->kspace;
}

void FixFeynman::clear_force_pointers()
{
  force->pair = NULL;
  force->kspace = NULL;
  update->integrate->init();
  atom->molecular=0;
}

void FixFeynman::load_force_pointers()
{
  force->pair = lmp_pair;
  force->kspace = lmp_kspace;
  update->integrate->init();
  atom->molecular=1;
}

/* ---------------------------------------------------------------------- */

double FixFeynman::compute_vector(int n)
{ 
  if(n==0) return static_cast<double>(label);
  else if(n==1) return compute_ring_temp();
  else if(n==2) return compute_ring_dev();
}

/* ---------------------------------------------------------------------- */

double FixFeynman::compute_scalar()
{
  if(last_compute < update->ntimestep) effective_force();  
  return pe_total;
}

int FixFeynman::modify_param(int narg, char** args)
{
  for(int i=0; i<narg-1; i++)
  {
    if(strcmp(args[i],"contract")==0) input.contract = atoi(args[++i]);
    else if(strcmp(args[i],"pscale")==0) input.pscale = atof(args[++i]);
    else if(strcmp(args[i],"rmix")==0) input.rmix = atof(args[++i]);
    else if(strcmp(args[i],"label")==0) label = atoi(args[++i]);
    else if(strcmp(args[i],"temp")==0) input.temp = atof(args[++i]);
  }

  const double Boltzmann = force->boltz;
  const double Plank     = force->hplanck;

  double hbar   = Plank / ( 2.0 * M_PI );
  double hbar2  = hbar * input.pscale;
  double beta   = 1.0 / ( Boltzmann * input.temp);
  beta_n = beta / (double)input.nbeads;
  double _fbond = 1.0 / (beta*beta*hbar2*hbar2) * input.nbeads;
  omega_nbeads = sqrt(input.nbeads) / (hbar * beta) * sqrt(force->mvv2e);
  omega_n = 1.0 / (beta_n*hbar2); // Check This! /Yining
  fbond = - _fbond * force->mvv2e;
  
  last_compute = 0;
}

/* ---------------------------------------------------------------------- */

void FixFeynman::initial_integrate(int vflag)
{
    langevin_update_v();
    if(input.nm) nm_transform(v_rp, M_vp2v, atom->nlocal);
    update_v_vv();
    if(input.nm) nm_transform(v_rp, M_v2vp, atom->nlocal);
   
    if(input.nm) nm_transform(x_rp, M_x2xp, atom->nlocal);

    update_x();
    compute_centroid(x_rp, atom->x, atom->nlocal);
    if(input.nm) nm_transform(x_rp, M_xp2x, atom->nlocal);
    if(input.nm) nm_transform(v_rp, M_vp2v, atom->nlocal);
}

/* ---------------------------------------------------------------------- */

void FixFeynman::initial_integrate_respa(int vflag, int ilevel, int iloop)
{
   
    dtv = step_respa[ilevel];
    dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
    dthalf = 0.5 * step_respa[ilevel];

    if (ilevel == nlevels_respa-1)
    {
       langevin_update_v();
       if(input.nm) nm_transform(v_rp, M_vp2v, atom->nlocal);
    }

    //update_v_vv(ilevel);
    update_v_vv();
    
    if (ilevel == 0)
    {
       if(input.nm) nm_transform(v_rp, M_v2vp, atom->nlocal);
       if(input.nm) nm_transform(x_rp, M_x2xp, atom->nlocal);
       update_x();
       compute_centroid(x_rp, atom->x, atom->nlocal);
       if(input.nm) nm_transform(x_rp, M_xp2x, atom->nlocal);
       if(input.nm) nm_transform(v_rp, M_vp2v, atom->nlocal);
    }
}

/* ---------------------------------------------------------------------- */

void FixFeynman::final_integrate()
{

  update_v_vv();

  if(input.nm) nm_transform(v_rp, M_v2vp, atom->nlocal);
  // for output
  langevin_update_v();
  compute_centroid(v_rp, atom->v, atom->nlocal);  
}

/* ---------------------------------------------------------------------- */

void FixFeynman::final_integrate_respa(int ilevel, int iloop)
{
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  dthalf = 0.5 * step_respa[ilevel];

  //update_v_vv(ilevel);
  update_v_vv();

  if (ilevel == nlevels_respa-1)
  {

     if(input.nm) nm_transform(v_rp, M_v2vp, atom->nlocal);
     // for output
     langevin_update_v();
     compute_centroid(v_rp, atom->v, atom->nlocal);  
  }

}

/* ---------------------------------------------------------------------- */

void FixFeynman::write_restart(FILE *fp)
{ 
  if(comm->me==0)
  {
    int size = sizeof(INP_option);
    fwrite(&size, sizeof(int), 1, fp);
    fwrite((void*)(&input), size, 1, fp);
  }
}

/* ---------------------------------------------------------------------- */

void FixFeynman::restart(char *buf)
{
  _ECHO "%s This is a job recovered from a RESTART file.\n", ERR_HEAD _END
  
  ms_restart = xv_restart = true;
  
  INP_option opt_save;
  memcpy(&opt_save, buf, sizeof(INP_option));
  
  if(opt_save.nbeads != input.nbeads)
  {
    _ECHO "%s WARNING: New setting for [nbeads], so ring-polymers will be recreated.\n", ERR_HEAD _END
    xv_restart = ms_restart = false;
  }
  
  if(opt_save.method != input.method || opt_save.mscale != input.mscale)
  {
    _ECHO "%s WARNING: New setting for [method] or [gamma2], so [rp_mass] will be recreated.\n", ERR_HEAD _END
    ms_restart = false;
  }
  
}

/* RESPA */
void FixFeynman::copy_f_flevel(int ilevel)
{
   //printf("force stored level: %d\n", ilevel);
   int n = atom->nlocal;
   for (int j=0; j<input.nbeads; j++)
     for (int i=0; i<n; i++)
     {
	f_rp_level[ilevel][j][i][0] = f_rp[j][i][0];
	f_rp_level[ilevel][j][i][1] = f_rp[j][i][1];	
        f_rp_level[ilevel][j][i][2] = f_rp[j][i][2];
     }
}

void FixFeynman::copy_flevel_f(int ilevel)
{
   //printf("force recovered! level: %d\n", ilevel);
   int n = atom->nlocal;

   for (int j=0; j<input.nbeads; j++)
     for (int i=0; i<n; i++)
     {
	f_rp[j][i][0] = f_rp_level[ilevel][j][i][0];
	f_rp[j][i][1] = f_rp_level[ilevel][j][i][1];
	f_rp[j][i][2] = f_rp_level[ilevel][j][i][2];
     }
}

/* --------------------FOR REMD----------------- */ 
void FixFeynman::scale_velocities(double sfactor)
{
   int nlocal = atom->nlocal;
   for (int i=0; i<input.nbeads; i++)
   {
	for (int j=0; j<nlocal; j++)
	{
		v_rp[i][j][0] = v_rp[i][j][0] * sfactor;
		v_rp[i][j][1] = v_rp[i][j][1] * sfactor;
		v_rp[i][j][2] = v_rp[i][j][2] * sfactor;
	}
   }
}

void FixFeynman::reset_target(double new_temp)
{
  input.temp = new_temp;
 
  const double Boltzmann = force->boltz;
  const double Plank     = force->hplanck;

  double hbar   = Plank / ( 2.0 * M_PI );
  double hbar2  = hbar * input.pscale;
  double beta   = 1.0 / ( Boltzmann * input.temp);
  double _fbond = 1.0 / (beta*beta*hbar2*hbar2) * input.nbeads;
  fbond = - _fbond * force->mvv2e;
  
  last_compute = 0;
}

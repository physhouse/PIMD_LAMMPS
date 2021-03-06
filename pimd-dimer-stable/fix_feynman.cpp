
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

#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "neighbor.h"
#include "integrate.h"
#include "timer.h"
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
 
  // spring
  input.k_spring  = 0.0;
  input.center	  = 0.0;
 
  label = 0;
  
  // Step 2: Check the number of parameters in command line
  
  if(narg!=5) 
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
          if (input.contract<1)
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

        CHECK_KEY(CENTER) {
	  input.center = atoi(value);
	}

        CHECK_KEY(K_SPRING) {
	  input.k_spring = atoi(value);
	  if (input.k_spring < 0)
	  {
	    ERROR_LINE;
	    _ECHO "%s K_spring less than 0 is not stable\n", ERR_HEAD _END
	    error->one(FLERR,ERR_LAST);
	  }
	}

      CHECK_END {
        _ECHO "%s Error: Unknown keyword \"%s\" at [%s:%d].\n",
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
  
  // Contracted Force Evaluation
  contract_force = (double*) memory->smalloc(sizeof(double)*input.nbeads, "fix_feynman:force_contract");
  memset(contract_force, 0.0, sizeof(double)*input.nbeads);

  if (comm->me==0)
  {
    FILE *f_p = fopen(arg[4],"r+");
    
    // If the file cannot be open, STOP the program
    
    if(!f_p)
    {
      _ECHO "%s Error: Failed to open force contraction file [%s].\n",ERR_HEAD,arg[4] _END
      error->one(FLERR,ERR_LAST);
    }
    
    char line[1000];
    fgets(line, 1000, f_p);

    char* key = strtok(line, " \t\n");
    char* value = NULL;

    for (int i=0; i<input.contract; i++)
    {
	char* value = strtok(NULL, " \t\n");
	contract_force[i] = atof(value);
    } 
    fclose(f_p);
  }  
  
  MPI_Bcast(contract_force, input.nbeads, MPI_DOUBLE, 0, world);

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
  //
  //  double nh_eta         +3 x nbeads x nh_chain
  //  double nh_eta_dot     +3 x nbeads x nh_chain
  //  double nh_eta_dotdot  +3 x nbeads x nh_chain
  //  double nh_eta_mass    +3 x nbeads x nh_chain
    
  size_peratom_cols = input.nbeads * (10 + 3) + 1;
  
  // In leap_frog we have extra v_rp eta, and nh_eta_dot
  
  comm_forward = comm_reverse = 3 * input.nbeads;
  
  /* Mission 4: Set up other variables */
  
  mass = NULL;
  x_rp = f_rp = v_rp = NULL;
  
  x_rp = (double***) malloc(sizeof(double**)*input.nbeads);
  f_rp = (double***) malloc(sizeof(double**)*input.nbeads);
  v_rp = (double***) malloc(sizeof(double**)*input.nbeads);
  mass_rp = (double**) malloc(sizeof(double*)*input.nbeads);
    
  memset(x_rp,0,sizeof(double***)*input.nbeads);
  memset(f_rp,0,sizeof(double***)*input.nbeads);
  memset(v_rp,0,sizeof(double***)*input.nbeads);
  memset(mass_rp,0,sizeof(double**)*input.nbeads);
    
  /* Mission 5: Initialize the memory for this fix object */
  
  atom->add_callback(0); // Call LAMMPS to allocate memory for per-atom array
  atom->add_callback(1); // Call LAMMPS to re-assign restart-data for per-atom array

  // setup the dump for the ring-polymer
  
  if(input.dump_rp) //&& universe->iworld==0)
  {
    char fr[10]; sprintf(fr,"%-d", input.dump_rp);
    char nb[10]; sprintf(nb,"%-d", input.nbeads);
    char rp[20]; sprintf(rp,"rp.xyz.%d", universe->iworld);
    char *arg[] = {(char*)"rp",(char*)"all",(char*)"feynman",fr,rp, nb};
    output->add_dump(6,arg);
    ((DumpFeynman*)(output->dump[output->ndump-1]))->fix_feynman = this;
  }
  // allocate per-atom array
  
  grow_arrays(atom->nmax);  
  
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
  memory->sfree(contract_force);

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
  mask |= FINAL_INTEGRATE;
  mask |= PRE_FORCE;
  mask |= POST_FORCE;
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
  timer->stamp(TIME_KSPACE);
  // call path-integral force calculator
  timer->stamp(); 
  effective_force();
  timer->stamp(TIME_BOND);
  // update force for ghost atoms
  timer->stamp();
  comm->reverse_comm_fix(this);
  timer->stamp(TIME_KSPACE);
  compute_centroid_force(f_rp, atom->f, atom->nlocal);
  //printf("atom->f[0] = %12.8lf\n", atom->f[105][0]);
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
  else if(n==1) return static_cast<double>(input.contract); //compute_ring_temp();
  else if(n==2)  return compute_ring_dev();
  /*else if(n==2) {
     for (int i=0; i<atom->nlocal; i++) {
	if (atom->tag[i] == 3540) {
	   return f_rp[0][i][2];
	}
     }
  }*/
}

/* ---------------------------------------------------------------------- */

double FixFeynman::compute_scalar()
{
  if(last_compute < update->ntimestep) effective_potential();  
  compute_centroid_force(f_rp, atom->f, atom->nlocal);  // Very Important!! It's this function to keep the forces unchanged
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
    else if(strcmp(args[i],"force")==0) 
	 {
	    memset(contract_force, 0.0, sizeof(double)*input.nbeads);
	    
	    if (comm->me==0)
	    { 
	      FILE* fp = fopen(args[++i], "r+");
              char line[1000];
    	      fgets(line, 1000, fp);

    	      char* key = strtok(line, " \t\n");
    	      char* value = NULL;

    	      for (int i=0; i<input.contract; i++)
              {
		char* value = strtok(NULL, " \t\n");
		contract_force[i] = atof(value);
    	      }
	      fclose(fp);
            }
            
            MPI_Bcast(contract_force, input.nbeads, MPI_DOUBLE, 0, world);
	 }
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
    //printf("feynman called\n");
    //printf("f[100] = %12.8lf %12.8lf v = %12.8lf\n", atom->f[105][0], f_rp[0][105][0], atom->v[105][0]);
    //printf("x[pi] = %12.8lf f = %12.8lf v = %12.8lf\n", x_rp[0][94][0], f_rp[0][94][0], v_rp[0][94][0]);
    compute_centroid(x_rp, atom->x, atom->nlocal);
    langevin_update_v();
    //printf("f[100] = %12.8lf %12.8lf v = %12.8lf\n", atom->f[105][0], f_rp[0][105][0], atom->v[105][0]);
    //printf("x[pi] = %12.8lf f = %12.8lf v = %12.8lf\n", x_rp[0][94][0], f_rp[0][94][0], v_rp[0][94][0]);
    if(input.nm) nm_transform(v_rp, M_vp2v, atom->nlocal);
    update_v_vv();
    if(input.nm) nm_transform(v_rp, M_v2vp, atom->nlocal);
   
    //printf("f[100] = %12.8lf %12.8lf v = %12.8lf\n", atom->f[105][0], f_rp[0][105][0], atom->v[105][0]);
    //printf("x[pi] = %12.8lf f = %12.8lf v = %12.8lf\n", x_rp[0][94][0], f_rp[0][94][0], v_rp[0][94][0]);
    if(input.nm) nm_transform(x_rp, M_x2xp, atom->nlocal);

    update_x();
    compute_centroid(x_rp, atom->x, atom->nlocal);
    if(input.nm) nm_transform(x_rp, M_xp2x, atom->nlocal);
    if(input.nm) nm_transform(v_rp, M_vp2v, atom->nlocal);
    reduce_bead(x_rp, atom->x, atom->nlocal);
    //printf("f[100] = %12.8lf %12.8lf v = %12.8lf\n", atom->f[105][0], f_rp[0][105][0], atom->v[105][0]);
    //printf("x[pi] = %12.8lf f = %12.8lf v = %12.8lf\n", x_rp[0][94][0], f_rp[0][94][0], v_rp[0][94][0]);
    //printf("finalpi: x = %12.8lf, f = %12.8lf, v = %12.8lf\n", atom->x[94][0], atom->f[94][0], atom->v[94][0]);
    //printf("final: x = %12.8lf, f = %12.8lf, v = %12.8lf\n", atom->x[105][0], atom->f[105][0], atom->v[105][0]); 
}

/* ---------------------------------------------------------------------- */

void FixFeynman::final_integrate()
{
    //printf("f[100] = %12.8lf %12.8lf v = %12.8lf\n", atom->f[105][0], f_rp[0][105][0], atom->v[105][0]);
    //printf("x[pi] = %12.8lf f = %12.8lf v = %12.8lf\n", x_rp[0][94][0], f_rp[0][94][0], v_rp[0][94][0]);
  update_v_vv();

    //printf("f[100] = %12.8lf %12.8lf v = %12.8lf\n", atom->f[105][0], f_rp[0][105][0], atom->v[105][0]);
    //printf("x[pi] = %12.8lf f = %12.8lf v = %12.8lf\n", x_rp[0][94][0], f_rp[0][94][0], v_rp[0][94][0]);
  if(input.nm) nm_transform(v_rp, M_v2vp, atom->nlocal);
  // for output
  langevin_update_v();
    //printf("f[100] = %12.8lf %12.8lf v = %12.8lf\n", atom->f[105][0], f_rp[0][105][0], atom->v[105][0]);
    //printf("x[pi] = %12.8lf f = %12.8lf v = %12.8lf\n", x_rp[0][94][0], f_rp[0][94][0], v_rp[0][94][0]);
  compute_centroid(v_rp, atom->v, atom->nlocal);  
  compute_centroid_force(x_rp, atom->x, atom->nlocal);
  reduce_bead(x_rp, atom->x, atom->nlocal);
    //printf("f[100] = %12.8lf %12.8lf v = %12.8lf\n", atom->f[105][0], f_rp[0][105][0], atom->v[105][0]);
    //printf("x[pi] = %12.8lf f = %12.8lf v = %12.8lf\n", x_rp[0][94][0], f_rp[0][94][0], v_rp[0][94][0]);
    //printf("finalpi: x = %12.8lf, f = %12.8lf, v = %12.8lf\n", atom->x[94][0], atom->f[94][0], atom->v[94][0]);
    //printf("final: x = %12.8lf, f = %12.8lf, v = %12.8lf\n", atom->x[105][0], atom->f[105][0], atom->v[105][0]); 
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



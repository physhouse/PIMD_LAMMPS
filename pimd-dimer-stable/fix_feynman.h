/* ----------------------------------------------------------------------
   
   Package      FixFeynman
   Purpose      Feynman Path Integral Algorithm for Quantum Chemistry
   Copyright    Voth Group @ University of Chicago
   Authors      Chris Knight & Yuxing Peng
   
   Updated      Oct-01-2011
   Version      1.0
   
------------------------------------------------------------------------- */


#ifdef FIX_CLASS

  FixStyle(feynman,FixFeynman)

#else

  #ifndef FIX_FEYNMAN_H
  #define FIX_FEYNMAN_H

  #include "fix.h"

/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/ 

namespace LAMMPS_NS {

  class FixFeynman : public Fix 
  {
    
    /* ----------------------------------------------*
       Internal member functions
         1. Constructor & Destructor
         2. Overloaded function from the [Fix] class 
    /* ----------------------------------------------*/
    
    public:
    
      /* Constructor & Destructor Functions */
    
      FixFeynman(class LAMMPS *, int, char **);
      virtual ~FixFeynman();
      
      /* Overloaded functions from the [Fix] class */
      
      void   init();
      void   setup(int);
      int    setmask();
      
      double memory_usage();

      void   grow_arrays(int);
      void   copy_arrays(int, int, int);

      int    pack_forward_comm(int, int *, double *, int, int *);
      void   unpack_forward_comm(int, int, double *);
      
      int    pack_reverse_comm(int, int, double *);
      void   unpack_reverse_comm(int, int*, double *);
      
      int    pack_exchange(int, double *);
      int    unpack_exchange(int, double *);

      void   restart(char *);
      void   write_restart(FILE *);

      int    pack_restart(int, double *);
      void   unpack_restart(int, int);

      int    size_restart(int);
      int    maxsize_restart();
      
      void   initial_integrate(int);
      void   pre_force(int);
      void   post_force(int);
      void   final_integrate();
      
      double compute_vector(int);    
      double compute_scalar();

    /* ----------------------------------------------*
       Self-defined member functions
    /* ----------------------------------------------*/
    
    public:
    
      void pbc_ring();
      void rp_init();
      void normal_modes_init();      
      void effective_force();
      void effective_potential();
      void apply_wall(); // Apply wall interactions to constrain particles
      double apply_spring(); // Umbrella Sampling

      void update_x();      
      void nm_transform(double***,double**, int);      
      void compute_centroid(double***, double**, int);
      void reduce_bead(double***, double**, int);
      void compute_centroid_force(double***, double**, int);
      void print_coordinates();
      void print_x();
      void print_velocities();
 
      double compute_ring_temp();
      double compute_ring_dev();
      double compute_cmd_temp();
      
      void update_v_vv();  // update v by velocity-verlet
      void langevin_update_v(); 
      void scale_velocities(double sfactor); // For Replica Exchange
      void reset_target(double new_temp); // For Replica Exchange
    /* ----------------------------------------------*
       For the purpose of skipping orginal force->compute
    /* ----------------------------------------------*/
    
    public:
     
      class Pair      *lmp_pair;
      class Bond      *lmp_bond;
      class Angle     *lmp_angle;
      class Dihedral  *lmp_dihedral;
      class Improper  *lmp_improper;
      class KSpace    *lmp_kspace;
     
      double evdwl, ecoul, ebond, eangle, edihed, eimp, ekspace, ekinetic;

      void save_force_pointers();
      void load_force_pointers();
      void clear_force_pointers();

    /* ----------------------------------------------*
       Member variables
    /* ----------------------------------------------*/
    
    public:

      struct INP_option
      {
        /* Settings for the path-integral job */
        
	int    method;
         
        int    nbeads;       /* Number of quasi-particles beads on the ring-polymer */
        int    dump_rp;      /* 0 - off, 1 - single files 2- multiple files */
	double mscale;
	double pscale;
	int    contract;     /* contraction scheme */
	double rmix;    
        double friction;      /* bath temperature */
	double temp;
	int    nprime;	      /* Contraction Scheme for Pair Potentials */	
       
        double k_spring;
	double center;

	bool   nm;
      
      } input;
      
      enum { PIMD, NMPIMD, CMD, RPMD };     

      double *contract_force;  // Mark the force evaluation scalings

      double spring_energy;
      
      /* Constants */
      
      double fbond;          /* force constant -k/m=-P/(beta*hbar)**2 for ring-polymer beads */
      double inverse_np;     /* 1.0 / nbeads */
      double omega_nbeads;
      double omega_n;
      double omgk;
      double beta_n;
      double *cosomg;
      double *sinomg;

      double dtv;         
      double dtf;     
      double dthalf;
      double dt4;     
      double dt8;     
      
      class RanPark *random; 
      /* For nomal-mode transform: see JCP_124_154103(2006) Eq. 27 */
      /* The matrix Kij = 2 * delta(i,j) - delta(i,j-1) - delta (i,j+1) */
      
      double **M_x2xp, **M_xp2x, **M_v2vp, **M_vp2v;  /* Unitary matrix of K */
      double **M_xp2x_C;  /* NM Matrix for Contracted Coordinates */
      double **T_x2p;     /* Transformation between coordinates and contracted coordinates */      

      double *lam;       /* Eigen-values of K */
      double *trans_buf; 
      
      /* For Ring-Polymer */
      
      double *mass;          /* saved masses for classical particles */
      
      /* For Restart */
      
      bool ms_restart, xv_restart, allset;
      
      /* For only single-partition use */
       
      double ***x_rp;        /* positions, forces and velocities for ring-polymer beads */
      double ***f_rp;        /* array[i][j][k]   i - beads  j - atoms  k - xyz          */
      double ***v_rp;        /* x and f are always stored in cartesian coord.           */
                             /* in normal-mode v will be stored in normal-mode coord    */
      double **mass_rp;
      
      /* For Nose-hooer chains */
      
      /* for replica-exchange */
      int label;
      bigint last_compute;
      double pe_total;
      virtual int modify_param(int, char**);
  };

/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/

}

  #endif
#endif

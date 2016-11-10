/* ----------------------------------------------------------------------   
   Package      FixFeynman
   Purpose      Feynman Path Integral Algorithm for Quantum Chemistry
   -------------------------------------------------------------------
   * For other informations please rper to fix_feynman.h
------------------------------------------------------------------------- */

#include "fix_feynman.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
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

#include <string.h>
#include <stdlib.h>
#include <math.h>

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixFeynman::memory_usage()
{
  double bytes = 0;
  bytes = atom->nmax * size_peratom_cols * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixFeynman::grow_arrays(int nmax)
{
  memory->grow(array_atom, nmax, size_peratom_cols, "fix_feynman:array_atom");
  memory->grow(mass,nmax,"fix_feynman:mass");
  
  int count = nmax*3;
  
  for(int i=0; i<input.nbeads; i++)
  {
    memory->grow(x_rp[i], nmax, 3, "fix_feynman:x_rp");
    memory->grow(f_rp[i], nmax, 3, "fix_feynman:f_rp");
    memory->grow(v_rp[i], nmax, 3, "fix_feynman:v_rp");
    memory->grow(mass_rp[i], nmax, "fix_feynman:mass_rp");
  }
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixFeynman::copy_arrays(int i, int j,int deflag)
{
  mass[j] = mass[i];
  
  int i_pos = i*3;
  int j_pos = j*3;
  
  for(int ii=0; ii<input.nbeads; ii++)
  {
    mass_rp[ii][j] = mass_rp[ii][i];
    x_rp[ii][j][0] = x_rp[ii][i][0];
    x_rp[ii][j][1] = x_rp[ii][i][1];
    x_rp[ii][j][2] = x_rp[ii][i][2];
    f_rp[ii][j][0] = f_rp[ii][i][0];
    f_rp[ii][j][1] = f_rp[ii][i][1];
    f_rp[ii][j][2] = f_rp[ii][i][2];
    v_rp[ii][j][0] = v_rp[ii][i][0];
    v_rp[ii][j][1] = v_rp[ii][i][1];
    v_rp[ii][j][2] = v_rp[ii][i][2];
  }
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixFeynman::pack_exchange(int i, double *buf)
{
  int offset=0;
  buf[offset++] = mass[i];
  
  int pos = i * 3;
  
  for(int ii=0; ii<input.nbeads; ii++)
  {
    buf[offset++] = mass_rp[ii][i];
    buf[offset++] = x_rp[ii][i][0];
    buf[offset++] = x_rp[ii][i][1];
    buf[offset++] = x_rp[ii][i][2];
    buf[offset++] = f_rp[ii][i][0];
    buf[offset++] = f_rp[ii][i][1];
    buf[offset++] = f_rp[ii][i][2];
    buf[offset++] = v_rp[ii][i][0];
    buf[offset++] = v_rp[ii][i][1];
    buf[offset++] = v_rp[ii][i][2];
  }
  
  return size_peratom_cols;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixFeynman::unpack_exchange(int nlocal, double *buf)
{
  int offset=0;
  mass[nlocal] = buf[offset++];
  
  int pos = nlocal*3;
  
  for(int ii=0; ii<input.nbeads; ii++)
  {
    mass_rp[ii][nlocal] = buf[offset++];
    x_rp[ii][nlocal][0] = buf[offset++];
    x_rp[ii][nlocal][1] = buf[offset++];
    x_rp[ii][nlocal][2] = buf[offset++];
    f_rp[ii][nlocal][0] = buf[offset++];
    f_rp[ii][nlocal][1] = buf[offset++];
    f_rp[ii][nlocal][2] = buf[offset++];
    v_rp[ii][nlocal][0] = buf[offset++];
    v_rp[ii][nlocal][1] = buf[offset++];
    v_rp[ii][nlocal][2] = buf[offset++];
  }

  return size_peratom_cols;
} 

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for File
------------------------------------------------------------------------- */

int FixFeynman::pack_restart(int i, double *buf)
{
  int offset=0;
  buf[offset++] = size_peratom_cols+1;
  buf[offset++] = mass[i];

  int pos = i * 3;

  for(int ii=0; ii<input.nbeads; ii++)
  {
    buf[offset++] = mass_rp[ii][i];
    buf[offset++] = x_rp[ii][i][0];
    buf[offset++] = x_rp[ii][i][1];
    buf[offset++] = x_rp[ii][i][2];
    buf[offset++] = f_rp[ii][i][0];
    buf[offset++] = f_rp[ii][i][1];
    buf[offset++] = f_rp[ii][i][2];
    buf[offset++] = v_rp[ii][i][0];
    buf[offset++] = v_rp[ii][i][1];
    buf[offset++] = v_rp[ii][i][2];
  }
  
  return size_peratom_cols+1;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixFeynman::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values
  
  int m = 0;
  for (int i=0; i<nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  mass[nlocal] = extra[nlocal][m++];
  
  for(int ii=0; ii<input.nbeads; ii++)
  {
    if(ms_restart) mass_rp[ii][nlocal] = extra[nlocal][m++];
    else m++;
    
    if(xv_restart)
    {
      x_rp[ii][nlocal][0] = extra[nlocal][m++];
      x_rp[ii][nlocal][1] = extra[nlocal][m++];
      x_rp[ii][nlocal][2] = extra[nlocal][m++];
      f_rp[ii][nlocal][0] = extra[nlocal][m++];
      f_rp[ii][nlocal][1] = extra[nlocal][m++];
      f_rp[ii][nlocal][2] = extra[nlocal][m++];
      v_rp[ii][nlocal][0] = extra[nlocal][m++];
      v_rp[ii][nlocal][1] = extra[nlocal][m++];
      v_rp[ii][nlocal][2] = extra[nlocal][m++]; 
    }
    else m+=9;
    
    int pos = nlocal * 3;
  }
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixFeynman::maxsize_restart()
{
  return size_peratom_cols+1;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixFeynman::size_restart(int nlocal)
{
  return size_peratom_cols+1;
}

/* ---------------------------------------------------------------------- */

int FixFeynman::pack_forward_comm(int n, int *list, double *buf,
			  int pbc_flag, int *pbc)
{
  int m = 0;
  double dx = 0.0, dy = 0.0, dz = 0.0;

  if(pbc_flag)
  {
    dx = domain->xprd * pbc[0];
    dy = domain->yprd * pbc[1];
    dz = domain->zprd * pbc[2];
  }

  for(int i= 0; i<n; i++) 
  {
    int j = list[i];

    for(int ii=0; ii<input.nbeads; ii++)
    {
      buf[m++] = x_rp[ii][j][0] + dx;
      buf[m++] = x_rp[ii][j][1] + dy;
      buf[m++] = x_rp[ii][j][2] + dz;
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void FixFeynman::unpack_forward_comm(int n, int first, double *buf)
{
  int m = 0;
  int last = first + n;
  
  for(int i=first; i<last; i++) 
  {
    for(int ii=0; ii<input.nbeads; ii++)
    {
      x_rp[ii][i][0] = buf[m++];
      x_rp[ii][i][1] = buf[m++];
      x_rp[ii][i][2] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixFeynman::pack_reverse_comm(int n, int first, double *buf)
{
  int m = 0;
  int last = first + n;
  
  for(int i=first; i<last; i++) 
  {
    for(int ii=0; ii<input.nbeads; ii++)
    {
      buf[m++] = f_rp[ii][i][0];
      buf[m++] = f_rp[ii][i][1];
      buf[m++] = f_rp[ii][i][2];
    }
  } 
  
  return m;
}

/* ---------------------------------------------------------------------- */

void FixFeynman::unpack_reverse_comm(int n, int *list, double *buf)
{
  int m = 0;
  
  for(int i= 0; i<n; i++) 
  {
    int j = list[i];

    for(int ii=0; ii<input.nbeads; ii++)
    {
      f_rp[ii][j][0] += buf[m++];
      f_rp[ii][j][1] += buf[m++];
      f_rp[ii][j][2] += buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

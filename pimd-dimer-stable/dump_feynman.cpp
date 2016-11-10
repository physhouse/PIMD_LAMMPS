/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "string.h"
#include "dump_feynman.h"
#include "fix_feynman.h"
#include "atom.h"
#include "group.h"
#include "error.h"
#include "memory.h"
#include "update.h"
#include "domain.h"

#include "stdlib.h"
#include "string.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

DumpFeynman::DumpFeynman(LAMMPS *lmp, int narg, char **arg) : Dump(lmp, narg, arg)
{
  if (narg != 6) error->all(FLERR,"Illegal dump feynman command");
  if (binary || multiproc) error->all(FLERR,"Invalid dump xyz_rpsp filename");
 
  _nbeads = atoi(arg[5]);
   size_one = 2 + 9 * _nbeads;  

  char *str;
  str = (char *) "%8d %12g %12g %12g %12g %12g %12g %12g %12g %12g";
  
  sort_flag = 1;
  sortcol = 0;

  int n = strlen(str) + 1;
  format_default = new char[n];
  strcpy(format_default,str);
  
  fix_feynman = NULL;
}

/* ---------------------------------------------------------------------- */

void DumpFeynman::init_style()
{
  delete [] format;
  char *str;
  if (format_user) str = format_user;
  else str = format_default;

  int n = strlen(str) + 2;
  format = new char[n];
  strcpy(format,str);
  strcat(format,"\n");

  // open single file, one time only

  if (multifile == 0) openfile();
}

/* ---------------------------------------------------------------------- */

void DumpFeynman::write_header(bigint n)
{
  if(me==0)
  {
    ntotal = n;
    fprintf(fp, BIGINT_FORMAT "\n", ntotal * _nbeads);
    fprintf(fp,"timestep " BIGINT_FORMAT "   lx %-lf   ly %-lf   lz %-lf\n",update->ntimestep, domain->xprd, domain->yprd, domain->zprd);
  }
}

/* ---------------------------------------------------------------------- */

int DumpFeynman::count()
{
  if (igroup == 0) return atom->nlocal;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int m = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) m++;
  return m;
}

/* ---------------------------------------------------------------------- */

void DumpFeynman::pack(int *ids)
{
  int m,n;

  int *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  double ***x = fix_feynman->x_rp;
  double ***v = fix_feynman->v_rp;
  double ***f = fix_feynman->f_rp;
  int nlocal = atom->nlocal;
  
  m = n = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      buf[m++] = tag[i];
      buf[m++] = type[i];
      
      for(int j=0; j<_nbeads; j++)
      {
        buf[m++] = x[j][i][0];
        buf[m++] = x[j][i][1];
        buf[m++] = x[j][i][2];
        buf[m++] = v[j][i][0];
        buf[m++] = v[j][i][1];
        buf[m++] = v[j][i][2];
        buf[m++] = f[j][i][0];
        buf[m++] = f[j][i][1];
        buf[m++] = f[j][i][2];
      }
      
      if (ids) ids[n++] = tag[i];
    }
}

/* ---------------------------------------------------------------------- */

void DumpFeynman::write_data(int n, double *mybuf)
{
  int m=0;

  for (int i=0; i<n; i++) 
  {
    for (int j =0; j<_nbeads; j++)
    {
      fprintf(fp,format, static_cast<int> (mybuf[m+1]),
        mybuf[m+j*9+2],mybuf[m+j*9+3],mybuf[m+j*9+4],
        mybuf[m+j*9+5],mybuf[m+j*9+6],mybuf[m+j*9+7],
        mybuf[m+j*9+8],mybuf[m+j*9+9],mybuf[m+j*9+10]);
    }
      
    m += size_one;
  }
}

/* ----------------------------------------------------------------------
Copyright (2010) Aram Davtyan and Garegin Papoian

Papoian's Group, University of Maryland at Collage Park
http://papoian.chem.umd.edu/

Last Update: 12/01/2010
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(qbias,FixQBias)

#else

#ifndef LMP_FIX_QBIAS_H
#define LMP_FIX_QBIAS_H

#include "fix.h"

//#include "smart_matrix_lib.h"

namespace LAMMPS_NS {

class FixQBias : public Fix {
 public:
  FixQBias(class LAMMPS *, int, char **);
  ~FixQBias();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double compute_scalar();
  double compute_vector(int);

private:
  double epsilon;
  double k_qbias;
  double k_qdiffbias;

  int force_flag;
  int nlevels_respa;
  bool allocated;

  int ntimestep;
  int n, nn;
  int *alpha_carbons;
  bool qbias_flag, qdiff_flag, qbias_exp_flag, qobias_flag, qobias_exp_flag, qinterfacebias_flag;
  int start1, end1, start2, end2;
  double q0, sigma, sigma_exp;
  double cutoff, min_sep;
  char rnative_file[100];
  char rnative2_file[100];
  double *sigma_sq;
  int l;
  double **rN, **rN2, **r, **q, **q2;
  double Q1, Q2;
  int *res_no, *res_info, *chain_no;
  double **x, **f;
  double **xca;
  int *image;
  double prd[3], half_prd[3];
  int *periodicity;

  enum ResInfo{NONE=0, LOCAL, GHOST, OFF};
  
  double energy[5], energy_all[5];
  enum EnergyTerms{ET_TOTAL=0, ET_QBIAS, nEnergyTerms};

private:
  void compute();
  void compute_qbias();
  double compute_qdiffconstant(double **rStructure);
  void compute_qdiffbias();
  void compute_qinterfacebias();

  double Sigma(int sep);

  void allocate();
  int Tag(int index);
  inline void Construct_Computational_Arrays();
  inline double PeriodicityCorrection(double d, int i);
  inline bool isFirst(int index);
  inline bool isLast(int index);
  inline void print_log(char *line);

  int Step;
  int sStep, eStep;
  FILE *fout;
  FILE *efile;
  FILE *ffile;
  void out_xyz_and_force(int coord=0);
};

}

#endif
#endif


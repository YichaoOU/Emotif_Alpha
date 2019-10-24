 #include <iostream>
#include <map>
#include <string>
#include <vector>
#include <set>
#include "glpk.h"
#include <list>
#include <cassert>
#include <cstdlib>
#include <fstream>
using namespace std;
string strip(string &x) {
  string t;
  for (int i=0;i<x.size();i++) {
    if (!isspace(x[i])) t+=x[i];
  }
  return t;
}
//This module is free software. You can redistribute it and/or modify it under 
//the terms of the MIT License, see the file COPYING included with this 
//distribution.
// We are going to solve the set cover problem using 
// the relaxed version of the integer linear programming
// formulation.
// In the language of GLPK, m = number of genes
// n = number of sets/number of motifs
// xi (1<= i <=m)  is going to be the sum of the selected motifs that 
// appear in each gene
// xi = ai1 x_m+1 + ai2 x_m+2 .... + ain x_m+n,
// where aij = 0 if the jth motif does not contain gene i
// and   aij = 1 if the jth motif contains gene i.
//
// Now, the matrix a is stored 
// AS A SPARSE MATRIX.
// The array ia gives the i index of each element.
// The array ja gives the j index of each element.
// The array ar give the value for each entry.
// In this case, the total number of elements is the sum of the sizes of the
// sets.
//
// Notice that we want gene i to be covered.  So, we want
// 1<= xi <= n
// Similarly, we want 0<= x_m+j <=1.
// 
// Overall, this is a straightforward relaxed ILP.
// Here, we are going to try to solve the MIXED Integer Linear Programming
// problem where all of the variables are boolean variables.
// Let's see what happens.
//
vector<double> solve(vector<set<int> > &adj_list, int gene_numb) {
  vector<double> x;
  int n;
  int m;
  n = adj_list.size();
  m = gene_numb;

  int ASIZE=0;
  for (int i=0;i<adj_list.size();i++) {
    ASIZE+=adj_list[i].size();
  }
  
  glp_prob *lp;

  int *ia;
  int *ja;
  double *ar;
  double Z;



  ia = new int [ASIZE+2];
  ja = new int [ASIZE+2];
  ar = new double [ASIZE + 2];

  int index_count = 0;
  for (int i=0;i<adj_list.size();i++) {
    for (set<int>::iterator p=adj_list[i].begin();p!=adj_list[i].end();++p) {
      index_count++;
      ia[index_count] = (*p)+1;
      ja[index_count] = i+1;
      ar[index_count] = 1.0;
    }
  }
  assert(index_count <ASIZE+2);

  lp = lpx_create_prob();
  glp_set_prob_name(lp,"SET COVER");
  glp_set_obj_dir(lp,GLP_MIN);

  //glp_set_int_parm(lp,LPX_K_MSGLEV,1);

  glp_add_rows(lp,m);
  for (int i=1;i<=m;i++) {
    glp_set_row_bnds(lp,i,GLP_DB,1.0,n);
  }
  glp_add_cols(lp,n);
  for (int i=1;i<=n;i++) {
    glp_set_col_bnds(lp,i,GLP_DB,0.0,1.0);
    glp_set_obj_coef(lp,i,1.0);
    glp_set_col_kind(lp, i, GLP_BV);
  }
  cout << n << " " << m << endl;
  
  glp_load_matrix(lp,index_count,ia,ja,ar);

 glp_iocp parm;

 glp_init_iocp(&parm);
 parm.presolve = GLP_ON;

 int err = glp_intopt(lp, &parm);
 cout << "Error Flag " << err << endl;
  double z = glp_mip_obj_val(lp);
  cout << "Optimum MIP Solution = " << z << endl;
  x.resize(n);
  for (int i=1;i<=n;i++) {
    double x1 = glp_mip_col_val(lp, i);    
    cout << x1 << endl;
    x[i-1] = x1;
  }

  glp_delete_prob(lp);
  delete [] ar;
  delete [] ja;
  delete [] ia;
  return x;
}

set<int> randomized(vector<double> &x, vector<set<int> > &adj_list, int m) {

  set<int> covered;
  set<int> cover;

  int tries = 0;
  while ((covered.size() < m) && (tries < 10)) {
    for (int i=0;i<x.size();i++) {
      double r = rand()/(RAND_MAX*1.0);
      //cout << r << " " << x[i] << endl;
      if ( r  < x[i]) {
	for (set<int>::iterator p = adj_list[i].begin();
	     p!=adj_list[i].end();++p) {
	 
	  covered.insert(*p);
	}
	cover.insert(i);
      }
    }
    tries++;
  }
  cout << "Tries = " << tries << endl;
  if (tries == 10) { cout << "Couldn't find one" << endl;}
  return cover;
}
	     
  

int main(int argc,char *argv[]) {
   time_t t1;

  srand(time(&t1)); // Seed the random number generator
                   // with the current number of seconds 
                   // since Jan 1, 1970


  //check for user arguments 
  if (argc < 2)
  {
	  cout << "Please input a motif seq file" << endl;
  } 
  
  string motifFileName = argv[1];
  cout << "motif file name:" << motifFileName <<endl;


  string line;
  int state = 1;
  int current_motif = 0;
  vector<string> motifs;
  map<string,int> motif_numb;
  vector<set<int> > adj_list;

  vector<string> genes;
  map<string,int> gene_numb;



	ifstream motifFile;
  motifFile.open(motifFileName.c_str());
  if (motifFile.is_open())
  {
    while ( getline (motifFile,line) )
	{
	  if (line[0]=='>') 
	  {
	  //cout << "Found Motif" << line << endl;
	  motif_numb[line] = current_motif;
	  current_motif++;
	  motifs.push_back(line);
	  adj_list.resize(current_motif);
	  state = 0;
	  } else 
	  {
	  int c;
	  if (gene_numb.count(line) > 0) 
	  {
	    c = gene_numb[line];
	  } else 
	  {
	    genes.push_back(line);
	    gene_numb[line] = genes.size()-1;
	    c = genes.size()-1;
	  }
	    adj_list[current_motif-1].insert(c);
	  }
	}	
    motifFile.close();
  } else cout << "Unable to open file"; 

  //while (!cin.eof()) {
    //getline(cin,line);

    //if (!cin.fail()) {
	//if (line[0]=='>') {
	  ////cout << "Found Motif" << line << endl;
	  //motif_numb[line] = current_motif;
	  //current_motif++;
	  //motifs.push_back(line);
	  //adj_list.resize(current_motif);
	  //state = 0;
	//} else {
	  //int c;
	  //if (gene_numb.count(line) > 0) {
	    //c = gene_numb[line];
	  //} else {
	    //genes.push_back(line);
	    //gene_numb[line] = genes.size()-1;
	    //c = genes.size()-1;
	  //}
	  //adj_list[current_motif-1].insert(c);
	//}
    //}
  //}
  cout << "Number of Genes = " << genes.size() << endl;
  int m = genes.size();
  vector<int> percentages;
  percentages.resize(motifs.size());
  //cout << "Sets " << endl;
  for (int i=0;i<adj_list.size();i++) {
    //cout << "{ " ;
    for (set<int>::iterator p = adj_list[i].begin();
	 p!=adj_list[i].end(); ++p) {
      //cout  << *p << ",";
    }
    //cout << "}" << endl;
  }
  vector<double> x;
  x = solve(adj_list,genes.size());

  set<int> t;
  set<int> r_cover;

  int count = 0;
  for (int i=0;i<x.size();i++) {
    //cout << x[i] << endl;
    if (x[i] >0.0) {
      count++;
      r_cover.insert(i);
      for (set<int>::iterator p = adj_list[i].begin();p!=adj_list[i].end();++p) {
	t.insert(*p);
      }
      percentages[i]=t.size();
      
    }
  }
  cout << "Direct = " << count << endl;
  cout << t.size() << endl;
  cout << "Found one of size " << r_cover.size() << endl;

  cout << "The results are:" << endl;
  cout << "The no. of features obtained are: " << r_cover.size() << endl;
  cout << "The no. of sequences covered are: " << m << endl;
  cout << "The maximum sequence coverage is: " << 100 << endl;
  int i=0;
  for (set<int>::iterator p=r_cover.begin();
       p!=r_cover.end(); ++p) {
       
    //cout << strip(motifs[*p]) << " " << percentages[*p] << " " << percentages[*p]/(1.0*m) << endl;
    //cout << motifs[*p] << " " << percentages[*p] << " " << percentages[*p]/(1.0*m) << endl;
    cout << motifs[*p] <<endl;
    i++;
  } 
  
  /* set<int> r_cover;
  set<int> b_cover;
  r_cover = randomized(x,adj_list,genes.size());
  int c1 = 0;
  do {
    b_cover = randomized(x,adj_list,genes.size());
    if (b_cover.size() < r_cover.size()) {
      r_cover = b_cover;
    }
    c1++;
  } while (c1 < 500);

  cout << "Found one of size " << r_cover.size() << endl;
  cerr << r_cover.size() << endl;

  for (set<int>::iterator p=r_cover.begin();
       p!=r_cover.end(); ++p) {
    cout << motifs[*p] << endl;
    } */

}

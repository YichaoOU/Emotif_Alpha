/*******************************************************************************************************\
* requried cover - heuristic selection of a set of covering strings					*
*													*
*   usage:	rc [-d <coverage_depth>] <input_file>							*
*													*
*		-d n  sets the "depth of coverage" (minimum number of times each target string is	*
*		      covered).  default is 1.								*
*		input file is in the following file format, each line containing one of the following:	*
*		    a string beginning with ">"								*
*			these lines represent the namea of new covering strings				*
*		    a string not beginning with ">"							*
*			subsequent lines until a new covering string are considered names of adjacent 	*
*			target strings									*
*													*
*  algorithmic notes:											*
*  + covering string names must be unique in appearance.						*
*  + the number of covering strings "required" (min depth on some target string) is reported		*
*  + after required strings, covering set is built by greedy selection on number of target strings	*
*    added to covered set										*
*													*
*  required cover reports various statisics about the data set, and the greedy minimized set of		*
*  strings.  the "percent strings fully covered" is the cumulative total coverage as the covering	*
*  set is built.  "percent adjacencies required" is applicable only to depth of coverage of 2 or	*
*  greater, and represents the percent completion of the cover to the required depth			*
*													*
*  the motivating application for development of this algorithm was bioinformatics network research.	*
*  the covering strings may be viewed as "motifs" or "words" representing TFBSs, while the target	*
*  strings would be genes (or promoter regions) in which the motif is found. 				*
*		
* 
* # This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.
* 											*
* Jeffrey S. Jones, 											*
* Bioinformatics Lab								*
* Ohio University											*
* Summer, 2013												*
\*******************************************************************************************************/
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string.h>
#include <cstdlib>
using namespace std;

#define SUCCESS	0
#define	FAIL	1

#define TEST

int main(int argc, char *argv[])
{
string	 infile_name;
ifstream infile;

string	 instring;
int	 num_cov_strings = 0;
int	 max_tar_strings = 0;

string	 *cov_names;
string	 *tar_names;
int	 cov_loaded;
int	 tar_loaded;
int	 tar_location;
int	 required_depth;

int	 max_cov_name_len;
bool	 **adjacency;
int	 current_cov_string;
int	 current_tar_string;
int	 *tar_cov_depth;
int	 *covering_strings;
int	 *covered_tar_strings;
int	 min_cov_depth;
bool	 *selected_cov_strings;
int	 total_adjacencies;
int	 *cov_string_weight;
int	 greedy_cover_string;
int	 num_selected_cover_strings;
int	 current_covered;
int	 max_covered_target;


/*******************************\
* assign positional parameters  *
\*******************************/
    required_depth = 1;					/* default */

    if ( argc == 4 )
    {
	if ( strncmp(argv[1], "-d", 2) != 0 )
	{
	    cerr << "\nusage:  rc [-d <coverage_depth>] <input_file>"  << endl;
            exit(FAIL);
	}
	required_depth = atoi(argv[2]);			/* should really check for valid int ... */
    }
    else if  ( argc != 2 )
    {
	cerr << "\nusage:  rc [-d <coverage_depth>] <input_file>"  << endl;
        exit(FAIL);
    }

    infile_name = argv[argc-1];
    cout << "solving for depth of coverage " << required_depth << endl;

/*******************************************************\
* count covering strings,				*
*   and maximum number possible for target strings	*
\*******************************************************/
    infile.open(infile_name.c_str());
    if ( infile.fail() )
    {
	cerr << "failed to open input file " << infile_name << endl;
	exit(FAIL);
    }

    while ( getline(infile, instring) )
    {
	if ( instring.c_str()[0] == '>' )
	{
	    num_cov_strings++;
	}
	else
	{
	    max_tar_strings++;
	}
    }
    if ( !infile.eof() )
    {
	cerr << "reading strings ended before EOF, exiting" << endl << endl;
	return(FAIL);
    }
    infile.close();

/*******************************************************\
* load the string names tables				*
\*******************************************************/
    infile.open(infile_name.c_str());
    cov_names = new string[num_cov_strings];
    tar_names = new string[max_tar_strings];

    cout << "loading string names tables\n";
    cov_loaded = 0;
    tar_loaded = 0;
    while ( getline(infile, instring) )
    {
	if ( instring.c_str()[0] == '>' )
	{
	    for ( int i = 0; i < cov_loaded; i++ )
	    {
		if ( instring == cov_names[i] )
		{
		    cerr << "repeated covering string, exting." << endl << endl;
		    return(FAIL);
		}
	    }
	    cov_names[cov_loaded] = instring;
	    cov_loaded++;
	}
	else
	{
	    tar_location = 0;
	    for ( int i = 0; i < tar_loaded; i++ )
	    {
		if ( instring > tar_names[i] ) tar_location++;
	    }
	    if ( tar_location == tar_loaded )		/* adding to end of array */
	    {
		tar_names[tar_location] = instring;
	        tar_loaded++;
	    }
	    else if ( instring != tar_names[tar_location] )
	    {
		for ( int i = tar_loaded; i > tar_location; i-- ) tar_names[i] = tar_names[i-1];
		tar_names[tar_location] = instring;
	        tar_loaded++;
	    }
	}
    }
    if ( !infile.eof() )
    {
	cerr << "string loading ended before EOF, exiting" << endl << endl;
	return(FAIL);
    }
    infile.close();
    cout << "number of covering strings: " << num_cov_strings  << endl;
    cout << "number of target strings:   " << tar_loaded       << endl;
    cout << "number of associations:     " << max_tar_strings  << endl;

/***********************************************************************\
* need max length of cover string name for formatting output later	*
\***********************************************************************/
    max_cov_name_len = cov_names[0].length();
    for (int i = 1; i < num_cov_strings; i++)
    {
	if ( cov_names[i].length() > max_cov_name_len ) max_cov_name_len = cov_names[i].length();
    }

/*******************************************************\
* create and populate an adjacency matrix		*
\*******************************************************/
    adjacency = new bool*[num_cov_strings];
    for (int i = 0; i < num_cov_strings; i++ )
    {
	adjacency[i] = new bool[tar_loaded];
	for (int j = 0; j < tar_loaded; j++)
	{
	    adjacency[i][j] = false;
	}
    }
    infile.open(infile_name.c_str());
    current_cov_string = -1;
    while ( getline(infile, instring) )
    {
	if ( instring.c_str()[0] == '>' )
	{
	    current_cov_string++;
	}
	else
	{
	    current_tar_string = 0;
	    while (  (current_tar_string < tar_loaded)
		  && (instring != tar_names[current_tar_string]) )
	    {
		current_tar_string++;
	    }
	    if ( current_tar_string >= tar_loaded)	/* should never happen */
	    {
		cerr << "loaded target string not found, exiting." << endl << endl;
		return(FAIL);
	    }
	    adjacency[current_cov_string][current_tar_string] = true;
	}
    }
    infile.close();
    cout << "adjacency matrix created" << endl;

/*******************************************************\
* identify required covering strings			*
*   by computing depth of coverage of target strings.	*
*   any with only required depth or less , will require	*
*   the associated covering string in the covering set	*
\*******************************************************/
    cout.setf(ios::fixed);
    cout.setf(ios::showpoint);
    cout.precision(2);
    tar_cov_depth = new int[tar_loaded];
    for (int i = 0; i < tar_loaded; i++)
    {
	tar_cov_depth[i] = 0;
    }
    for (int i = 0; i < num_cov_strings; i++)
    {
	for (int j = 0; j < tar_loaded; j++)
	{
	    if ( adjacency[i][j] )
	    {
		tar_cov_depth[j]++;
	    }
	}
    }

    min_cov_depth = tar_cov_depth[0];
    for (int i = 1; i < tar_loaded; i++)
    {
	if ( tar_cov_depth[i] <  min_cov_depth ) min_cov_depth = tar_cov_depth[i];
    }
    if (min_cov_depth <= 0)		/* should never happen */
    {
	cerr << "impossible zero or negative depth of coverage, exiting." << endl << endl;
	return(FAIL);
    }
    cout << "best possible complete coverage depth is: " << min_cov_depth << endl;

    cout << "                                      ";
    for (int slen = 0; slen < max_cov_name_len; slen++) cout << " ";
    cout << "               percent" << endl;
    cout << "                                      ";
    for (int slen = 0; slen < max_cov_name_len; slen++) cout << " ";
    cout << "   percent     strings" << endl;
    cout << "                                      covering";
    for (int slen = 8; slen < max_cov_name_len; slen++) cout << " ";
    cout << " adjacencies    fully" << endl;
    cout << "                                       string ";
    for (int slen = 8; slen < max_cov_name_len; slen++) cout << " ";
    cout << "  required     covered" << endl;
    cout << "                                      ========";
    for (int slen = 8; slen < max_cov_name_len; slen++) cout << " ";
    cout << " ===========   =======" << endl;

    covered_tar_strings = new int[tar_loaded];
    for (int i = 0; i < tar_loaded; i++)
    {
	covered_tar_strings[i] = 0;
    }

    selected_cov_strings = new bool[num_cov_strings];
    for (int i = 0; i < num_cov_strings; i++) selected_cov_strings[i] = false;

    num_selected_cover_strings = 0;
    for (int i = 0; i < tar_loaded; i++)
    {
	/***************************************************************\
	* if a target string is covered required_depth or fewer times	*
	*   then select all strings which cover it			*
	\***************************************************************/
	if ( tar_cov_depth[i] <= required_depth )
	{
	    for (int j = 0; j < num_cov_strings; j++)
	    {
		if (  adjacency[j][i]
		   && !selected_cov_strings[j] )
		{
		    cout << setw(4) << (num_selected_cover_strings+1) << ":   ";
		    cout << "adding required cover string: " << cov_names[j];
		    for (int slen = cov_names[j].length(); slen < max_cov_name_len; slen++) cout << " ";
		    total_adjacencies = 0;
		    for (int k = 0; k < tar_loaded; k++)
		    {
			if ( adjacency[j][k] )
			{
			    covered_tar_strings[k]++;
			}
			total_adjacencies +=  min(covered_tar_strings[k], required_depth);
		    }
		    selected_cov_strings[j] = true;
		    num_selected_cover_strings++;
		    current_covered = 0;
		    for (int k = 0; k < tar_loaded; k++)
		    {
			if ( covered_tar_strings[k] >= required_depth ) current_covered++;
		    }
		    cout << ":    " << setw(6);
		    cout << ( (float)total_adjacencies / (float)(required_depth*tar_loaded) *100.0);
		    cout << ":    " << setw(6);
		    cout << ( (float)current_covered / (float)(tar_loaded) *100.0) << endl;
		}
	    }
	}
    }

/*******************************************************\
* score each potential covering vertex			*
*   add to selected cov strings until cannot improve	*
*   coverage						*
\*******************************************************/
    cov_string_weight = new int[num_cov_strings];

    for (int i = 0; i < num_cov_strings; i++)
    {
	cov_string_weight[i] = 0;
	if ( !selected_cov_strings[i] )
	{
	    for (int j = 0; j < tar_loaded; j++)
	    {
		if (  adjacency[i][j]
		   && (covered_tar_strings[j] < required_depth) )
	        {
		    cov_string_weight[i]++;
		}
	    }
	}
    }
    greedy_cover_string = 0;
    for (int i = 1; i < num_cov_strings; i++)
    {
	if ( cov_string_weight[i] > cov_string_weight[greedy_cover_string] )
	{
	    greedy_cover_string = i;
	}
    }

    while ( cov_string_weight[greedy_cover_string] > 0 )
    {
    /*******************************************\
    * add the greedy cover			*
    \*******************************************/
	cout << setw(4) << (num_selected_cover_strings+1) << ":   ";
	cout << "adding greedy   cover string: " << cov_names[greedy_cover_string];
	for (int slen = cov_names[greedy_cover_string].length(); slen < max_cov_name_len; slen++) cout << " ";
	total_adjacencies = 0;
	for (int i = 0; i < tar_loaded; i++)
	{
	    if ( adjacency[greedy_cover_string][i] ) covered_tar_strings[i]++;
	    total_adjacencies += min(covered_tar_strings[i], required_depth);
	}
	selected_cov_strings[greedy_cover_string] = true;
	num_selected_cover_strings++;
	current_covered = 0;
	for (int k = 0; k < tar_loaded; k++)
	{
	    if ( covered_tar_strings[k] >= required_depth ) current_covered++;
	}
	cout << ":    " << setw(6);
	cout << ( (float)total_adjacencies / (float)(required_depth*tar_loaded) *100.0);
	cout << ":    " << setw(6);
	cout << ( (float)current_covered / (float)(tar_loaded) *100.0) << endl;

    /*******************************************\
    * re-score and select next cover string	*
    \*******************************************/
	for (int i = 0; i < num_cov_strings; i++)
	{
	    cov_string_weight[i] = 0;
	    if ( !selected_cov_strings[i] )
	    {
		for (int j = 0; j < tar_loaded; j++)
		{
		    if (  adjacency[i][j]
		       && (covered_tar_strings[j] < required_depth) )
	            {
			cov_string_weight[i]++;
		    }
	        }
	    }
	}
	greedy_cover_string = 0;
	for (int i = 1; i < num_cov_strings; i++)
	{
	    if ( cov_string_weight[i] > cov_string_weight[greedy_cover_string] )
	    {
		greedy_cover_string = i;
	    }
	}
    }
    cout << "\n" << setw(6) << ( (float)current_covered / (float)(tar_loaded) *100.0);
    cout << " coverage achieved with " << num_selected_cover_strings << " covering strings" << endl;
    cout << "average coverage density         : " << ((float)tar_loaded/(float)num_selected_cover_strings) << endl;
    
    max_covered_target = covered_tar_strings[0];
    for (int i = 1; i < tar_loaded; i++)
    {
	if ( covered_tar_strings[i] > max_covered_target ) max_covered_target = covered_tar_strings[i];
    }
    cout << "maximum coverage of single target: " << max_covered_target << endl;

    delete [] cov_names;
    delete [] tar_names;
    for (int i = 0; i < num_cov_strings; i++ ) delete [] adjacency[i];
    delete [] adjacency;
    delete [] tar_cov_depth;
    delete [] covered_tar_strings;
    delete [] selected_cov_strings;
    delete [] cov_string_weight;
    return(SUCCESS);
}

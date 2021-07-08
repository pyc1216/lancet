#ifndef VARIANT_HH
#define VARIANT_HH 1

/****************************************************************************
** Variant.hh
**
** Class for storing basic variant information storage
**
**  Authors: Giuseppe Narzisi
**    Date: November 5, 2015
**
*****************************************************************************/

/************************** COPYRIGHT ***************************************
**
** New York Genome Center
**
** SOFTWARE COPYRIGHT NOTICE AGREEMENT
** This software and its documentation are copyright (2016) by the New York
** Genome Center. All rights are reserved. This software is supplied without
** any warranty or guaranteed support whatsoever. The New York Genome Center
** cannot be responsible for its use, misuse, or functionality.
**
** Version: 1.0.0
** Author: Giuseppe Narzisi
**
*************************** /COPYRIGHT **************************************/

#include <string>
#include <iostream>
#include <limits>
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <array>

#include "util.hh"
#include "FET.hh"

using namespace std;

struct Filters
{
	// filter thresholds]
	double minPhredFisherSTR;
	double minPhredFisher;
	double maxVafNormal;
	double minVafTumor;
	int minCovNormal;
	int maxCovNormal;
	int minCovTumor;
	int maxCovTumor;
	int minAltCntTumor;
	int maxAltCntNormal;
	int minStrandBias;
};

class Variant_t
{
public:

	/*
	Type				size	range
	short int			2		-32,768 to 32,767
	unsigned short int	2		0 to 65,535
	unsigned int		4		0 to 4,294,967,295
	int					4		-2,147,483,648 to 2,147,483,647
	long int			4		-2,147,483,648 to 2,147,483,647
	unsigned long int	4		0 to 4,294,967,295
	*/
		
	bool LR_MODE;
	bool isSomatic;
	unsigned short kmer;
	// similar variants(have same getSignature and isSomatic) from different window(default 600bp, step=100bp)
	unsigned short similar_variants_count;
	string chr;
	int pos;
	char type;
	unsigned short len;
	
	string ref;
	string alt;
	string str;
	
	char status; // T=somatic, S=shared
	unsigned short ref_cov_normal_fwd;
	unsigned short ref_cov_normal_rev;
	unsigned short ref_cov_tumor_fwd;
	unsigned short ref_cov_tumor_rev;
	unsigned short alt_cov_normal_fwd;
	unsigned short alt_cov_normal_rev;
	unsigned short alt_cov_tumor_fwd;
	unsigned short alt_cov_tumor_rev;
	
	std::array<unsigned short, 3> HPRN;
	std::array<unsigned short, 3> HPRT; 
	std::array<unsigned short, 3> HPAN;
	std::array<unsigned short, 3> HPAT;
	
	string bxset_ref_N;
	string bxset_ref_T;
	string bxset_alt_N;
	string bxset_alt_T;
	
	char prev_bp_ref; // base-pair preceding the mutation in reference
	char prev_bp_alt; // base-pair preceding the mutation in alternative
	
	Variant_t(bool mode, string chr_, int pos_, string ref_, string alt_, bool isS,
		unsigned short svc,
		const pair <int,int> & RCN, const pair <int,int> & RCT,
		const pair <int,int> & ACN, const pair <int,int> & ACT,
		const array<unsigned short,3> & HPRN_, const array<unsigned short,3> & HPRT_, 
		const array<unsigned short,3> & HPAN_, const array<unsigned short,3> & HPAT_,
		char prev_bp_ref_, char prev_bp_alt_, int k, string str_, char code,
		const string & bxRN, const string & bxRT, const string & bxAN, const string & bxAT)
	{ 				
		LR_MODE = mode;
		kmer = k;

		str = str_;		
		//str = new char [str_.length()+1];
		//strcpy(str, str_.c_str());
		
		chr = chr_;
		//chr = new char [chr_.length()+1];
		//strcpy(chr, chr_.c_str());
		
		pos = pos_;
		isSomatic = isS;
		similar_variants_count = svc;
		
		/*
		if(ref_.at(0) == '-') { type = 'I'; ref_ = ""; len = alt_.length(); }  // deletion
		if(alt_.at(0) == '-') { type = 'D'; alt_ = ""; len = ref_.length(); }  // insertion
		if(ref_.size()==1 && alt_.size()==1 && ref_.at(0)!='-' && alt_.at(0)!='-') { type = 'S'; pos++; } // snp 
		*/
		
		if(code == '^') { type = 'I'; ref_ = ""; len = alt_.length(); }  // deletion
		if(code == 'v') { type = 'D'; alt_ = ""; len = ref_.length(); }  // insertion
		if(code == 'x') { type = 'S'; pos++; } // snp 
		if(code == 'c') { // complex 
			type = 'C'; 
			ref_.erase(remove(ref_.begin(), ref_.end(), '-'), ref_.end()); 			
			alt_.erase(remove(alt_.begin(), alt_.end(), '-'), alt_.end());
			
			unsigned short rl = ref_.length();
			unsigned short al = alt_.length();
			if(rl == al)  { len = al; }
			else if(rl > al) {len = rl-al;}
			else { len = al-rl; }
			//cout << "R:" << ref_.length() << " A:" << alt_.length() << " LEN:" << len << endl;
		}
		
		if(type != 'S') {
			ref = prev_bp_alt_ + ref_;
			alt = prev_bp_alt_ + alt_;
		}
		else { alt = alt_; ref = ref_; len = 1; }
		
		ref_cov_normal_fwd = RCN.first;
		ref_cov_normal_rev = RCN.second;
		ref_cov_tumor_fwd  = RCT.first;
		ref_cov_tumor_rev  = RCT.second;
		
		alt_cov_normal_fwd = ACN.first;
		alt_cov_normal_rev = ACN.second;
		alt_cov_tumor_fwd  = ACT.first;
		alt_cov_tumor_rev  = ACT.second;
		
		prev_bp_ref = prev_bp_ref_;
		prev_bp_alt = prev_bp_alt_;
		
		HPRN = HPRN_;
		HPRT = HPRT_; 
		HPAN = HPAN_;
		HPAT = HPAT_;
		
		bxset_ref_N = bxRN;
		bxset_ref_T = bxRT;
		bxset_alt_N = bxAN;
		bxset_alt_T = bxAT;
		
		//compute genotype
		//reGenotype();
	}
	
	string printVCF(Filters * fs);
	string printVcfWithoutFilters();
	string genotype(int R, int A);
	char bestState(int Rn, int An, int Rt, int At);
	string getSignature() const;
	string getPosition();
	double compute_FET_score();
	double compute_SB_score();
	double compute_HP_score(int hpr1, int hpr2, int hpa1, int hpa2);
	void printSizeOfElements();
};

#endif

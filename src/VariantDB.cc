#include <vector>

#include "VariantDB.hh"

/****************************************************************************
** VariantDB.cc
**
** Class for storing multiple variants (DB-style)
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

// add variant to DB and update counts per position
void VariantDB_t::addVar(const Variant_t & v) {
	
	string key = itos(v.isSomatic) + sha256(v.getSignature());
	// string key = sha256(v.getSignature());
    map<string,Variant_t>::iterator it_v = DB.find(key);	
	// a ratio when somatic status dose not match.
	// float cover_ratio = 1.5;
	
	//bool flag = false;
	
	if (it_v != DB.end()) {
		// keep variant at location with highest total coverage (tumor + normal)
		
		// bool old_isSomatic = it_v->second.isSomatic;
		int old_ref_cov_normal = it_v->second.ref_cov_normal_fwd + it_v->second.ref_cov_normal_rev;
		int old_ref_cov_tumor = it_v->second.ref_cov_tumor_fwd + it_v->second.ref_cov_tumor_rev;
		int old_alt_cov_normal = it_v->second.alt_cov_normal_fwd + it_v->second.alt_cov_normal_rev;
		int old_alt_cov_tumor = it_v->second.alt_cov_tumor_fwd + it_v->second.alt_cov_tumor_rev;
		
		// bool new_isSomatic = v.isSomatic;
		int new_ref_cov_normal = v.ref_cov_normal_fwd + v.ref_cov_normal_rev;
		int new_ref_cov_tumor = v.ref_cov_tumor_fwd + v.ref_cov_tumor_rev;
		int new_alt_cov_normal = v.alt_cov_normal_fwd + v.alt_cov_normal_rev;
		int new_alt_cov_tumor = v.alt_cov_tumor_fwd + v.alt_cov_tumor_rev;
		
		int old_tot_cov = old_ref_cov_normal + old_ref_cov_tumor + old_alt_cov_normal + old_alt_cov_tumor;
		int new_tot_cov = new_ref_cov_normal + new_ref_cov_tumor + new_alt_cov_normal + new_alt_cov_tumor;
		
		// update count
		// unsigned short svc = it_v->second.similar_variants_count;
		it_v->second.similar_variants_count ++;
		if(old_tot_cov < new_tot_cov) {
			/*
			// if somatic status not match and less than cover_ratio, give up update variant.
			cerr << "old_tot_cov * cover_ratio: " << old_tot_cov * cover_ratio << endl;
			cerr << "new_tot_cov: " << new_tot_cov << endl;
			if (old_isSomatic != new_isSomatic && old_tot_cov * cover_ratio > new_tot_cov) {
				cerr << "*** Variant replacement failed ***" << endl;
				cerr << "Key: " << v.getSignature() << endl;
				cerr << "Old var: " << it_v->second.printVcfWithoutFilters() << endl;
				Variant_t tmp = v;
				cerr << "New var: " << tmp.printVcfWithoutFilters() << endl;
				cerr << "***************************" << endl;
			}
			*/
			cerr << "*** Variant replacement ***" << endl;
			cerr << "Signature: " << v.getSignature() << ", isSomatic: " << itos(v.isSomatic) << endl;
			cerr << "Old var: " << it_v->second.printVcfWithoutFilters() << endl;
			Variant_t tmp = v;
			cerr << "New var: " << tmp.printVcfWithoutFilters() << endl;
			cerr << "Current similar variants counts: " << itos(it_v->second.similar_variants_count) <<endl;
			cerr << "***************************" << endl;
			// it_v->second.isSomatic = v.isSomatic;
			it_v->second.kmer = v.kmer;
			
			it_v->second.ref_cov_normal_fwd = v.ref_cov_normal_fwd;
			it_v->second.ref_cov_normal_rev = v.ref_cov_normal_rev;
			it_v->second.ref_cov_tumor_fwd  = v.ref_cov_tumor_fwd;
			it_v->second.ref_cov_tumor_rev  = v.ref_cov_tumor_rev;
			it_v->second.alt_cov_normal_fwd = v.alt_cov_normal_fwd; 
			it_v->second.alt_cov_normal_rev = v.alt_cov_normal_rev;
			it_v->second.alt_cov_tumor_fwd  = v.alt_cov_tumor_fwd;
			it_v->second.alt_cov_tumor_rev  = v.alt_cov_tumor_rev;
			
			it_v->second.HPRN = v.HPRN;
			it_v->second.HPRT = v.HPRT;
			it_v->second.HPAN = v.HPAN;
			it_v->second.HPAT = v.HPAT;
			
			if(LR_MODE) {
				it_v->second.bxset_ref_N = v.bxset_ref_N;
				it_v->second.bxset_ref_T = v.bxset_ref_T;
				it_v->second.bxset_alt_N = v.bxset_alt_N;
				it_v->second.bxset_alt_T = v.bxset_alt_T;
			
				//it_v->second.reGenotype(); // recompute genotype
			}
		}
	}
	else { 
		DB.insert(pair<string,Variant_t>(key,v));
	}
}

// select variant supported by most windows
void VariantDB_t::selectVar() {
	// a ratio when somatic status dose not match && have same similar_variants_count .
	// float cover_ratio = 1.5;
	for (map<string, Variant_t>::iterator it_v = DB.begin(); it_v != DB.end(); ++it_v) {
		string key = it_v->first;
		string isSomatic = key.substr(0, 1);
		string signature = key.substr(1);
		string isSomatic2 = (isSomatic=="1") ? "0": "1";
		string key2 = isSomatic2 + signature;
		map<string, Variant_t>::iterator it_v2 = DB.find(key2);
		if (it_v2 == DB.end())
			continue;
		bool select_success = 1;
		// map<string, Variant_t>::iterator select_v; // = DB.begin();
		// map<string, Variant_t>::iterator discard_v; // = DB.begin();
		string select_reason;
		Variant_t select_v = it_v->second;
		Variant_t discard_v = it_v2->second;
		
		int svc = it_v->second.similar_variants_count;
		int svc2 = it_v2->second.similar_variants_count;
		if (svc > svc2) {
			// Variant_t select_v = it_v->second;
			// Variant_t discard_v = it_v2->second;
			select_reason = "Select by similar-variants-count. Select:" + itos(svc) + ", Discard:" + itos(svc2);
			DB.erase(it_v2);
		}
		else if (svc < svc2) {
			Variant_t tmp = select_v;
			select_v = discard_v;
			discard_v = tmp;
			select_reason = "Select by similar-variants-count. Select:" + itos(svc2) + ", Discard:" + itos(svc);
			DB.erase(it_v++);
		}
		else {
			int ref_cov_normal_1 = it_v->second.ref_cov_normal_fwd + it_v->second.ref_cov_normal_rev;
			int ref_cov_tumor_1 = it_v->second.ref_cov_tumor_fwd + it_v->second.ref_cov_tumor_rev;
			int alt_cov_normal_1 = it_v->second.alt_cov_normal_fwd + it_v->second.alt_cov_normal_rev;
			int alt_cov_tumor_1 = it_v->second.alt_cov_tumor_fwd + it_v->second.alt_cov_tumor_rev;
			
			int ref_cov_normal_2 = it_v2->second.ref_cov_normal_fwd + it_v2->second.ref_cov_normal_rev;
			int ref_cov_tumor_2 = it_v2->second.ref_cov_tumor_fwd + it_v2->second.ref_cov_tumor_rev;
			int alt_cov_normal_2 = it_v2->second.alt_cov_normal_fwd + it_v2->second.alt_cov_normal_rev;
			int alt_cov_tumor_2 = it_v2->second.alt_cov_tumor_fwd + it_v2->second.alt_cov_tumor_rev;
			
			int tot_cov_1 = ref_cov_normal_1 + ref_cov_tumor_1 + alt_cov_normal_1 + alt_cov_tumor_1;
			int tot_cov_2 = ref_cov_normal_2 + ref_cov_tumor_2 + alt_cov_normal_2 + alt_cov_tumor_2;
			select_reason = "Same similar-variants-count("+ itos(svc) + ", " + itos(svc2) + "). ";
			if (tot_cov_1 > tot_cov_2) {
				// Variant_t select_v = it_v->second;
				// Variant_t discard_v = it_v2->second;
				select_reason += "Select by  total-coverage. Select:" + itos(tot_cov_1) + ", Discard:" + itos(tot_cov_2);
				DB.erase(it_v2);
			}
			else if (tot_cov_1 < tot_cov_2) {
				// Variant_t select_v = it_v2->second;
				// Variant_t discard_v = it_v->second;
				Variant_t tmp = select_v;
				select_v = discard_v;
				discard_v = tmp;
				select_reason += "Select by  total-coverage. Select:" + itos(tot_cov_2) + ", Discard:" + itos(tot_cov_1);
				DB.erase(it_v++);
			}
			else {
				// retain all
				select_reason = "Have same similar-variants-count(" + itos(svc) + ", " + itos(svc2) + ") and total-coverage(" + itos(tot_cov_1) + ", " + itos(tot_cov_2) + ")";
				select_success = 0;
			}
		}
		if (select_success)	{
			cerr << "***** Variant select *****" << endl;
			cerr << "Signature: " << select_v.getSignature() << endl;
			cerr << "Select variant: " << select_v.printVcfWithoutFilters() << endl;
			cerr << "               Somatic Status: " << itos(select_v.isSomatic) << endl;
			cerr << "Discard variant: " << discard_v.printVcfWithoutFilters() << endl;
			cerr << "               Somatic Status: " << itos(discard_v.isSomatic) << endl;
			cerr << "Reason: " << select_reason << endl;
			cerr << "***************************" << endl;
		}
		else {
			cerr << "** Variant select failed **" << endl;
			cerr << "Signature: " << it_v->second.getSignature() << endl;
			cerr << "Retain both variants: " << endl;
			cerr << "First variant: " << it_v->second.printVcfWithoutFilters() << endl;
			cerr << "               Somatic Status: " << itos(it_v->second.isSomatic) << endl;
			cerr << "Second variant: " << it_v2->second.printVcfWithoutFilters() << endl;
			cerr << "               Somatic Status: " << itos(it_v2->second.isSomatic) << endl;
			cerr << "Reason: " << select_reason;
			cerr << "***************************" << endl;
		}
		// avoid it_v out of range
		if (it_v == DB.end())
			break;
		
		
	}
}


void VariantDB_t::printHeader(const string version, const string reference, char * date, Filters &fs, string &sample_name_N, string &sample_name_T) {
	
    stringstream hdr;
	
	hdr << "##fileformat=VCFv4.2\n"
			"##fileDate=" << date << ""
			"##source=lancet " << version << "\n"
			"##cmdline=" << command_line << "\n"
			"##reference=" << reference << "\n"
			"##INFO=<ID=FETS,Number=1,Type=Float,Description=\"Phred-scaled p-value of the Fisher's exact test for tumor-normal allele counts\">\n"
			"##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic mutation\">\n"
			"##INFO=<ID=SHARED,Number=0,Type=Flag,Description=\"Shared mutation betweem tumor and normal\">\n"
			"##INFO=<ID=NORMAL,Number=0,Type=Flag,Description=\"Mutation present only in the normal\">\n"
			"##INFO=<ID=NONE,Number=0,Type=Flag,Description=\"Mutation not supported by data\">\n"
			"##INFO=<ID=KMERSIZE,Number=1,Type=Integer,Description=\"K-mer size used to assemble the locus\">\n"
			"##INFO=<ID=SB,Number=1,Type=Float,Description=\"Strand bias score: phred-scaled p-value of the Fisher's exact test for the forward/reverse read counts in the tumor\">\n"
			"##INFO=<ID=MS,Number=1,Type=String,Description=\"Microsatellite mutation (format: #LEN#MOTIF)\">\n"
			"##INFO=<ID=LEN,Number=1,Type=Integer,Description=\"Variant size in base pairs\">\n"
			"##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Variant type (snv, del, ins, complex)\">\n";
	
	if(LR_MODE)	{
		hdr << "##INFO=<ID=HPS,Number=1,Type=Float,Description=\"Haplotype score for the T/N pair: phred-scaled p-value of the Fisher's exact test of the total counts of the two haplotype in the tumor-normal pair\">\n"
			   "##INFO=<ID=HPSN,Number=1,Type=Float,Description=\"Normal haplotype score: phred-scaled p-value of the Fisher's exact test for ref/alt haplotype counts in the normal\">\n"
			   "##INFO=<ID=HPST,Number=1,Type=Float,Description=\"Tumor haplotype score: phred-scaled p-value of the Fisher's exact test for ref/alt haplotype counts in the tumor\">\n";
	}
				
	hdr <<	"##FILTER=<ID=LowCovNormal,Description=\"Low coverage in the normal (<" << fs.minCovNormal << ")\">\n"
			"##FILTER=<ID=HighCovNormal,Description=\"High coverage in the normal (>" << fs.maxCovNormal << ")\">\n"
			"##FILTER=<ID=LowCovTumor,Description=\"Low coverage in the tumor (<" << fs.minCovTumor << ")\">\n"
			"##FILTER=<ID=HighCovTumor,Description=\"High coverage in the tumor (>" << fs.maxCovTumor << ")\">\n"
			"##FILTER=<ID=LowVafTumor,Description=\"Low variant allele frequency in the tumor (<" << fs.minVafTumor << ")\">\n"
			"##FILTER=<ID=HighVafNormal,Description=\"High variant allele frequency in the normal (>" << fs.maxVafNormal << ")\">\n"
			"##FILTER=<ID=LowAltCntTumor,Description=\"Low alternative allele count in the tumor (<" << fs.minAltCntTumor << ")\">\n"
			"##FILTER=<ID=HighAltCntNormal,Description=\"High alternative allele count in the normal (>" << fs.maxAltCntNormal << ")\">\n"
			"##FILTER=<ID=LowFisherScore,Description=\"Low Fisher's exact test score for tumor-normal allele counts (<" << fs.minPhredFisher << ")\">\n"
			"##FILTER=<ID=LowFisherSTR,Description=\"Low Fisher's exact test score for tumor-normal STR allele counts (<" << fs.minPhredFisherSTR << ")\">\n"
			"##FILTER=<ID=StrandBias,Description=\"Strand bias: # of non-reference reads in either forward or reverse strand below threshold (<" << fs.minStrandBias << ")\">\n"
			"##FILTER=<ID=STR,Description=\"Microsatellite mutation\">\n";
	
	if(LR_MODE)	{
		hdr << "##FILTER=<ID=MultiHP,Description=\"Supporting reads from multiple haplotypes based on linked-reads analysis\">\n";
	}	
				
	hdr << 	"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
			"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
			"##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allele depth: # of supporting ref,alt reads at the site\">\n"
			"##FORMAT=<ID=SR,Number=.,Type=Integer,Description=\"Strand counts for ref: # of supporting forward,reverse reads for reference allele\">\n"
			"##FORMAT=<ID=SA,Number=.,Type=Integer,Description=\"Strand counts for alt: # of supporting forward,reverse reads for alterantive allele\">\n";
	
	if(LR_MODE)	{
		hdr << "##FORMAT=<ID=BX,Number=.,Type=String,Description=\"Barcodes supporting ref and alt alleles\">\n"
		       "##FORMAT=<ID=HPR,Number=.,Type=Integer,Description=\"Haplotype counts for ref: # of reads supporting reference allele in haplotype 1, 2, and 0 respectively (0 = unassigned)\">\n"
			   "##FORMAT=<ID=HPA,Number=.,Type=Integer,Description=\"Haplotype counts for alt: # of reads supporting alternative allele in haplotype 1, 2, and 0 respectively (0 = unassigned)\">\n";
	}
	
	hdr << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << sample_name_N << "\t" << sample_name_T << "\n";
	
	cout << hdr.str();
}

// print variant in VCF format
void VariantDB_t::printToVCF(const string version, const string reference, char * date, Filters &fs, string &sample_name_N, string &sample_name_T) {
	
	cerr << "Export variants to VCF file" << endl;
	
	printHeader(version,reference,date,fs,sample_name_N,sample_name_T);
	
	// dump map content to vector for custom sorting
	vector< pair<string,Variant_t> > myVec(DB.begin(), DB.end());
	// sort based on chromosome location
	sort(myVec.begin(),myVec.end(),byPos());

	vector< pair<string,Variant_t> >::iterator it;
	for (it=myVec.begin(); it!=myVec.end(); ++it) {
		//cerr << it->first << "\t";
		
		//cerr << "Size of " << it->first << " : " << sizeof(it->second) << endl;
		
		//it->second.printSizeOfElements();
				
		//string pos = (it->second).getPosition();
	    //unordered_map<string,int>::iterator itp = nCNT.find(pos);
		//if (itp == nCNT.end()) { // print variant if no muations in the normal at locus
			cout << it->second.printVCF(filters);
		//}
	}	
}

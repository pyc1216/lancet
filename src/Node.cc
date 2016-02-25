#include "Node.hh"

/****************************************************************************
** Node.cc
**
** Implementation of a node of the de Bruijn graph
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

// isTandem
//////////////////////////////////////////////////////////////

bool Node_t::isTandem()
{
	for (unsigned int i = 0; i < edges_m.size(); i++)
	{
		if (edges_m[i].nodeid_m == nodeid_m)
		{
			return true;
		}
	}

	return false;
}


// addEdge
//////////////////////////////////////////////////////////////

void Node_t::addEdge(Mer_t nodeid, Edgedir_t dir, ReadId_t readid)
{
	if (readid != -1)
	{
		reads_m.insert(readid);
	}

	int edgeid = -1;

	for (unsigned int i = 0; i < edges_m.size(); i++)
	{
		if ((edges_m[i].nodeid_m == nodeid) &&
			(edges_m[i].dir_m == dir))
		{
			if (readid != -1)
			{
				edges_m[i].readids_m.push_back(readid);
			}

			edgeid = i;
			break;
		}
	}

	if (edgeid == -1)
	{
		Edge_t ne(nodeid, dir);

		if (readid != -1)
		{
			ne.readids_m.push_back(readid);
		}
		edges_m.push_back(ne);
	}
}


// updateEdge
//////////////////////////////////////////////////////////////

void Node_t::updateEdge(const Mer_t & oldid, Edgedir_t olddir, 
	const Mer_t & newid, Edgedir_t newdir)
{
	bool found = false;
	for (unsigned int i = 0; i < edges_m.size(); i++)
	{
		if (edges_m[i].nodeid_m == oldid &&
			edges_m[i].dir_m == olddir)
		{
			edges_m[i].nodeid_m = newid;
			edges_m[i].dir_m = newdir;
			return;
		}
	}

	if (!found)
	{
		cerr << "ERROR updating " << nodeid_m << ": didn't find: " << Edge_t::toString(olddir) << ":" << oldid << endl;
		cerr << "wanted to replace with " << Edge_t::toString(newdir) << ":" << newid << endl;
		print(cerr);
	}

	assert(found);
}

// removeEdge
//////////////////////////////////////////////////////////////

void Node_t::removeEdge(const Mer_t & nodeid, Edgedir_t dir)
{
	bool found = false;

	for (unsigned int i = 0; i < edges_m.size(); i++)
	{
		if (edges_m[i].nodeid_m == nodeid &&
			edges_m[i].dir_m == dir)
		{
			edges_m.erase(edges_m.begin()+i);
			return;
		}
	}

	if (!found)
	{
		cerr << "ERROR removing " << nodeid_m << ": didn't find: " << dir << ":" << nodeid << endl;
	}

	assert(found);
}


// getBuddy: If node has a single edge in the dir direction, return edge id, else -1
//////////////////////////////////////////////////////////////

int Node_t::getBuddy(Ori_t dir)
{
	int retval = -1;

	//if (isRef()) { return retval; }
	if (isSpecial()) { return retval; }

	for (unsigned int i = 0; i < edges_m.size(); i++)
	{
		if (edges_m[i].isDir(dir))
		{
			//cerr << edges_m[i] << endl;
			if (retval != -1)
			{
				return -1;
			}

			retval = i;
		}
	}

	// self loop
	if (retval != -1)
	{
		if (edges_m[retval].nodeid_m == nodeid_m)
		{
			return -1;
		}
	}

	return retval;
}

// markRef
//////////////////////////////////////////////////////////////

int Node_t::markRef(Ref_t * ref, int K)
{
	// handle special source/sink nodes
	if (isSource() || isSink()) { return 1; }

	CanonicalMer_t cmer;

	touchRef_m = false;

	for (unsigned int i = 0; i < str_m.length()-K+1; i++)
	{
		cmer.set(str_m.substr(i, K));

		if (ref->hasMer(cmer.mer_m))
		{
			touchRef_m = true;
			//isRef_m = true;
			return 1;
		}
	}

	return 0;
}


// degree
//////////////////////////////////////////////////////////////

int Node_t::degree(Ori_t dir)
{
	int retval = 0;
	for (unsigned int i = 0; i < edges_m.size(); i++)
	{
		if (edges_m[i].isDir(dir))
		{
			retval++;
		}
	}

	return retval;
}

// print
//////////////////////////////////////////////////////////////

ostream & Node_t::print(ostream & out) const
{
	out << nodeid_m;

	out << "\t*s\t" << str_m;
	out << "\t*c\t" << cov_tmr_m_fwd;
	out << "\t*c\t" << cov_nml_m_fwd;
	out << "\t*c\t" << cov_tmr_m_rev;
	out << "\t*c\t" << cov_nml_m_rev;
	out << "\t*r\t" << isRef();

	for (unsigned int i = 0; i < edges_m.size(); i++)
	{
		out << "\t" << edges_m[i];
	}

	return out;
}

// strlen
//////////////////////////////////////////////////////////////

int Node_t::strlen() const
{
	//if (isRef()) { return 0; }
	if (isSpecial()) { return 0; }
	return str_m.length();
}


// addReadStart
//////////////////////////////////////////////////////////////

void Node_t::addReadStart(ReadId_t readid, int nodeoffset, int trim5, Ori_t ori)
{
	readstarts_m.push_back(ReadStart_t(readid, nodeoffset, trim5, ori));
}

// revreads
//////////////////////////////////////////////////////////////

void Node_t::revreads()
{
	//assert(!isRef());
	assert(!isSpecial());

	int len = strlen();

	for (unsigned int i = 0; i < readstarts_m.size(); i++)
	{
		ReadStart_t & rs = readstarts_m[i];
		rs.nodeoffset_m = len - 1 - rs.nodeoffset_m;
		rs.ori_m = Edge_t::flipdir(rs.ori_m);
	}
}

// sortReadStarts
//////////////////////////////////////////////////////////////

void Node_t::sortReadStarts()
{
	sort(readstarts_m.begin(), readstarts_m.end(), ReadStart_t::cmpstarts);
}


// addContigLink
//////////////////////////////////////////////////////////////

void Node_t::addContigLink(Mer_t contigid, ReadId_t rid)
{
	ContigLinkMap_t::iterator mi = contiglinks_m.find(contigid);

	if (mi == contiglinks_m.end())
	{
		mi = contiglinks_m.insert(make_pair(contigid, new ContigLinkList_t)).first;
	}

	mi->second->addLink(rid);
}

// cntReadCode
//////////////////////////////////////////////////////////////

int Node_t::cntReadCode(char code)
{
	int retval = 0;

	unordered_set<ReadId_t>::iterator si;
	for (si = reads_m.begin();
	si != reads_m.end();
	si++)
	{
		if ( (*si) < 0 || (*si) > (int)readid2info->size() ) { continue; } // skip over illegal values
		if (readid2info->at(*si).code_m == code)
		{
			retval++;
		}
	}

	return retval;
}

// updateCovDistrTmr
// updated the coverage distribution along the node string
//////////////////////////////////////////////////////////////
void Node_t::updateCovDistr(int c, unsigned int strand, char sample) 
{
	vector<int> * cov_distr = NULL;
	
	if(sample == 'T') {
		if(strand == FWD) { cov_distr = &cov_distr_tmr_fwd; }
		if(strand == REV) { cov_distr = &cov_distr_tmr_rev; }
	}
	else if(sample == 'N') {
		if(strand == FWD) { cov_distr = &cov_distr_nml_fwd; }
		if(strand == REV) { cov_distr = &cov_distr_nml_rev; }
	}
	
	for (unsigned int i = 0; i < cov_distr->size(); i++) {
		(*cov_distr)[i] = c; 
	}
}

// updateCovDistrMinQv
// updated the coverage distribution along the node with 
// base-quality value greater than MIN_QUAL
//////////////////////////////////////////////////////////////
void Node_t::updateCovDistrMinQV(const string & qv, unsigned int strand, char sample) {
	
	vector<int> * cov_distr = NULL;

	if(sample == 'T') {
		if(strand == FWD) { cov_distr = &cov_distr_tmr_minqv_fwd; }
		if(strand == REV) { cov_distr = &cov_distr_tmr_minqv_rev; }
	}
	else if(sample == 'N') {
		if(strand == FWD) { cov_distr = &cov_distr_nml_minqv_fwd; }
		if(strand == REV) { cov_distr = &cov_distr_nml_minqv_rev; }
	}

	unsigned int i = 0;
	for ( string::const_iterator it=qv.begin(); it!=qv.end(); ++it) {
		if(*it >= MIN_QUAL) { (*cov_distr)[i]++; }
		i++;
	}
}

// avgCovDistr
// average coverage of non-zoero elements
//////////////////////////////////////////////////////////////
int Node_t::avgCovDistr(char sample)
{
	vector<int> cov_distr_fwd;
	vector<int> cov_distr_rev;
	
	if(sample == 'T') { 
		cov_distr_fwd = cov_distr_tmr_fwd; 
		cov_distr_rev = cov_distr_tmr_rev;
	}
	else if(sample == 'N') { 
		cov_distr_fwd = cov_distr_nml_fwd; 
		cov_distr_rev = cov_distr_nml_rev; 
	}
	
	int sum = 0;
	int cnt = 0;
	for (unsigned int i = 0; i < cov_distr_fwd.size(); i++)
	{
		int totcov = cov_distr_fwd[i] + cov_distr_rev[i];
		if(totcov !=0) {
			sum += cov_distr_fwd[i] + cov_distr_rev[i];
			cnt++;
		}
	}
	//cerr << "SPANNER: (" << sum << "," << cnt << ")" << endl; 

	int ans = 0;
	if (sum>0) { ans = floor(float(sum)/float(cnt)); }
	
	return ans;
}

// revCovDistr
// reverse the coverage distribution along the node string
//////////////////////////////////////////////////////////////
void Node_t::revCovDistr() 
{
	int i=0;
	int j=cov_distr_tmr_fwd.size()-1;
	while(i<j){
		
		swap(cov_distr_tmr_fwd[i], cov_distr_tmr_fwd[j]);
		swap(cov_distr_tmr_rev[i], cov_distr_tmr_rev[j]);
		swap(cov_distr_nml_fwd[i], cov_distr_nml_fwd[j]);
		swap(cov_distr_nml_rev[i], cov_distr_nml_rev[j]);
		
		swap(cov_distr_tmr_minqv_fwd[i], cov_distr_tmr_minqv_fwd[j]);
		swap(cov_distr_tmr_minqv_rev[i], cov_distr_tmr_minqv_rev[j]);
		swap(cov_distr_nml_minqv_fwd[i], cov_distr_nml_minqv_fwd[j]);
		swap(cov_distr_nml_minqv_rev[i], cov_distr_nml_minqv_rev[j]);
		
		i++;j--;
	}
}

// minCov
// retunr the minimum coverage along the node
//////////////////////////////////////////////////////////////
int Node_t::minNon0Cov(char sample) 
{
	
	vector<int> cov_distr_fwd;
	vector<int> cov_distr_rev;

	if(sample == 'T') { 
		cov_distr_fwd = cov_distr_tmr_fwd; 
		cov_distr_rev = cov_distr_tmr_rev; 
	}
	else if(sample == 'N') { 
		cov_distr_fwd = cov_distr_nml_fwd; 
		cov_distr_rev = cov_distr_nml_rev; 	
	}
	
	int min = 10000000;
	for (unsigned int i = 0; i < cov_distr_fwd.size(); i++)
	{
		int totcov = cov_distr_fwd[i] + cov_distr_rev[i];
		if( (totcov > 0) && (totcov < min) ) { min = totcov; }
	}
	
	return min;
}


// minCov
// retunr the minimum coverage along the node
//////////////////////////////////////////////////////////////
int Node_t::minCov() 
{
	int min = 10000000;
	for (unsigned int i = 0; i < cov_distr_tmr_fwd.size(); i++)
	{
		int totcov = cov_distr_tmr_fwd[i] + cov_distr_tmr_rev[i] + cov_distr_nml_fwd[i] + cov_distr_nml_rev[i];;
		if(totcov < min) { min = totcov; } 
	}
	return min;
}

// readOverlaps
// return the size of the reads overlap
//////////////////////////////////////////////////////////////
int Node_t::readOverlaps(const Node_t & other)
{
	int retval = 0;
	unordered_set<ReadId_t>::const_iterator it;
	
	for (auto it = other.reads_m.begin(); it != other.reads_m.end(); it++) {

		if (reads_m.find(*it) != reads_m.end())
		{
			retval++;
		}
	}

	return retval;
}

// return tumor coverage on the input strand
//////////////////////////////////////////////////////////////
float Node_t::getTmrCov(unsigned int strand) { 
	float ans = 0;
	
	if(strand == FWD) { ans = cov_tmr_m_fwd; } 
	else if(strand == REV) { ans = cov_tmr_m_rev; } 
	
	return ans;
}

// return normal coverage on the input strand
//////////////////////////////////////////////////////////////
float Node_t::getNmlCov(unsigned int strand) { 
	float ans = 0;
	
	if(strand == FWD) { ans = cov_nml_m_fwd; } 
	else if(strand == REV) { ans = cov_nml_m_rev; } 
	
	return ans;
}

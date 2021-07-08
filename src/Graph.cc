#include <vector>

#include "Graph.hh"

/****************************************************************************
** Graph.cc
**
** Class for representing and storing a bi-directed de Bruijn graph  
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

// clear edge flags
//////////////////////////////////////////////////////////////
void Graph_t::clear(bool flag)
{
	if(flag) {
		//readid2info.clear();
		vector<ReadInfo_t>().swap(readid2info); 
		is_ref_added = false; // reference reads was in readid2info and removed
	}
	totalreadbp_m = 0;

	MerTable_t::iterator mi;
	for (mi = nodes_m.begin(); mi != nodes_m.end(); ++mi)
	{
		delete mi->second;
	}
	nodes_m.clear();
	unordered_map<Mer_t, Node_t *>().swap(nodes_m);	
	//sparse_hash_map<Mer_t, Node_t *>().swap(nodes_m);
	//hopscotch_map<Mer_t, Node_t *, hash<Mer_t>, equal_to<Mer_t>, allocator<pair<Mer_t,Node_t *>>, 30, true>().swap(nodes_m);

	bx_table_tmr.clear(); unordered_map<Mer_t,set<string>>().swap(bx_table_tmr);
	bx_table_nml.clear(); unordered_map<Mer_t,set<string>>().swap(bx_table_nml);

	source_m = NULL;
	sink_m = NULL;

	if (ref_m && flag == true)
	{
		ref_m->clear();
		delete ref_m;
		ref_m = NULL;
	}	
}


// addBX
//////////////////////////////////////////////////////////////
void Graph_t::addBX(const string & bx, Mer_t & mer, int sample) {
		
	unordered_map<Mer_t,set<string>> * map = NULL;
	
	if(sample == TMR) { map = &bx_table_tmr; }
	if(sample == NML) { map = &bx_table_nml; } 
	
	auto got = map->find(mer);
	if ( got != map->end() ) {
		(got->second).insert(bx);
	}
	else { // new BX
		(*map)[mer].insert(bx);
	}
}

// getBXsetAt
//////////////////////////////////////////////////////////////
string Graph_t::getBXsetAt(int start, int end, string & seq, int sample) {
	
	CanonicalMer_t cmer;
	string result;
	set<string> bxset;
	unordered_map<Mer_t,set<string>> * map;
	
	if(sample == TMR) { map = &bx_table_tmr; }
	else if(sample == NML) { map = &bx_table_nml; }
	else { cerr << "Error: unrecognized sample " << sample << endl; }
		
	for (int i=start; i<=end; i++) {
		
		cmer.set(seq.substr(i,K));
				
		auto it = map->find(cmer.mer_m);
		if(it != map->end()) {
			bxset.insert((it->second).begin(),(it->second).end());
		}
	}
	
	for ( auto itv = bxset.begin(); itv != bxset.end(); ++itv ) { 
		if (next(itv) == bxset.end()) { result += *itv; }
		else { result += *itv + ";"; }
	}
		
	if(result == "") { result = "."; }
	
	//cerr << "BX set: " << result << endl;
	
	return result;
}

// loadSequence
//////////////////////////////////////////////////////////////

void Graph_t::loadSequence(int readid, const string & seq, const string & qv, bool isRef, int trim5, unsigned int strand)
{	
	if (!isRef)
	{
		totalreadbp_m += seq.length();
	}
		
	CanonicalMer_t uc;
	CanonicalMer_t vc;
	string uc_qv;
	string vc_qv;
	
	Node_t * unode = NULL;
	Node_t * vnode = NULL;
		
	MerTable_t::iterator ui;
	MerTable_t::iterator vi;

	unordered_set<Mer_t> readmers;

	int sample = readid2info[readid].label_m; // TMR or NML

	int end = seq.length() - K;
	int offset = 0;
	for (; offset < end; ++offset)
	{			
		if (offset == 0) {
			uc.set(seq.substr(offset,   K));
			vc.set(seq.substr(offset+1, K));	
			uc_qv = qv.substr(offset,   K);
			vc_qv = qv.substr(offset+1, K);
			if (uc.ori_m == R) { reverse(uc_qv); }
			if (vc.ori_m == R) { reverse(vc_qv); }
		}
		else {
			uc = vc; 
			uc_qv = vc_qv; 
			vc.set(seq.substr(offset+1, K));
			vc_qv = qv.substr(offset+1, K);
			if (vc.ori_m == R) { reverse(vc_qv); }
		}

		//cerr << readid << "\t" << offset << "\t" << uc << "\t" << vc << endl;
		
		if (offset == 0) {
			ui = nodes_m.find(uc.mer_m);
			vi = nodes_m.find(vc.mer_m);
		}
		else {
			ui = vi;
			vi = nodes_m.find(vc.mer_m);
		}
		
		bool ui_found = false;
		bool vi_found = false;
		
		if (ui != nodes_m.end()) { ui_found = true; unode = ui->second; }
		if (vi != nodes_m.end()) { vi_found = true; vnode = vi->second; }

		//if (ui == nodes_m.end())
		if(!ui_found)
		{
			ui = nodes_m.insert(make_pair(uc.mer_m, new Node_t(uc.mer_m))).first;
			unode = ui->second;
			unode->setMinQV(MIN_QUAL_CALL);
			unode->setK(K);			
			//ui->second->setRead2InfoList(&readid2info);
		}
		//else { unode = ui->second; }
				
		//if (vi == nodes_m.end())
		if(!vi_found)
		{
			vi = nodes_m.insert(make_pair(vc.mer_m, new Node_t(vc.mer_m))).first;
			vnode = vi->second;
			vnode->setMinQV(MIN_QUAL_CALL);
			vnode->setK(K);
			//vi->second->setRead2InfoList(&readid2info);
		}
		//else { vnode = vi->second; }

		// always set node label for normal reads 
		// even if kmer has low quality bases 
		if(sample == NML) {
			unode->setIsNormal();
			vnode->setIsNormal();
			unode->updateCovStatus('N');
			vnode->updateCovStatus('N');
		}
						
		if( (seqAboveQual(uc_qv,MIN_QUAL_CALL) && seqAboveQual(vc_qv,MIN_QUAL_CALL)) ) {
			
			// set node label
			if(sample == TMR) {
				unode->setIsTumor();
				vnode->setIsTumor();
				unode->updateCovStatus('T');
				vnode->updateCovStatus('T');
			}
			/*
			else if(readid2info[readid].label_m == NML) {
				ui->second->setIsNormal();
				vi->second->setIsNormal();
				ui->second->updateCovStatus('N');
				vi->second->updateCovStatus('N');
			}
			*/
		}
		//ui->second->appendRefFlag(isRef);
		//vi->second->appendRefFlag(isRef);
				
		// add mate name info to the nodes
		// (used to check for overlapping mates)
		unode->addMateName(readid2info[readid].readname_m, readid2info[readid].mate_order_m);
		vnode->addMateName(readid2info[readid].readname_m, readid2info[readid].mate_order_m);
		
		bool isOvlMate = false;
		//bool bxovl_u = false;
		//bool bxovl_v = false;
				 
		if(LR_MODE) {
										
			if (offset == 0) {
								
				if (readid2info[readid].BX != "null") { // skip over null barcodes
					addBX(readid2info[readid].BX, uc.mer_m, sample);
					ref_m->addBX(readid2info[readid].BX, uc.mer_m, sample);
				}
				
				if(!(unode->hasBX(readid2info[readid].BX, sample))) { // update only if BX not already present in this node
					unode->addBX(readid2info[readid].BX, strand, sample);
					unode->addHP(readid2info[readid].HP, sample);
				}
			}	
			
			if (readid2info[readid].BX != "null") { // skip over null barcodes
				addBX(readid2info[readid].BX, vc.mer_m, sample);
				ref_m->addBX(readid2info[readid].BX, vc.mer_m, sample);
			}
			
			if(!(vnode->hasBX(readid2info[readid].BX, sample))) { // update only if BX not already present in this node
				vnode->addBX(readid2info[readid].BX, strand, sample);
				vnode->addHP(readid2info[readid].HP, sample);
			}
		}
		
		if (!isRef)
		{		
			if (offset == 0) 
			{ 
				isOvlMate = (unode->hasOverlappingMate(readid2info[readid].readname_m, readid2info[readid].mate_order_m)); //kmer from overlapping mates
				
				//if(isOvlMate) { cerr << "Overlapping mates for fragment:" << readid2info[readid].readname_m << endl; }
									
				if( !isOvlMate) { // do not update coverage for overlapping mates
										
					unode->incCov(strand, sample);
					
					if(LR_MODE) { 
						unode->updateCovDistr(unode->BXcnt(strand,sample), uc_qv, strand, sample); 
						unode->updateHPCovDistr(unode->HPcnt(0,sample), unode->HPcnt(1,sample), unode->HPcnt(2,sample), uc_qv, sample); 
						ref_m->updateCoverage(uc.mer_m, unode->BXcnt(strand,sample), strand, sample); // update reference k-mer coverage
						ref_m->updateHPCoverage(uc.mer_m, unode->HPcnt(0,sample), unode->HPcnt(1,sample), unode->HPcnt(2,sample), sample); 
					}
					else { 
						unode->updateCovDistr((int)(unode->getCov(strand,sample)), uc_qv, strand, sample); 
						ref_m->updateCoverage(uc.mer_m, (int)(unode->getCov(strand,sample)), strand, sample); // update reference k-mer coverage
					}

					if (uc.ori_m == F)
					{
						unode->addReadStart(readid, 0, trim5, uc.ori_m);
					}
					else
					{
						unode->addReadStart(readid, K-1, trim5, uc.ori_m);
					}
				}
			}
			
			isOvlMate = (vnode->hasOverlappingMate(readid2info[readid].readname_m, readid2info[readid].mate_order_m));

			//if(isOvlMate) { cerr << "Overlapping mates for fragment:" << readid2info[readid].readname_m << endl; }

			if( !isOvlMate ) { // do not update coverage for overlapping mates
								
					vnode->incCov(strand, sample);
					
					if(LR_MODE) { 
						vnode->updateCovDistr(vnode->BXcnt(strand,sample), vc_qv, strand, sample); 
						vnode->updateHPCovDistr(vnode->HPcnt(0,sample), vnode->HPcnt(1,sample), vnode->HPcnt(2,sample), vc_qv, sample); 
						ref_m->updateCoverage(vc.mer_m, vnode->BXcnt(strand,sample), strand, sample); // update reference k-mer coverage
						ref_m->updateHPCoverage(vc.mer_m, vnode->HPcnt(0,sample), vnode->HPcnt(1,sample), vnode->HPcnt(2,sample), sample);
					}
					else { 
						vnode->updateCovDistr((int)(vnode->getCov(strand,sample)), vc_qv, strand, sample); 
						ref_m->updateCoverage(vc.mer_m, (int)(vnode->getCov(strand,sample)), strand, sample); // update reference k-mer coverage
					}
			}
		}
		
		Edgedir_t fdir = FF;
		Edgedir_t rdir = FF;

		if      (uc.ori_m == F && vc.ori_m == F) { fdir = FF; rdir = RR; }
		else if (uc.ori_m == F && vc.ori_m == R) { fdir = FR; rdir = FR; }
		else if (uc.ori_m == R && vc.ori_m == F) { fdir = RF; rdir = RF; }
		else if (uc.ori_m == R && vc.ori_m == R) { fdir = RR; rdir = FF; }

		readmers.insert(uc.mer_m);

		if (readmers.find(vc.mer_m) != readmers.end())
		{
			if (VERBOSE) { cerr << "cycle detected in read " << readid << " offset: " << offset << " : " << seq << endl; }

			if (readid > -1)
			{
				if (readCycles == 0)
				{
					if(verbose) { cerr << "WARNING: Cycles detected in the reads" << endl << endl; }
				}

				++readCycles;
			}

			//readid = -1;
		}
		unode->addEdge(vc.mer_m, fdir, readid);
		vnode->addEdge(uc.mer_m, rdir, readid);
	}
}


// trim
//////////////////////////////////////////////////////////////

void Graph_t::trim(int readid, const string & seq, const string & qv, bool isRef)
{
	int trim3 = 0;
	int trim5 = 0;
	int len = seq.length();
	string cseq = seq;

	for (int i = 0; i < len; ++i) { cseq[i] = toupper(cseq[i]); }

	while ((!isDNA(seq[trim5]) || (qv[trim5] < MIN_QUAL_TRIM)) && (trim5 < len)) { ++trim5; }

	if (trim5 < len) {
		while ((!isDNA(seq[len-1-trim3]) || (qv[len-1-trim3] < MIN_QUAL_TRIM)) && (trim3 < len)) { ++trim3; }

		readid2info[readid].isjunk = false;
		
		for (int i = trim5; i < len-trim3; ++i)
		{
			if (!isDNA(seq[i])) {
				// skip the junk
				readid2info[readid].isjunk = true;
				break;
			}
		}
	}
	else { readid2info[readid].isjunk = true; }
	
	readid2info[readid].trm5 = trim5;
	readid2info[readid].trm3 = trim3;
}

// trimAndLoad
//////////////////////////////////////////////////////////////

void Graph_t::trimAndLoad(int readid, const string & seq, const string & qv, bool isRef, unsigned int strand)
{
	int len = seq.length();
	string cseq = seq;

	for (int i = 0; i < len; ++i) { cseq[i] = toupper(cseq[i]); }

	int trim5 = 0;
	while ((!isDNA(seq[trim5]) || (qv[trim5] < MIN_QUAL_TRIM)) && (trim5 < len)) { ++trim5; }

	if (trim5 < len)
	{
		int trim3 = 0;
		while ((!isDNA(seq[len-1-trim3]) || (qv[len-1-trim3] < MIN_QUAL_TRIM)) && (trim3 < len)) { ++trim3; }

		bool cleanRead = true;

		for (int i = trim5; i < len-trim3; ++i)
		{
			if (!isDNA(seq[i]))
			{
				// skip the junk
				cleanRead = false;
				break;
			}
		}

		if (cleanRead)
		{
			if (trim5 || trim3) { cseq = seq.substr(trim5, len-trim5-trim3); }
			loadSequence(readid, cseq, qv, isRef, trim5, strand);
		}
	}
}

int Graph_t::countBastardReads()
{
	int bastards = 0;

	for (unsigned int i = 0; i < readid2info.size(); ++i)
	{
		if (readid2info[i].code_m == CODE_BASTARD)
		{
			++bastards;
		}
	}

	return bastards;
}

int Graph_t::countMappedReads()
{
	int mapped = 0;

	for (unsigned int i = 0; i < readid2info.size(); ++i)
	{
		if (readid2info[i].code_m == CODE_MAPPED)
		{
			++mapped;
		}
	}

	return mapped;
}

// addRead
////////////////////////////////////////////////////////////////

ReadId_t Graph_t::addRead(const string & set, const string & readname, const string & seq, const string & qv, char code, int label, unsigned int strand, int mate_order, const string & bx, const int hp)
{
	ReadId_t retval = readid2info.size();
	readid2info.push_back(ReadInfo_t(label, set, readname, seq, qv, code, strand, mate_order, bx, hp));
	return retval;
}

// addMates
////////////////////////////////////////////////////////////////

void Graph_t::addMates(ReadId_t r1, ReadId_t r2)
{
	readid2info[r1].mateid_m = r2;
	readid2info[r2].mateid_m = r1;
}

// printReads
////////////////////////////////////////////////////////////////

void Graph_t::printReads()
{
	for (unsigned int i = 0; i < readid2info.size(); ++i)
	{
		cout << i << "\t" << readid2info[i].readname_m << "\t" << readid2info[i].set_m << endl;
	}
}

// addpaired
//////////////////////////////////////////////////////////////

void Graph_t::addAlignment(const string & set,
	const string & readname,
	const string & seq,
	const string & qv,
	const int mate_id,
	char code,
	int label,
	unsigned int strand,
	const string & bx,
	const int hp)
	
{
	int rd = addRead(set, readname, seq, qv, code, label, strand, mate_id, bx, hp);
	trim(rd, seq, qv, false);
}

// addpaired
//////////////////////////////////////////////////////////////
/*
void Graph_t::addpaired(const string & set,
	const string & readname,
	const string & seq,
	const string & qv,
	const int mate_id,
	char code,
	int label,
	unsigned int strand)
{
	int rd;
	if(mate_id != 0) {
		rd = addRead(set, readname, seq, qv, code, label, strand, mate_id);
		trim(rd, seq, qv, false);
	}
	else { // no mate
		rd = addRead(set, readname, seq, qv, code, label, strand, 0);
		trim(rd, seq, qv, false);
	}	
}
*/

// build graph by trimming and loading the reads
//////////////////////////////////////////////////////////////

void Graph_t::buildgraph(Ref_t * refinfo)
{
	ref_m = refinfo;
	
	int refid = 0; 
	if(!is_ref_added) {
		string qv ((ref_m->rawseq).size(), 'K'); // create base-quality value string for reference
		refid = addRead("ref", ref_m->hdr, ref_m->rawseq, qv, 'R', REF, FWD, 0, "null", 0);
		is_ref_added = true;
		if (VERBOSE) { cerr << "refid: " << refid << endl; }
	}	
		
	for (unsigned int i = 0; i < readid2info.size(); ++i)
	{
		if ( !(readid2info[i].isjunk) ) { // skip junk (not A,C,G,T)
			string seq; // = readid2info[i].seq_m;
			string qv; // = readid2info[i].qv_m;
			int len = (readid2info[i].seq_m).length();
			int t5 = readid2info[i].trm5;
			int t3 = readid2info[i].trm3;
			unsigned int strand = readid2info[i].strand;
			//string cseq = seq;
			//string cqv = qv;
			if (t5 || t3) { 
				seq = (readid2info[i].seq_m).substr(t5, len-t5-t3); 
				qv = (readid2info[i].qv_m).substr(t5, len-t5-t3); 
			}
			else {
				seq = readid2info[i].seq_m;
				qv = readid2info[i].qv_m;
			}
																			
			if(readid2info[i].label_m == REF) {
				loadSequence(i, seq, qv, true, t5, strand);
			}
			else {
				loadSequence(i, seq, qv, false, t5, strand);
			}			
		}
	}
	
	// precompute min coverage values for each node
	MerTable_t::iterator mi;
	for (mi = nodes_m.begin(); mi != nodes_m.end(); ++mi) {
		(mi->second)->computeMinCov();
		
		sort((mi->second)->mate1_name.begin(), (mi->second)->mate1_name.end()); // sort mate1 names
		sort((mi->second)->mate2_name.begin(), (mi->second)->mate2_name.end()); // sort mate2 names
		
		(mi->second)->mate1_name.erase(unique((mi->second)->mate1_name.begin(), (mi->second)->mate1_name.end()), (mi->second)->mate1_name.end()); // remove duplicates
		(mi->second)->mate2_name.erase(unique((mi->second)->mate2_name.begin(), (mi->second)->mate2_name.end()), (mi->second)->mate2_name.end()); // remove duplicates
	}
	
	//cerr << "# of nodes: " << nodes_m.size() << endl;
	
	ref_m->computeCoverage(TMR);
	ref_m->computeCoverage(NML);
	//if (verbose) { ref_m->printKmerCoverage(NML); }
	//if (verbose) { ref_m->printKmerCoverage(TMR); }
}

// dfs to detect cycles
//////////////////////////////////////////////////////////////
bool Graph_t::hasCycle() {
	
	//cout << "Check for cycles (kmer = " << K << ")..." << endl;
	
	bool ans1 = false; 
	bool ans2 = false; 
	bool ans = false;
	
	if ( (source_m != NULL) && (sink_m != NULL) ) {
		
		MerTable_t::iterator mi;
	
		for (mi = nodes_m.begin(); mi != nodes_m.end(); ++mi) {
			Node_t * node = mi->second;

			//if (node->isRef())    { continue; }
			if (node->isSpecial())    { continue; }
			//if (node->touchRef_m) { continue; }
		
			node->setColor(WHITE);
		}
		
		hasCycleRec(source_m, F, &ans1);
		hasCycleRec(source_m, R, &ans2);
		ans = ans1 || ans2;
	}
	
	if(ans) {
		if(verbose) { cerr << "Cycle found in the graph (kmer = " << K << ")!" << endl; }
	}
	
	return ans;
}

// check if there is any node/kmer specific only to the tumor
//////////////////////////////////////////////////////////////
bool Graph_t::hasTumorOnlyKmer() {
	
	//cout << "Check for cycles (kmer = " << K << ")..." << endl;
	bool ans = false;
	
	if ( (source_m != NULL) && (sink_m != NULL) ) {
		
		MerTable_t::iterator mi;
	
		for (mi = nodes_m.begin(); mi != nodes_m.end(); ++mi) {
			Node_t * node = mi->second;

			if (node->isTumor() && !node->isNormal()) { 
				ans = true; 
				break; // exit as soon as tumor spcific node is found
			}
		}
	}
	
	return ans;
}

void Graph_t::hasCycleRec(Node_t * node, Ori_t dir, bool *ans) {
		
	if(node != NULL && !(*ans)) {
		
		node->setColor(GREY);
			
		for (unsigned int i = 0; i < node->edges_m.size(); ++i) {
						
			Edge_t & edge = node->edges_m[i];
			if (edge.isDir(dir)) {

				Node_t * other = getNode(edge);
				
				//if(other->isRef()) { continue; }
				if(other->isSpecial()) { continue; }
				
				//if (DFS_VERBOSE) { cerr << "     ==> " << other->nodeid_m << endl; }

				if(other->getColor() == GREY) { // cycle!
					//cout << "cycle found!" << endl;
					*ans = true;
					break;
				}
				if(other->getColor() == WHITE) { 
					hasCycleRec(other,  edge.destdir(), ans);
				}
			}
		}
		node->setColor(BLACK);
	}
}

// Edmonds–Karp style algorithm to enumarate the minimum number of 
// paths (source-to-sink) that cover every edge of the graph
//////////////////////////////////////////////////////////////
bool Graph_t::findRepeatsInGraphPaths(Node_t * source, Node_t * sink, Ori_t dir, Ref_t * ref)
{
	if (verbose) { cerr << endl << "looking for near-perfect repeats:" << endl; }
	if ( (source_m == NULL) || (sink_m == NULL) ) { 
		if(verbose) { cerr << "Missing source or sink" << endl; }
		return false; 
	}
	
	if(verbose) { cerr << endl << "searching from " << source->nodeid_m << " to " << sink->nodeid_m << " dir: " << dir << endl; }
		
	int complete = 0;
	bool answer = false;
	vector<Edge_t *> edges;
	
	while(true) {
		
		Path_t * path = bfs(source, sink, dir, ref);
				
		if (path == NULL) { break; }
		
		++complete;
		
		if(isAlmostRepeat(path->str(), K, MAX_MISMATCH)) {
			answer = true;
			if(verbose) { cerr << "Near-perfect repeat in assembled sequence for kmer " << K << endl; }
			break;
		}
		
		for (unsigned int i = 0; i < path->edges_m.size(); ++i) {
			(path->edges_m[i])->setFlag(1);
			edges.push_back(path->edges_m[i]);
		}

		path->reset();
		delete path;
	}
	
	// clear edge flags for next call to eka graph traversal
	for (unsigned int i = 0; i < edges.size(); ++i) {
		edges[i]->setFlag(0);
	}
	//edges.clear();
	
	return answer;
}

// printAlignment
// print alignment
void Graph_t::printAlignment(const string &ref_aln, const string &path_aln, Path_t * path) 
{
	cerr << "r':" << ref_aln << endl;
	cerr << "p':" << path_aln << " " << path->cov('A') << " [" << path->mincov('A') << " - " << path->maxcov('A') << "]" << endl;	
	cerr << "d':"; 
	for (unsigned int i = 0; i < ref_aln.length(); ++i)
	{
		if (ref_aln[i] == path_aln[i]) { ++(path->match_bp);  cerr << ' '; }
		else if (ref_aln[i] == '-')    { ++(path->ins_bp);    cerr << '^'; }
		else if (path_aln[i] == '-')   { ++(path->del_bp);    cerr << 'v'; }
		else                           { ++(path->snp_bp);    cerr << 'x'; }
	}
	cerr << "\n"; 
}

void Graph_t::printVerticalAlignment(const string &ref_aln, const string &path_aln, Path_t * path, vector<cov_t> & covN, vector<cov_t> & covT, vector<cov_t> & refcovN, vector<cov_t> & refcovT) 
{
	cerr << "Pos\tr\tp\td\tan+\tan-\tat+\tat-\trn+\trn-\trt+\trt-" << endl;
	
	int k = 0;
	int j = 0;
	for (unsigned int i = 0; i < ref_aln.length(); ++i)
	{
		cerr << k << "\t" << ref_aln[i] << "\t" << path_aln[i] << "\t";
		
		bool ins = false;  
		bool del = false;
		if (ref_aln[i] == path_aln[i]) { ++(path->match_bp);  cerr << ' '; }
		else if (ref_aln[i] == '-')    { ++(path->ins_bp);    cerr << '^'; ins = true; }
		else if (path_aln[i] == '-')   { ++(path->del_bp);    cerr << 'v'; del = true; }
		else                           { ++(path->snp_bp);    cerr << 'x'; }
		
		if (ins) {
			cerr << "\t" << covN[k].fwd << "\t" << covN[k].rev << "\t" << covT[k].fwd << "\t" << covT[k].rev;
			cerr << "\t\t\t\t";
			k++;
		}
		else if (del) {
			cerr << "\t\t\t\t";
			cerr << "\t" << refcovN[j].fwd << "\t" << refcovN[j].rev << "\t" << refcovT[j].fwd << "\t" << refcovT[j].rev;
			j++;
		}
		else {
			cerr << "\t" << covN[k].fwd << "\t" << covN[k].rev << "\t" << covT[k].fwd << "\t" << covT[k].rev;
			cerr << "\t" << refcovN[j].fwd << "\t" << refcovN[j].rev << "\t" << refcovT[j].fwd << "\t" << refcovT[j].rev;
			k++; j++;
		}
		cerr << "\n";
	} 
}

// processPath
// align path sequence to reference and parse aligment to extract mutations
//////////////////////////////////////////////////////////////
void Graph_t::processPath(Path_t * path, Ref_t * ref, FILE * fp, bool printPathsToFile,
	int &complete, int &perfect, int &withsnps, int &withindel, int &withmix) {
	
	//const string & refseq = ref->seq;
	string refseq = ref->seq;
	int HD_DISTANCE_CUTOFF = 5;
		
	path->match_bp = 0;
	path->snp_bp = 0;
	path->ins_bp = 0;
	path->del_bp = 0;
	
	if (VERBOSE) { cerr << "alignment" << endl; }
	if (VERBOSE) { cerr << "r:  " << refseq << endl; }
	if (VERBOSE) { cerr << "p:  " << path->str() << endl; }
	
	// Get alignment
	string ref_aln;
	string path_aln;
	//vector<int> cov_ref_aln;
	//vector<cov_t> cov_path_aln;
	
	vector<cov_t> coverageN = path->covDistr('N');
	vector<cov_t> coverageT = path->covDistr('T');
	
	assert(coverageN.size() == coverageT.size());
	
	string pathseq = path->str();
	
	// Run global align if strings have different length or large hamming distance	
	int hd = HammingDistance(refseq,pathseq); 
	//cerr << "HD=" << hd << endl;	
	if( (hd == -1) || (hd > HD_DISTANCE_CUTOFF) ) { 
		global_align_aff(refseq, pathseq, ref_aln, path_aln, 0, 0); 
	}
	else {
		ref_aln = refseq;
		path_aln = pathseq;
	}

	/*
	global_cov_align_aff(refseq, path->str(), coverageT_fwd, ref_aln, path_aln, cov_path_aln, 0, 0);
	coverageT = cov_path_aln;
	*/
	assert(ref_aln.length() == path_aln.length());	
	//if(verbose) { printAlignment(ref_aln, path_aln, path); }
		
	if(verbose) { 
		vector<cov_t> refcovN = (ref->getNormalCoverage());
		vector<cov_t> refcovT = (ref->getTumorCoverage());
		
		cerr << "Array size test: alt(" << coverageN.size() << ") =? ref(" << refcovN.size() << ")" << endl;
			
		printVerticalAlignment(ref_aln, path_aln, path, coverageN, coverageT, refcovN, refcovT); 
	
	}
	
	unsigned short similar_variants_count = 1;
	
	
	try {
		// scan aligned sequences for differences

		unsigned int pos_in_ref = 0; 
		unsigned int refpos  = 0; 
		unsigned int pathpos = 0;

		Node_t * spanner;

		char code = '?';
		char prev_code = '?';

		vector<Transcript_t> transcript;
		
		// cov_window keeps track of the minimum (non-0) coverage in a window of size K
		/*
		CoverageWindow_t cov_window_N;
		CoverageWindow_t cov_window_T;
	 			
		int end = min( (int)(K-1), (int)(coverageN.size()-1) );
		assert(end >= 0);
		assert(end < (int)coverageN.size());
		for (int t=0; t<end; t++) { 
			if (coverageN[t]>0) { cov_window_N.insert(coverageN[t]); }
			if (coverageT[t]>0) { cov_window_T.insert(coverageT[t]); }
		}
		*/
		
		for (unsigned int i = 0; i < ref_aln.length(); ++i) {	
						
			/*
			int toadd = min( (int)(pathpos+K-1), (int)(coverageN.size()-1));
			assert(toadd >= 0);
			assert(toadd < (int)coverageN.size());
			if (coverageN[toadd]>0) { cov_window_N.insert(coverageN[toadd]); }
			if (coverageT[toadd]>0) { cov_window_T.insert(coverageT[toadd]); }
			unsigned int old_pathpos = pathpos;
			*/
			
			prev_code = code;
			if (ref_aln[i] == '-') {
				code = '^'; // insertion
				pos_in_ref = refpos; // save value of position in reference before increment
				++pathpos;           
			}
			else if (path_aln[i] == '-') { 
				code = 'v'; // deletion
				pos_in_ref = refpos; // save value of position in reference before increment
				++refpos;            
			}
			else { 
				code = '=';
				if (ref_aln[i] != path_aln[i]) { code = 'x'; }
				pos_in_ref = refpos; // save value of position in reference before increment
				++refpos;
				++pathpos;
			}

			// pathpos is a 1-based coordinate, need to subtract 1 to correctly search 
			// into the sequence path which is indexed using 0-based coordinates
			spanner = path->pathcontig(pathpos);
			if (spanner == NULL) { cerr << "Error: path position out of range: " << pathpos << endl; break; }
			spanner->setRead2InfoList(&readid2info);
			
			bool within_tumor_node = false;
			
			if (spanner->isStatusCnt('T')) { 
			//if (spanner->isTumor() && !spanner->isNormal()) { 
				//cerr << "Within tumor only node" << endl;
				within_tumor_node = true;
				
				// print read ids
				/*
				unordered_set<ReadId_t>::const_iterator it;
				for (auto it = spanner->reads_m.begin(); it != spanner->reads_m.end(); it++) {
					cerr << readid2info[*it].readname_m.c_str() << endl;
				}
				*/
			}			

			assert(pathpos <= coverageN.size());
			assert(pathpos <= coverageT.size());
			
			int P = pathpos-1;
			/*
			if (code == 'x') { P = pathpos-1; } // snv
			else if(code == 'v') { P = pathpos-1; } // del
			else if(code == '^') { P = pathpos-1; } //ins
			*/
			cov_t COVn = coverageN[P];
			cov_t COVt = coverageT[P];
			
			cov_t REFn = ref->getCovStructAt(pos_in_ref+ref->trim5, NML);
			cov_t REFt = ref->getCovStructAt(pos_in_ref+ref->trim5, TMR);
			//cov_t REFn = ref->getCovStructAt(pos_in_ref, NML);
			//cov_t REFt = ref->getCovStructAt(pos_in_ref, TMR);
			/*
			if(old_pathpos != pathpos) { // remove old coverage only if change in path position
				assert(pathpos > 0);
				assert(pathpos <= coverageN.size());
				
				cov_window_N.remove(coverageN[pathpos-1]);
				cov_window_T.remove(coverageT[pathpos-1]);
			}
			*/
								
			if (code != '=')
			{ 
				if(verbose) {
					
					if (code == 'x') {
						cerr << (ref_aln[i] == '-' ? '*' : ref_aln[i]) << " " << (path_aln[i] == '-' ? '*' : path_aln[i]) << " " << code 
						<< " " << pos_in_ref + ref->refstart + ref->trim5 << " " << P << " " << pos_in_ref
						<< " " << spanner->nodeid_m << " " << spanner->getTotCov() << " (" <<  spanner->avgCovDistr('N') << "," << spanner->avgCovDistr('T') << ") A:("
						<< COVn.minqv_fwd << "+," << COVn.minqv_rev << "-)n,(" << COVt.minqv_fwd << "+," << COVt.minqv_rev << "-)t R:(" 
						<< REFn.fwd << "+," << REFn.rev << "-)n,(" << REFt.fwd << "+," << REFt.rev << "-)t " << spanner->reads_m.size() << " " << spanner->cntReadCode(CODE_BASTARD)							
						<< " HPalt(" << COVn.hp0 << "," << COVn.hp1 << "," << COVn.hp2 << ")n,HPalt(" << COVt.hp0 << "," << COVt.hp1 << "," << COVt.hp2 << ")t"
						<< " HPref(" << REFn.hp0 << "," << REFn.hp1 << "," << REFn.hp2 << ")n,HPref(" << REFt.hp0 << "," << REFt.hp1 << "," << REFt.hp2 << ")t"

						<< endl;						
					}
					else {
						cerr << (ref_aln[i] == '-' ? '*' : ref_aln[i]) << " " << (path_aln[i] == '-' ? '*' : path_aln[i]) << " " << code 
						<< " " << pos_in_ref + ref->refstart + ref->trim5 << " " << P << " " << pos_in_ref
						<< " " << spanner->nodeid_m << " " << spanner->getTotCov() << " (" <<  spanner->avgCovDistr('N') << "," << spanner->avgCovDistr('T') << ") A:("
						<< COVn.fwd << "+," << COVn.rev << "-)n,(" << COVt.fwd << "+," << COVt.rev << "-)t R:(" 
						<< REFn.fwd << "+," << REFn.rev << "-)n,(" << REFt.fwd << "+," << REFt.rev << "-)t " << spanner->reads_m.size() << " " << spanner->cntReadCode(CODE_BASTARD)
						<< " HPalt(" << COVn.hp0 << "," << COVn.hp1 << "," << COVn.hp2 << ")n,HPalt(" << COVt.hp0 << "," << COVt.hp1 << "," << COVt.hp2 << ")t"
						<< " HPref(" << REFn.hp0 << "," << REFn.hp1 << "," << REFn.hp2 << ")n,HPref(" << REFt.hp0 << "," << REFt.hp1 << "," << REFt.hp2 << ")t"
						<< endl;
					}
				}
				unsigned int rrpos = pos_in_ref+ref->refstart+ref->trim5;
				unsigned int ts = transcript.size();
				
				// compute previous base to the event for both reference and alternative sequences
				// [required for VCF output format]
				int pr=i-1; // referecne index
				assert(pr >= 0);
				int pa=i-1; // alternative index
				assert(pa >= 0);
				while( (ref_aln[pr] != 'A') && (ref_aln[pr] != 'C') && (ref_aln[pr]) != 'G' && (ref_aln[pr] != 'T') ) { --pr; }
				while( (path_aln[pa] != 'A') && (path_aln[pa] != 'C') && (path_aln[pa] != 'G') && (path_aln[pa] != 'T') ) { --pa; }
				
				if ( (ts > 0) && (prev_code != '=') ) // combination of multiple consecutive SNV,ins,del
				{
					if(within_tumor_node) { 
						transcript[ts-1].isSomatic = true; 
						//transcript[ts-1].nodesize = spanner->getSize();
					}
					transcript[ts-1].ref += ref_aln[i];
					transcript[ts-1].qry += path_aln[i];
					transcript[ts-1].end_pos = P; // update end position (in the path)
					transcript[ts-1].ref_end_pos = pos_in_ref; // update end position (in the ref)
					
					if ( (code == '^') && (transcript[ts-1].code == code) && (transcript[ts-1].pos == rrpos) ) { // extending insertion
						transcript[ts-1].addAltCovNml(COVn);
						transcript[ts-1].addAltCovTmr(COVt);
					} 
					else if ( (code == 'v') && (transcript[ts-1].code == code) && ((transcript[ts-1].pos + transcript[ts-1].ref.length()) == rrpos) ) { // extending deletion
						transcript[ts-1].addRefCovNml(REFn);
						transcript[ts-1].addRefCovTmr(REFt);
					} 
					else if ( (code == 'x') || (transcript[ts-1].code != code) ) { // extending complex replacement event
						transcript[ts-1].code = 'c'; // complex transcript
						transcript[ts-1].addAltCovNml(COVn);
						transcript[ts-1].addAltCovTmr(COVt);
						transcript[ts-1].addRefCovNml(REFn);
						transcript[ts-1].addRefCovTmr(REFt);
					}
				}
				
				else {
					// create new transcript for mutation
					transcript.push_back(Transcript_t(rrpos, pos_in_ref, P+1, code, 
						ref_aln[i], path_aln[i], 
						COVn, COVt, REFn, REFt,
						ref_aln[pr], path_aln[pa], 
						P, pos_in_ref, within_tumor_node));
				}
			}
		}

		if(verbose) {
			cerr << ">p_" << ref->refchr << ":" << ref->refstart << "-" << ref->refend << "_" << complete
			<< " cycle: " << path->hasCycle_m
			<< " match: " << path->match_bp
			<< " snp: "   << path->snp_bp
			<< " ins: "   << path->ins_bp
			<< " del: "   << path->del_bp;
		}
		for (unsigned int ti = 0; ti < transcript.size(); ++ti)
		{
			// if the alignment left-shifts the indel, coverage and alignment can be out of sinc. 
			// Fix: add coverage for K-1 bp after variant end position
			if (transcript[ti].code != 'x') { // only do this for indels (skip SNVs)
				if (VERBOSE) {
					cerr << endl << "coverageN.size()=" << coverageN.size() << endl; // debug
					cerr << endl << "K=" << K << endl; // debug
				}
				for (int j=0; j<=K; ++j) {
					unsigned int idx1 = transcript[ti].end_pos + j; 
					// add coverage
					if (VERBOSE) {
						cerr << "idx1=" << idx1 << endl; // debug
					}
					if (idx1 < coverageN.size()) { // check for out of range
				
						// chech if within tumor only node
						spanner = path->pathcontig(idx1);
						if (VERBOSE) {
							spanner->printx(cerr);
						}
						if (spanner == NULL) { cerr << "Error: path position out of range: " << idx1 << endl; break; }
						if (spanner->isStatusCnt('T')) { 
							//if (spanner->isTumor() && !spanner->isNormal()) 
							//cerr << "Within tumor only node" << endl;
							transcript[ti].isSomatic = true;
							//transcript[ti].nodesize = spanner->getSize();
							
							// cerr << "isSomatic" << endl;
							
							// print read ids
							//unordered_set<ReadId_t1>::const_iterator it;
							//for (auto it = spanner->reads_m.begin(); it != spanner->reads_m.end(); it++) {
							//cerr << readid2info[*it].readname_m.c_str() << endl;
							//}
						}
						// else {
						// 	spanner->print(cerr);
						// }
						// debug start
						if (VERBOSE) {
							cerr << 'N' << endl;
							spanner->isStatusCnt('N');
							cerr << 'B' << endl;
							spanner->isStatusCnt('B');
							cerr << 'E' << endl;
							spanner->isStatusCnt('E');
							cov_t c_n = coverageN[idx1];
							cov_t c_t = coverageT[idx1];
							cerr << "c_n.fwd=" << c_n.fwd << endl;
							cerr << "c_n.rev=" << c_n.rev << endl;
							cerr << "c_t.fwd=" << c_t.fwd << endl;
							cerr << "c_t.rev=" << c_t.rev << endl;
						}
						// debug end
						transcript[ti].addAltCovNml(coverageN[idx1]);
						transcript[ti].addAltCovTmr(coverageT[idx1]);
					}
					unsigned int idx2 = transcript[ti].ref_end_pos + ref->trim5 + j;
					// debug start
					if (VERBOSE) {
						cov_t r_n = ref->getCovStructAt(idx2, NML);
						cov_t r_t = ref->getCovStructAt(idx2, TMR);
						cerr << "r_n.fwd=" << r_n.fwd << endl;
						cerr << "r_n.rev=" << r_n.rev << endl;
						cerr << "r_t.fwd=" << r_t.fwd << endl;
						cerr << "r_t.rev=" << r_t.rev << endl;
					}
					// debug end
					transcript[ti].addRefCovNml(ref->getCovStructAt(idx2, NML));
					transcript[ti].addRefCovTmr(ref->getCovStructAt(idx2, TMR));
				}
			}
			
			// updated transcript stats
			transcript[ti].updateStats();
			
			// select coverage according to mutation type and sample
			unsigned short RCNF = transcript[ti].getMinRefCovNfwd(); // ref cov normal fwd
			unsigned short RCNR = transcript[ti].getMinRefCovNrev(); // ref cov normal rev

			unsigned short RCTF = transcript[ti].getMinRefCovTfwd(); // ref cov tumor fwd
			unsigned short RCTR = transcript[ti].getMinRefCovTrev(); // ref cov tumor rev

			unsigned short ACNF = transcript[ti].getMinCovNfwd(); // alt normal cov fwd
			unsigned short ACNR = transcript[ti].getMinCovNrev(); // alt normal cov rev			
			
			
			// NOT skip zero 
			/*
			if (transcript[ti].code != 'x') { // for indels skip over zero coverage values (due to the +k coverage values added after variant end postion
				ACNF = transcript[ti].getMinNon0CovNfwd(); // alt normal cov fwd
				ACNR = transcript[ti].getMinNon0CovNrev(); // alt normal cov rev
			}
			*/
			
			unsigned short ACTF = transcript[ti].getMinCovTfwd(); // alt tumor cov fwd
			unsigned short ACTR = transcript[ti].getMinCovTrev(); // alt tumor cov rev
			
			unsigned short HP0RN = transcript[ti].getMinRefCovNhp0();
			unsigned short HP1RN = transcript[ti].getMinRefCovNhp1();
			unsigned short HP2RN = transcript[ti].getMinRefCovNhp2();
			
			unsigned short HP0RT = transcript[ti].getMinRefCovThp0();
			unsigned short HP1RT = transcript[ti].getMinRefCovThp1();
			unsigned short HP2RT = transcript[ti].getMinRefCovThp2();
			
			unsigned short HP0AN = transcript[ti].getMinCovNhp0();
			unsigned short HP1AN = transcript[ti].getMinCovNhp1();
			unsigned short HP2AN = transcript[ti].getMinCovNhp2();
			
			unsigned short HP0AT = transcript[ti].getMinCovThp0();
			unsigned short HP1AT = transcript[ti].getMinCovThp1();
			unsigned short HP2AT = transcript[ti].getMinCovThp2();
			
			if(transcript[ti].isSomatic) {
				RCNF = transcript[ti].getAvgRefCovNfwd(); // ref cov normal fwd
				RCNR = transcript[ti].getAvgRefCovNrev(); // ref cov normal rev

				RCTF = transcript[ti].getAvgRefCovTfwd(); // ref cov tumor
				RCTR = transcript[ti].getAvgRefCovTrev(); // ref cov tumor
				//ACNF = transcript[ti].getMinCovNfwd(); // alt normal cov fwd
				//ACNR = transcript[ti].getMinCovNrev(); // alt normal cov rev
				ACNF = 0; // alt normal cov fwd
				ACNR = 0; // alt normal cov rev
				//ACTF = transcript[ti].getMedianCovTfwd(); // alt tumor cov fwd
				//ACTR = transcript[ti].getMedianCovTrev(); // alt tumor cov rev
								
				HP0RT = transcript[ti].getAvgRefCovThp0();
				HP1RT = transcript[ti].getAvgRefCovThp1();
				HP2RT = transcript[ti].getAvgRefCovThp2();
			
				HP0RN = transcript[ti].getAvgRefCovNhp0();
				HP1RN = transcript[ti].getAvgRefCovNhp1();
				HP2RN = transcript[ti].getAvgRefCovNhp2();
					
				HP0AN = 0; HP1AN = 0; HP2AN = 0;
			}
			
			//int ACTF = (transcript[ti].code=='x') ? transcript[ti].getMinCovTfwd() : transcript[ti].getMedianCovTfwd(); // alt tumor cov fwd
			//int ACTR = (transcript[ti].code=='x') ? transcript[ti].getMinCovTrev() : transcript[ti].getMedianCovTrev(); // alt tumor cov rev
					
			if(verbose) { cerr << " " << transcript[ti].isSomatic << " " << transcript[ti].pos << ":" << transcript[ti].ref << "|" << transcript[ti].qry << "|R:(" << 
						RCNF << "+," << RCNR << "-)n,(" << 
						RCTF << "+," << RCTR << "-)t|A:(" << 								
						ACNF << "+," << ACNR << "-)n,(" << 
						ACTF << "+," << ACTR << "-)t|HPref(" <<
							
						HP0RN << "," << HP1RN << "," << HP2RN << ")n,(" << 
						HP0RT << "," << HP1RT << "," << HP2RT << ")t|HPalt(" << 
						HP0AN << "," << HP1AN << "," << HP2AN << ")n,(" << 
						HP0AT << "," << HP1AT << "," << HP2AT << ")t|" << 
						
						transcript[ti].prev_bp_ref << "|" << transcript[ti].prev_bp_alt; 
			}
			
			// save variant into DB
			if( (ACNF > 0) || (ACNR > 0) || (ACTF > 0) || (ACTR > 0) ) {
				
				// annotate if in STR
				int LEN;
				string MOTIF = "";
				stringstream STR;
					
				bool ans = findTandems(pathseq, "test", MAX_UNIT_LEN, MIN_REPORT_UNITS, MIN_REPORT_LEN, DIST_FROM_STR, transcript[ti].start_pos, LEN, MOTIF);
				if(ans) { 
					STR << LEN << MOTIF;
					//cerr << "STR = " << STR.str() << endl; 					
				}
				
				pair <unsigned short,unsigned short> RCN (RCNF,RCNR);   // reference coverage normal
				pair <unsigned short,unsigned short> RCT (RCTF,RCTR);   // reference coverage tumor        
				pair <unsigned short,unsigned short> ACN (ACNF,ACNR);   // alternative coverage normal        
				pair <unsigned short,unsigned short> ACT (ACTF,ACTR);   // alternative coverage tumor
				
				array<unsigned short,3> HPRN = {HP1RN, HP2RN, HP0RN}; // reference HP normal
				array<unsigned short,3> HPRT = {HP1RT, HP2RT, HP0RT}; // reference HP tumor
				array<unsigned short,3> HPAN = {HP1AN, HP2AN, HP0AN}; // alternative HP normal
				array<unsigned short,3> HPAT = {HP1AT, HP2AT, HP0AT}; // alternative HP tumor
				
				string bxset_ref_N = "";
				string bxset_ref_T = "";
				string bxset_alt_N = "";
				string bxset_alt_T = "";
					
				if (LR_MODE) {
					bxset_ref_N = ref->getBXsetAt(transcript[ti].ref_pos-1, transcript[ti].ref_end_pos-1, refseq, NML);
					bxset_ref_T = ref->getBXsetAt(transcript[ti].ref_pos-1, transcript[ti].ref_end_pos-1, refseq, TMR);			
					bxset_alt_N = getBXsetAt(transcript[ti].start_pos-2, transcript[ti].end_pos-1, pathseq, NML);
					bxset_alt_T = getBXsetAt(transcript[ti].start_pos-2, transcript[ti].end_pos-1, pathseq, TMR);
				}
				
				cerr << "==================s" << endl;
 				cerr << LR_MODE << " " << ref->refchr << " " << transcript[ti].pos-1 << " " << transcript[ti].ref << " " << transcript[ti].qry << endl;
				cerr << RCN.first << " " << RCN.second << " " << RCT.first << " " << RCT.second << " " << ACN.first << " " << ACN.second << " " << ACT.first << " " << ACT.second << endl;
				cerr << "==================e" << endl;

				vDB->addVar(Variant_t(LR_MODE, ref->refchr, transcript[ti].pos-1, transcript[ti].ref, transcript[ti].qry, transcript[ti].isSomatic,
					similar_variants_count,
					RCN, RCT, ACN, ACT,
					HPRN, HPRT, HPAN, HPAT,
					transcript[ti].prev_bp_ref, transcript[ti].prev_bp_alt, K, STR.str(), transcript[ti].code,
					bxset_ref_N, bxset_ref_T, bxset_alt_N, bxset_alt_T));
			}
		}
		if(verbose) { cerr << endl; }

		if      ((path->snp_bp + path->ins_bp + path->del_bp) == 0) { ++perfect;   }
		else if ((path->snp_bp) == 0)                               { ++withindel; }
		else if ((path->ins_bp + path->del_bp) == 0)                { ++withsnps;  }
		else                                                        { ++withmix;   }

		if(printPathsToFile) {
			fprintf(fp,  ">p_%s:%d-%d_%d len=%d cov=%0.2f mincov=%0.2f maxcov=%0.2f pathlen=%d hasCycle=%d match=%d snp=%d ins=%d del=%d pathstr=%s\n",
				ref->refchr.c_str(), ref->refstart, ref->refend, complete, 
				path->strlen(), path->cov('A'), path->mincov('A'), path->maxcov('A'), path->pathlen(), 
				path->hasCycle_m, path->match_bp, path->snp_bp, path->ins_bp, path->del_bp, path->pathstr().c_str());
			
			fprintf(fp, "%s\n", path->str().c_str());
		}

		for (unsigned int i = 0; i < path->nodes_m.size(); ++i)
		{
			Node_t * cur = path->nodes_m[i];
			++(cur->onRefPath_m);
		}
	}
	catch(std::out_of_range& e) {
   		cerr << "An exception occurred: " << e.what( ) << endl;
 	}
	catch (...) { 
		cout << "default exception" << endl; 
	}
	
}

// processShortPath
// compute poartial alignment to reference and parse aligment to extract mutations
//////////////////////////////////////////////////////////////
void Graph_t::processShortPath(Node_t * source, Ref_t * ref, FILE * fp, bool printPathsToFile,
	int &complete, int &perfect, int &withsnps, int &withindel, int &withmix) {

	const string & refseq = ref->seq;
	
	Node_t * node = getNode(source->edges_m[0]);
	++(node->onRefPath_m);

	string str = node->str_m;

	if (source->edges_m[0].destdir() == R)
	{
		str = CanonicalMer_t::rc(str);
	}


	VERBOSE = 1;

	if (VERBOSE) { cerr << "partial align" << endl; }

	if (VERBOSE) { cerr << "alignment" << endl; }
	if (VERBOSE) { cerr << "r:  " << refseq << endl; }
	if (VERBOSE) { cerr << "p:  " << str << endl; }

	// Get alignment
	string ref_aln;
	string path_aln;
	string cov_aln;

	global_align_aff(refseq, str, ref_aln, path_aln, 1, 0);

	if (VERBOSE) 
	{ 
		cerr << "r': " << ref_aln << endl;  
		cerr << "p': " << path_aln << endl; 
	}

	int match_bp = 0;
	int snp_bp = 0;
	int ins_bp = 0;
	int del_bp = 0;

	assert(ref_aln.length() == path_aln.length());

	for (unsigned int i = 0; i < ref_aln.length(); ++i)
	{
		if (ref_aln[i] == path_aln[i]) { ++match_bp; }
		else if (ref_aln[i] == '-')    { ++ins_bp; }
		else if (path_aln[i] == '-')   { ++del_bp; }
		else                           { ++snp_bp; }
	}

	cerr << ">sp" << complete
		<< " match: " << match_bp
		<< " snp: " << snp_bp
		<< " ins: " << ins_bp
		<< " del: " << del_bp
		<< endl;

	if      ((snp_bp + ins_bp + del_bp) == 0) { ++perfect;   }
	else if ((snp_bp) == 0)                   { ++withindel; }
	else if ((ins_bp + del_bp) == 0)          { ++withsnps;  }
	else                                      { ++withmix;   }

	if(printPathsToFile) {
		fprintf(fp,  ">p_%d len=%d cov=%0.2f mincov=%0.2f maxcov=%0.2f pathlen=%d match=%d snp=%d ins=%d del=%d pathstr=%s\n",
			complete, (int) str.length(), node->getTotCov(), node->getTotCov(), node->getTotCov(), 1, match_bp, snp_bp, ins_bp, del_bp, "shortmatch");
		
		fprintf(fp, "%s\n", str.c_str());
	}
}

// bfs
//////////////////////////////////////////////////////////////
Path_t * Graph_t::bfs(Node_t * source, Node_t * sink, Ori_t dir, Ref_t * ref)
{
	const string & refseq = ref->seq;
	
	int complete = 0;
	int toolong = 0;
	int deadend = 0;
	int shortpaths = 0;
	int allcycles = 0;
	int visit = 0;
	
	int reflen = refseq.length();

	deque<Path_t *> Q;
	
	Path_t * path = new Path_t(K);
	path->nodes_m.push_back(source);
	path->dir_m = dir;
	path->len_m = K;
	//Path_t * best = new Path_t(path,K);
	Path_t * best = NULL;

	Q.push_back(path);
	
	while (!Q.empty())
	{
		++visit;

		if ((DFS_LIMIT) && (visit > DFS_LIMIT)) {
			if(verbose) { cerr << "WARNING: DFS_LIMIT (" << DFS_LIMIT << ") exceeded" << endl; }
			break;
		}

		path = Q.front();
		Q.pop_front();

		Node_t * cur = path->curNode();

		// if sink is found and at lest one of the edges in the path was not discovered (flag == 0)
		if ( (cur == sink) && (path->flag == 0) )
		{
			// success!
			++complete;
			if (best == NULL) { best = new Path_t(path,K); }
			else if(path->score > best->score) { 
				Path_t * old_best = best;
				best = new Path_t(path,K); // keep best path
				delete old_best;
			} 
			//break;
		}
		else if (path->len_m > reflen + MAX_INDEL_LEN)
		{
			// abort
			++toolong;
			//cerr << "too long: " << path->pathstr() << " " << path->str() << endl;
		}
		else
		{
			int tried = 0;

			for (unsigned int i = 0; i < cur->edges_m.size(); ++i)
			{
				Edge_t * edge = &(cur->edges_m[i]);	
				
				if (edge->isDir(path->dir_m))
				{
					++tried;

					Node_t * other = getNode(*edge);

					if (!path->hasCycle_m && path->hasCycle(other))
					{
						++allcycles;
						//cerr << "Cycle detected in BFS!!" << endl;
					}

					Path_t * newpath = new Path_t(path,K);

					newpath->nodes_m.push_back(other);
					newpath->edges_m.push_back(edge);
					newpath->edgedir_m.push_back(edge->dir_m);
					newpath->dir_m = edge->destdir();
					newpath->len_m = path->len_m + other->strlen() - K + 1;
					newpath->flag = path->flag * edge->getFlag(); // update flag;
					if(edge->getFlag() == 0) { newpath->score = path->score + 1; }// update score;
					
					Q.push_back(newpath);
				}
			}

			if (tried == 0)
			{
				++deadend;
				//cerr << "deadend: " <<  cur->nodeid_m << endl;
			}
		}

		path->reset();
		delete path;
	}

	while (!Q.empty())
	{
		Path_t * path = Q.front();
		delete path;
		Q.pop_front();
	}

	if (complete == 0)
	{
		// didn't find an end-to-end path

		if (visit == 2)
		{
			// source to single node
			assert(source->edges_m.size() == 1);

			++complete;
			++shortpaths;			
		}
		//path = NULL;
		best = NULL;
	}
	
	return best;
}

// Edmonds–Karp style algorithm to enumarate the minimum number of 
// paths (source-to-sink) that cover every edge of the graph
//////////////////////////////////////////////////////////////
void Graph_t::eka(Node_t * source, Node_t * sink, Ori_t dir, 
	Ref_t * ref, FILE * fp, bool printPathsToFile)
{
	if(verbose) { cerr << endl << "searching from " << source->nodeid_m << " to " << sink->nodeid_m << " dir: " << dir << endl; }
	
	int complete = 0;
	//int toolong = 0;
	//int deadend = 0;
	//int visit = 0;
	//int shortpaths = 0;
	int allcycles = 0;

	int perfect   = 0;
	int withsnps  = 0;
	int withindel = 0;
	int withmix   = 0;
	
	while(true) {
		
		Path_t * path = bfs(source, sink, dir, ref);
				
		if (path == NULL) { break; }
				
		if (path->hasCycle_m) { ++allcycles; }
		++complete;
				
		//if(path->hasTumorOnlyNode()) {
			processPath(path, ref, fp, printPathsToFile, complete, perfect, withsnps, withindel, withmix);
		//}
			
		for (unsigned int i = 0; i < path->edges_m.size(); ++i) {
			(path->edges_m[i])->setFlag(1);
		}

		path->reset();
		delete path;
	}
	
	int withmixindel = withmix + withindel;
	//int withmixsnp   = withmix + withsnps;
	//int withvar      = withsnps + withindel + withmix;

	//cerr << " refcomp: "    << ref_m->refcomp
	//	<< " refnodes: "   << ref_m->refnodes-2
	//	<< " visit: "      << visit
	//	<< " complete: "   << complete
	//	<< " allcycles: "  << allcycles
	//	<< " shortpaths: " << shortpaths
	//	<< " toolong: "    << toolong
	//	<< " deadend: "    << deadend << endl;
	
	if(verbose) {
		cerr << " refcomp: "   << ref_m->refcomp
		<< " refnodes: "   << ref_m->refnodes-2
		<< " complete: "   << complete 
		<< " allcycles: "  << allcycles << endl;

		cerr << " perfect: "     << perfect
		<< " withsnps: "     << withsnps
		<< " withindel: "    << withindel
		<< " withmix: "      << withmix 
		<< " withmixindel: " << withmixindel
		<< endl;
	}

	//if(printPathsToFile) {
	//	fprintf(fp, ">stats\treflen=%d\tnumreads=%d\tcov=%0.02f\ttrim5=%d\ttrim3=%d\tnodes=%d\trefnodes=%d\tcomp=%d\trefcomp=%d\tvisit=%d\tcomplete=%d\tallcycles=%d\tshortpath=%d\ttoolong=%d\tdeadend=%d\tperfect=%d\twithsnps=%d\twithindel=%d\twithmix=%d\twithmixindel=%d\twithmixsnp=%d\twithvar=%d\n",
	//		(int) ref_m->seq.length(), (int) readid2info.size(), ((float)totalreadbp_m / (float) ref_m->seq.length()),
	//		ref_m->trim5, ref_m->trim3, (int) nodes_m.size()-2, ref_m->refnodes-2, ref_m->allcomp, ref_m->refcomp,
	//		visit, complete, allcycles, shortpaths, toolong, deadend, perfect, withsnps, withindel, withmix, withmixindel, withmixsnp, withvar);
	//}
}

// dfs
//////////////////////////////////////////////////////////////

void Graph_t::dfs(Node_t * source, Node_t * sink, Ori_t dir, 
	Ref_t * ref, FILE * fp, bool printPathsToFile)
{
	cerr << endl << "searching from " << source->nodeid_m << " to " << sink->nodeid_m << " dir: " << dir << endl;

	const string & refseq = ref->seq;

	int complete = 0;
	int toolong = 0;
	int deadend = 0;
	int visit = 0;
	int shortpaths = 0;
	int allcycles = 0;

	int perfect   = 0;
	int withsnps  = 0;
	int withindel = 0;
	int withmix   = 0;

	int reflen = refseq.length();

	deque<Path_t *> Q;

	Path_t * start = new Path_t(K);
	start->nodes_m.push_back(source);
	start->dir_m = dir;
	start->len_m = K;

	Q.push_back(start);

	bool DFS_VERBOSE = 0;
	bool OLD_VERBOSE = VERBOSE;
	
	while (!Q.empty())
	{
		++visit;

		if ((DFS_LIMIT) && (visit > DFS_LIMIT))
		{
			if(verbose) { cerr << "WARNING: DFS_LIMIT (" << DFS_LIMIT << ") exceeded" << endl; }
			break;
		}

		Path_t * path = Q.front();
		Q.pop_front();

		Node_t * cur = path->curNode();

		if (DFS_VERBOSE) { cerr << " --> " << cur->nodeid_m << " [" << Q.size() << "]" << endl; }

		if (cur == sink)
		{
			// success!
			++complete;
			processPath(path, ref, fp, printPathsToFile, complete, perfect, withsnps, withindel, withmix);
		}
		else if (path->len_m > reflen + MAX_INDEL_LEN)
		{
			// abort
			++toolong;
			//cerr << "too long: " << path->pathstr() << " " << path->str() << endl;
		}
		else
		{
			int tried = 0;

			for (unsigned int i = 0; i < cur->edges_m.size(); ++i)
			{
				Edge_t & edge = cur->edges_m[i];
				if (edge.isDir(path->dir_m))
				{
					++tried;

					Node_t * other = getNode(edge);
					if (DFS_VERBOSE) { cerr << "     ==> " << other->nodeid_m << endl; }

					if (!path->hasCycle_m && path->hasCycle(other))
					{
						++allcycles;
					}

					Path_t * newpath = new Path_t(path,K);

					newpath->nodes_m.push_back(other);
					newpath->edgedir_m.push_back(edge.dir_m);
					newpath->dir_m = edge.destdir();
					newpath->len_m = path->len_m + other->strlen() - K + 1;

					Q.push_back(newpath);
				}
			}

			if (tried == 0)
			{
				++deadend;
				//cerr << "deadend: " <<  cur->nodeid_m << endl;
			}
		}

		delete path;
	}

	while (!Q.empty())
	{
		Path_t * path = Q.front();
		delete path;
		Q.pop_front();
	}

	if (complete == 0)
	{
		// didn't find an end-to-end path

		if (visit == 2)
		{
			// source to single node
			assert(source->edges_m.size() == 1);

			++complete;
			++shortpaths;
			
			processShortPath(source, ref, fp, printPathsToFile, complete, perfect, withsnps, withindel, withmix);
		}
	}

	VERBOSE = OLD_VERBOSE;

	int withmixindel = withmix + withindel;
	int withmixsnp   = withmix + withsnps;
	int withvar      = withsnps + withindel + withmix;

	cerr << " refcomp: "    << ref_m->refcomp
		<< " refnodes: "   << ref_m->refnodes-2
		<< " visit: "      << visit
		<< " complete: "   << complete
		<< " allcycles: "  << allcycles
		<< " shortpaths: " << shortpaths
		<< " toolong: "    << toolong
		<< " deadend: "    << deadend << endl;

	cerr << " perfect: "      << perfect
		<< " withsnps: "     << withsnps
		<< " withindel: "    << withindel
		<< " withmix: "      << withmix 
		<< " withmixindel: " << withmixindel
		<< endl;

	if(printPathsToFile) {
		fprintf(fp, ">stats\treflen=%d\tnumreads=%d\tcov=%0.02f\ttrim5=%d\ttrim3=%d\tnodes=%d\trefnodes=%d\tcomp=%d\trefcomp=%d\tvisit=%d\tcomplete=%d\tallcycles=%d\tshortpath=%d\ttoolong=%d\tdeadend=%d\tperfect=%d\twithsnps=%d\twithindel=%d\twithmix=%d\twithmixindel=%d\twithmixsnp=%d\twithvar=%d\n",
			(int) ref_m->seq.length(), (int) readid2info.size(), ((float)totalreadbp_m / (float) ref_m->seq.length()),
			ref_m->trim5, ref_m->trim3, (int) nodes_m.size()-2, ref_m->refnodes-2, ref_m->allcomp, ref_m->refcomp,
			visit, complete, allcycles, shortpaths, toolong, deadend, perfect, withsnps, withindel, withmix, withmixindel, withmixsnp, withvar);
	}
}

string Graph_t::nodeColor(Node_t * cur, string & who)
{
	string color = COLOR_ALL;

	//if      (cur->nodeid_m == "source") { color = COLOR_SOURCE; }
	//else if (cur->nodeid_m == "sink")   { color = COLOR_SINK; }
	
	if      (cur->isSource()) { color = COLOR_SOURCE; }
	else if (cur->isSink())   { color = COLOR_SINK; }
	else if (cur->touchRef_m) { color = COLOR_TOUCH; }

	stringstream whostr;

	map<string, int> whocnt;

	unordered_set<ReadId_t>::const_iterator si;
	for (si = cur->reads_m.begin(); si != cur->reads_m.end(); ++si)
	{
		++(whocnt[readid2info[*si].set_m]);
	}

	bool isTumor = cur->isTumor();
	bool isNormal = cur->isNormal();
	
	map<string, int>::iterator mi;
	for (mi = whocnt.begin(); mi != whocnt.end(); ++mi)
	{
		if (mi != whocnt.begin()) { whostr << " "; }
		whostr << mi->first << ":" << mi->second;
		//cout << whostr.str() << endl;
	}

	double avgcov = ((double) totalreadbp_m) / ((double) ref_m->rawseq.length());

	if ((color == COLOR_ALL) && ((cur->getMinCov() <= LOW_COV_THRESHOLD) || (cur->getMinCov() <= (MIN_COV_RATIO*avgcov))) ) 
	{
		color = COLOR_LOW;
	}
	else if (isTumor && !isNormal)
	{
		color = COLOR_TUMOR;
	}
	else if (!isTumor && isNormal)
	{
		color = COLOR_NORMAL;
	}
	else if (isTumor && isNormal)
	{
		color = COLOR_SHARED;
	}
	
	who = whostr.str();

	return color;
}

string Graph_t::edgeColor(Node_t * cur, Edge_t & e)
{
	string color = COLOR_ALL;

	string w;
	string c1 = nodeColor(cur, w);
	string c2 = nodeColor(getNode(e), w);

	if ((c1 == COLOR_LOW) || (c2 == COLOR_LOW))
	{
		color = COLOR_LOW;
	}
	else if ((c1 == COLOR_NOVO) && (c2 == COLOR_NOVO))
	{
		color = COLOR_NOVO;
	}

	return color;
}



// printDot
//////////////////////////////////////////////////////////////
void Graph_t::printDot(const string & filename, int compid)
{
	if(verbose) { cerr << "Saving graph: " << filename << endl; }

	FILE * fp = xfopen(filename, "w");

	if (PRINT_DOT_READS)
	{
		for (unsigned int i = 0; i < readid2info.size(); ++i)
		{
			fprintf(fp, "// %s %d %s -> %d (%s)\n",
				readid2info[i].set_m.c_str(),
				i, 
				readid2info[i].readname_m.c_str(),
				readid2info[i].mateid_m, 
				readid2info[i].contigid_m.c_str());
		}

		fprintf(fp, "\n\n");
	}
	
	fprintf(fp, "digraph structs{\n");

    int printstrings = 0;

    if (printstrings)
    {
	  fprintf(fp, "  graph [bgcolor=white,size=\"8.5,11\",ratio=fill,center=true]\n");
	  fprintf(fp, "  node [shape=record];\n");
	  fprintf(fp, "  rankdir=LR\n");
    }
    else
    {
	  fprintf(fp, "  graph [bgcolor=black,size=\"80,110\",layout=neato,ratio=compress,rankdir=LR]\n");
	  fprintf(fp, "  node [shape=circle,style=filled,fontsize=1,fixedsize=true,hight=1,width=1];\n");
	  fprintf(fp, "  edge [fixedsize=true,len=1.2];\n");
    }

	int nodes = 0;

	MerTable_t::iterator mi;
	for (mi = nodes_m.begin(); mi != nodes_m.end(); ++mi)
	{
		// only print nodes in the selected component
		if(mi->second->component_m != compid) { continue; }
		//if( (mi->second->component_m != compid) && !(mi->second->isSource()) && !(mi->second->isSink()) ) { continue; }
		 
		++nodes;

		Node_t * cur = mi->second;
		cur->setRead2InfoList(&readid2info);

		string who;
		string color = nodeColor(cur, who);
		
		string shape = "circle";
		if ( (cur->isSource()) || (cur->isSink()) ) { shape = "diamond"; }

		if (NODE_STRLEN == 0)
		{
			fprintf(fp, "  %s [label=\". | <F> . | <R> .\" color=\"%s\"]\n",
				cur->nodeid_m.c_str(),
				color.c_str());
		}
		else if (cur->strlen() > NODE_STRLEN)
		{
			fprintf(fp, "  %s [label=\"%d:%s | <F> %s... | <R> len=%d cov=%0.02f rd:%d B:%d %s\" color=\"%s\" shape=\"%s\"]\n",
				cur->nodeid_m.c_str(),
				nodes,
				cur->nodeid_m.c_str(),
				cur->str_m.substr(0, NODE_STRLEN).c_str(),
				cur->strlen(),
				cur->getTotCov(),
				(int) cur->reads_m.size(),
				cur->cntReadCode(CODE_BASTARD),
				who.c_str(),
				color.c_str(),
				shape.c_str());
		}
		else
		{
			fprintf(fp, "  %s [label=\"%d:%s | <F> %s | <R> len=%d cov=%0.02f rd:%d B:%d %s\" color=\"%s\" shape=\"%s\"]\n",
				cur->nodeid_m.c_str(),
				nodes,
				cur->nodeid_m.c_str(),
				cur->str_m.substr(0, NODE_STRLEN).c_str(),
				cur->strlen(),
				cur->getTotCov(),
				(int) cur->reads_m.size(),
				cur->cntReadCode(CODE_BASTARD),
				who.c_str(),
				color.c_str(),
				shape.c_str());
		}

		if (PRINT_DOT_READS)
		{
			fprintf(fp, "  //reads:");

			unordered_set<ReadId_t>::const_iterator ri;

			for (ri = cur->reads_m.begin(); ri != cur->reads_m.end(); ++ri)
			{
				//fprintf(fp, " %s", readid2info[e.readids_m[j]].readname_m.c_str());
				fprintf(fp, " %d", *ri);
			}

			fprintf(fp, "\n");

			fprintf(fp, "  //readstarts:");

			for (unsigned int i = 0; i < cur->readstarts_m.size(); ++i)
			{
				fprintf(fp, " %d:%c%d",
					cur->readstarts_m[i].readid_m,
					(cur->readstarts_m[i].ori_m == F) ? '+' : '-',
					cur->readstarts_m[i].nodeoffset_m);
			}

			fprintf(fp, "\n");

			fprintf(fp, "  //links:");

			ContigLinkMap_t::iterator li;
			for (li = cur->contiglinks_m.begin();
			li != cur->contiglinks_m.end();
			++li)
			{
				fprintf(fp, " %s(%d)", li->first.c_str(), li->second->linkCnt());
			}

			fprintf(fp, "\n");
		}

		for (unsigned int i = 0; i < cur->edges_m.size(); ++i)
		{
			Edge_t & e = cur->edges_m[i];

			bool printed = 0;

			if (e.dir_m != RR)
			{
				if ((e.dir_m == FF) || (cur->nodeid_m <= e.nodeid_m))
				{
					string ecolor = edgeColor(cur, e);

					printed = 1;
					fprintf(fp, "    %s:%c -> %s:%c [arrowhead=\"normal\" arrowtail=\"normal\" color=\"%s\"]\n",
						cur->nodeid_m.c_str(), Edge_t::toString(e.startdir()),
						e.nodeid_m.c_str(),    Edge_t::toString(e.destdir()),
						ecolor.c_str());
				}
			}

			if (!printed)
			{
				fprintf(fp, "    //%s:%c -> %s:%c\n",
					cur->nodeid_m.c_str(), Edge_t::toString(e.startdir()),
					e.nodeid_m.c_str(),    Edge_t::toString(e.destdir()));
			}


			if (PRINT_DOT_READS)
			{
				fprintf(fp, "    //reads:");

				for (unsigned int j = 0; j < e.readids_m.size(); ++j)
				{
					//fprintf(fp, " %s", readid2info[e.readids_m[j]].readname_m.c_str());
					fprintf(fp, " %d", e.readids_m[j]);
				}

				fprintf(fp, "\n");
			}
		}

		fprintf(fp, "\n\n");
	}

	fprintf(fp, "}\n");

	xfclose(fp);
}


// printFasta: print all contigs
//////////////////////////////////////////////////////////////

void Graph_t::printFasta(const string & filename)
{
	cerr << "Saving fasta: " << filename << endl;

	FILE * fp = xfopen(filename, "w");

	int nodes = 0;

	MerTable_t::iterator mi;
	for (mi = nodes_m.begin(); mi != nodes_m.end(); ++mi)
	{
		++nodes;
		Node_t * cur = mi->second;

        int fdeg = 0;
        int rdeg = 0;

        for (unsigned int i = 0; i < cur->edges_m.size(); ++i)
        {
            Edge_t & edge = cur->edges_m[i];
            if (edge.isDir(F)) { ++fdeg; }
            else               { ++rdeg; }
        }

		fprintf(fp, ">%d:%s len=%d cov=%0.2f fdeg=%d rdeg=%d\n", nodes, cur->nodeid_m.c_str(), cur->strlen(), cur->getTotCov(), fdeg, rdeg);
		fprintf(fp, "%s\n",  cur->str_m.c_str());
	}

	xfclose(fp);
}

// printPairs: print all pairs of neighboring contigs
//////////////////////////////////////////////////////////////

void Graph_t::printPairs(const string & filename)
{
	cerr << "Saving pairs fasta: " << filename << endl;

	FILE * fp = xfopen(filename, "w");

	int nodes = 0;

	Path_t pairpath(K);

	MerTable_t::iterator mi;
	for (mi = nodes_m.begin(); mi != nodes_m.end(); ++mi)
	{
		Node_t * cur = mi->second;

		if (cur->edges_m.size() == 0)
		{
		// print isolated contigs

			++nodes;
			fprintf(fp, ">%d:%s len=%d cov=%0.2f\n", nodes, cur->nodeid_m.c_str(), cur->strlen(), cur->getTotCov());
			fprintf(fp, "%s\n",  cur->str_m.c_str());
		}
		else
		{
			// print contig pairs, making sure a given pair is only printed once

			for (unsigned int i = 0; i < cur->edges_m.size(); ++i)
			{
				Edge_t & e = cur->edges_m[i];

				if (e.dir_m != RR) 
				{
					if ((e.dir_m == FF) || (cur->nodeid_m < e.nodeid_m))
					{
						++nodes;

						Node_t * other = getNode(e);

						pairpath.reset();
						pairpath.nodes_m.push_back(cur);
						pairpath.edgedir_m.push_back(e.dir_m);
						pairpath.nodes_m.push_back(other);

						string pathstr = pairpath.pathstr();
						string path    = pairpath.str();
						float  cov     = pairpath.cov('T');

						int    strlen  = path.length();

						fprintf(fp, ">%d:%s len=%d cov=%0.2f\n", nodes, pathstr.c_str(), strlen, cov);
						fprintf(fp, "%s\n",  path.c_str());
					}
				}
			}
		}
	}

	xfclose(fp);
}


// mark source and sink nodes
//////////////////////////////////////////////////////////////

void Graph_t::markRefEnds(Ref_t * refinfo, int compid)
{
	if (VERBOSE) 
	{ 
		cerr << "Looking at " << refinfo->hdr << " " << refinfo->seq.length() 
			<< " " << refinfo->seq << endl;
	}

	//ref_m = refinfo;
	int refid = 0; 
	if(!is_ref_added) {
		string qv ((ref_m->rawseq).size(), 'K'); // create base-quality value string for reference
		refid = addRead("ref", ref_m->hdr, ref_m->rawseq, qv, 'R', REF, FWD, 0, "null", 0);
		is_ref_added = true;
		if (VERBOSE) { cerr << "refid: " << refid << endl; }
	}

	//loadSequence(refid, ref_m->seq, 1, 0);

	CanonicalMer_t source_mer;
	CanonicalMer_t sink_mer;
	CanonicalMer_t tmp_mer;
	
	Node_t * source_tmp;
	Node_t * sink_tmp;

	int source_offset = -1;
	int sink_offset = -1;
	int offset = -1;

	ref_m->trim5 = -1;
	ref_m->trim3 = -1;

	source_m = NULL;
	sink_m = NULL;
	
	bool ambiguous_match;
	
	// Find the first matching mer with sufficient coverage
	source_tmp = NULL;
	ambiguous_match = false;
	for (offset = 0; offset < (int) ref_m->rawseq.length(); ++offset)
	{	
		tmp_mer.set(ref_m->rawseq.substr(offset, K));
		source_tmp = getNode(tmp_mer);

		if ((source_tmp) && (source_tmp->getTotCov() >= COV_THRESHOLD) && (source_tmp->component_m == compid) )
		{ 
			if(source_m == NULL) { // found 1st match
				source_m = source_tmp;
				source_mer = tmp_mer;
				source_offset = offset;
			}
			else { // check if there is another match for the same mer
				if (source_m == source_tmp) { //found another identical match
					source_m = NULL;
					ambiguous_match = true;
					break;
				}
			}
			//source_m = source_tmp;
			//break; 
		}
	}

	if(ambiguous_match) {
		if(verbose) { cerr << "Ambiguous match to reference for source" << endl; }
		return;
	}
	
	if (!source_m) {
		if(verbose) { cerr << "No match to reference for source" <<  endl; }
		return;
	}

	// Find the last matching mer
	sink_tmp = NULL;
	ambiguous_match = false;
	for (offset = ref_m->rawseq.length()-K; offset >= 0; --offset)
	{
		tmp_mer.set(ref_m->rawseq.substr(offset, K));
		sink_tmp = getNode(tmp_mer);

		if ((sink_tmp) && (sink_tmp->getTotCov() >= COV_THRESHOLD) && (sink_tmp->component_m == compid))
		{
			if(sink_m == NULL) { // found 1st macth
				sink_m = sink_tmp;
				sink_mer = tmp_mer;
				sink_offset = offset;
			}
			else { // check if there is another match for the same mer
				if (sink_m == sink_tmp) { //found another identical match
					sink_m = NULL;
					ambiguous_match = true;
					break;
				}
			}
			//sink_m = sink_tmp;
			//break; 
		}
	}
	
	if(ambiguous_match) {
		if(verbose) { cerr << "Ambiguous match to reference for sink" << endl; }
		return;
	}
	
	if (!sink_m) {
		if(verbose) { cerr << "No match to reference for sink" << endl; }
		return;
	}

	int ref_dist = sink_offset - source_offset + K;
	sink_offset = ref_m->rawseq.length() - sink_offset - K;
	ref_m->seq = ref_m->rawseq.substr(source_offset, ref_dist);

	if(verbose) {
		cerr << "ref trim5: "     << source_offset 
		<< " trim3: "     << sink_offset
		<< " uncovered: " << source_offset + sink_offset
		<< " ref_dist: "  << ref_dist << endl;
	}

	ref_m->trim5 = source_offset;
	ref_m->trim3 = sink_offset;

	//cerr << " searching " << source_m->nodeid_m << " to " << sink_m->nodeid_m << endl;

	// Add the fake source node	
	std::stringstream sourceid;
	sourceid << "source" << compid;
	Node_t * newsource = new Node_t(sourceid.str());
	newsource->component_m = compid;

	Edgedir_t sourcedir = FF;
	if (source_mer.ori_m == R) { sourcedir = FR; }

	bool CLIP_REF_ENDS = true;
	if (CLIP_REF_ENDS)
	{
		if (VERBOSE) { cerr << "checking for other source edges" << endl; }
		for (int i = source_m->edges_m.size()-1; i >= 0; --i)
		{
			if (Edge_t::edgedir_start(source_m->edges_m[i].dir_m) == Edge_t::flipdir(source_mer.ori_m))
			{
				Node_t * other = getNode(source_m->edges_m[i]);
				
				if ( (other != NULL) && (other != source_m) ) { 

					if (VERBOSE) { cerr << "  removing node before source: " << other->nodeid_m << endl; }

					other->removeEdge(source_m->nodeid_m, Edge_t::fliplink(source_m->edges_m[i].dir_m));
					source_m->edges_m.erase(source_m->edges_m.begin() + i);
				}
			}
		}
	}

	newsource->addEdge(source_mer.mer_m, sourcedir, refid);
	newsource->setIsSource();
	source_m->addEdge(newsource->nodeid_m, Edge_t::fliplink(sourcedir), refid);
	source_m = newsource;

	nodes_m.insert(make_pair(newsource->nodeid_m, newsource));

	// Add the fake sink node	
	std::stringstream sinkid;
	sinkid << "sink" << compid;
	Node_t * newsink = new Node_t(sinkid.str());
	newsink->component_m = compid;
	
	Edgedir_t sinkdir = RR;
	if (sink_mer.ori_m == R) { sinkdir = FF; } 

	if (CLIP_REF_ENDS) 
	{
		if (VERBOSE) { cerr << "checking for other sink edges" << endl; }
		for (int i = sink_m->edges_m.size()-1; i >= 0; --i)
		{
			if (Edge_t::edgedir_start(sink_m->edges_m[i].dir_m) == sink_mer.ori_m)
			{
				Node_t * other = getNode(sink_m->edges_m[i]);
				
				if ( (other != NULL) && (other != sink_m) ) { 
				
					if (VERBOSE) { cerr << "  removing node after sink: " << other->nodeid_m << endl; }

					other->removeEdge(sink_m->nodeid_m, Edge_t::fliplink(sink_m->edges_m[i].dir_m));
					sink_m->edges_m.erase(sink_m->edges_m.begin() + i);
				}
			}
		}
	}

	newsink->addEdge(sink_mer.mer_m, sinkdir, refid);
	newsink->setIsSink();
	sink_m->addEdge(newsink->nodeid_m, Edge_t::fliplink(sinkdir), refid);
	sink_m = newsink;

	nodes_m.insert(make_pair(newsink->nodeid_m, newsink));
}

// markRefNodes
//////////////////////////////////////////////////////////////

void Graph_t::markRefNodes()
{
	if(verbose) { cerr << endl << "mark refnodes" << endl; }
	int nodes = 0;
	int refnodes = 0;

	MerTable_t::iterator mi;
	for (mi = nodes_m.begin(); mi != nodes_m.end(); ++mi)
	{
		++nodes;
		refnodes += mi->second->markRef(ref_m, K);
		mi->second->component_m = 0;
	}
	
	if(verbose) { cerr << " nodes: " << nodes << " refnodes: " << refnodes << endl; }
}

// numConnectedComponents
//////////////////////////////////////////////////////////////
int Graph_t::markConnectedComponents()
{
	if(verbose) { cerr << endl << "connected components" << endl; }
	int nodes = 0;
	int refnodes = 0;

	ref_m->refcompids.clear();

	MerTable_t::iterator mi;
	for (mi = nodes_m.begin(); mi != nodes_m.end(); ++mi)
	{
		++nodes;
		mi->second->component_m = 0;
	}

	int comp = 0;
	int refcomp = 0;

	for (mi = nodes_m.begin(); mi != nodes_m.end(); ++mi)
	{
		if (mi->second->component_m != 0) { continue; }

		++comp;

		deque<Node_t *> Q;
		Q.push_back(mi->second);

		int touches = 0;

		while (!Q.empty())
		{
			Node_t * cur = Q.front();
			Q.pop_front();

			if (cur->component_m == 0)
			{
				cur->component_m = comp;

				if (cur->touchRef_m) { ++touches; }

				for (unsigned int i = 0; i < cur->edges_m.size(); ++i)
				{
					Node_t * next = getNode(cur->edges_m[i]);
					Q.push_back(next);
				}
			}
			else
			{
				if (cur->component_m != comp)
				{
					cerr << "ERROR: cur->comp: " << cur->component_m << " != comp: " << comp << endl;
					exit(1);
				}
			}
		}

		if (touches)
		{
			++refcomp;
			ref_m->refcompids.insert(comp);
		}
	}


	ref_m->refnodes = refnodes;
	ref_m->refcomp  = refcomp;
	ref_m->allcomp  = comp;

	if (verbose) {
		cerr << " nodes: "    << nodes
		<< " refnodes: " << refnodes
		<< " comp: "     << comp
		<< " refcomp: "  << refcomp
		<< " refcompids: ";
	}

	set<int>::iterator si;
	for (si = ref_m->refcompids.begin(); si != ref_m->refcompids.end(); ++si)
	{
		if(verbose) { cerr << " " << *si; }
	}
	if(verbose) { cerr << endl; }
	
	return comp;
}

// denovoNodes
//////////////////////////////////////////////////////////////

void Graph_t::denovoNodes(const string & filename, const string & refname)
{
	cerr << endl << "finding denovo nodes" << endl;

	markRefNodes();

	FILE * fp = xfopen(filename, "w");

	MerTable_t::iterator mi;
	for (mi = nodes_m.begin(); mi != nodes_m.end(); ++mi)
	{
		Node_t * cur = mi->second;

		map<string, int> who;

		unordered_set<ReadId_t>::const_iterator si;
		for (si = cur->reads_m.begin(); si != cur->reads_m.end(); ++si)
		{
			string & set = readid2info[*si].set_m;
			if (set != "ref")
			{
				++(who[set]);
			}
		}

		if (who.size() == 1)
		{
			map<string, int>::iterator mi = who.begin();

			//char isRef = ref_m->isRefComp(cur->component_m) ? 'R' : 'N';

			if (VERBOSE)
			{
				fprintf(stderr, ">%s_%s len=%d cov=%0.02f comp=%d\n%s\n",
					mi->first.c_str(), 
					cur->nodeid_m.c_str(),
					cur->strlen(), 
					cur->getTotCov(),
					cur->component_m,
					cur->str_m.c_str());

			}

			fprintf(fp, ">%s_%s len=%d cov=%0.02f\n%s\n",
				mi->first.c_str(), 
				cur->nodeid_m.c_str(),
				cur->strlen(), 
				cur->getTotCov(),
				cur->str_m.c_str());
		}
	}

	xfclose(fp);
}

// alignRefNodes
//////////////////////////////////////////////////////////////

void Graph_t::alignRefNodes()
{
	int refpathnodes = 0;

	MerTable_t::iterator mi;
	for (mi = nodes_m.begin(); mi != nodes_m.end(); ++mi)
	{
		Node_t * cur = mi->second;

		if (cur->onRefPath_m)
		{
			++refpathnodes;
		}
	}

	if(verbose) { cerr << " Found " << refpathnodes << " on ref path" << endl; }
}

// countRefPath
//////////////////////////////////////////////////////////////

void Graph_t::countRefPath(const string & filename, const string & refname, bool printPathsToFile)
{
	if (source_m)
	{
		FILE * fp = NULL;
		
		if(printPathsToFile) {
			FILE * fp = xfopen(filename, "w");
		
			fprintf(fp, ">ref_raw\n%s\n",
				ref_m->rawseq.c_str());

			fprintf(fp, ">ref_trim %s trim5:%d trim3:%d\n%s\n", 
				refname.c_str(), ref_m->trim5, ref_m->trim3, ref_m->seq.c_str());
		}

		if (source_m != NULL && sink_m != NULL) {
			//dfs(source_m, sink_m, F, ref_m, fp, printPathsToFile);
			eka(source_m, sink_m, F, ref_m, fp, printPathsToFile);
		}
		alignRefNodes();

		if(printPathsToFile) { xfclose(fp); }
		
	}
}



// getNode
//////////////////////////////////////////////////////////////

Node_t * Graph_t::getNode(Mer_t nodeid)
{
	MerTable_t::iterator ni = nodes_m.find(nodeid);

	if (ni == nodes_m.end()) { return NULL; }

	return ni->second;
}

// getNode
//////////////////////////////////////////////////////////////

Node_t * Graph_t::getNode(CanonicalMer_t mer)
{
	return getNode(mer.mer_m);
}


// getNode
//////////////////////////////////////////////////////////////

Node_t * Graph_t::getNode(Edge_t & edge)
{
	MerTable_t::iterator ni = nodes_m.find(edge.nodeid_m);

	if (ni == nodes_m.end()) { return NULL; }

	return ni->second;
}


// compressNode
//////////////////////////////////////////////////////////////

void Graph_t::compressNode(Node_t * node, Ori_t dir)
{
	bool cnVERBOSE = 0;

	if (cnVERBOSE) { cerr << "Compress " << node->nodeid_m << ":" << dir << endl; }

	while (true)
	{
		if (cnVERBOSE) { cerr << node << endl; }

		int uniqueid = node->getBuddy(dir);
		if (uniqueid == -1) { return; }
		if (node->isTandem()) { return; }

		// make sure they are mutual buddies
		Edgedir_t edir = node->edges_m[uniqueid].dir_m;

		if (cnVERBOSE) { cerr << " --> " << node->edges_m[uniqueid] << endl; }

		Ori_t bdir = F;

		if (edir == FF || edir == RF) { bdir = R; }

		Node_t * buddy = getNode(node->edges_m[uniqueid]);
		if (buddy->isTandem()) { return; }

		if (!buddy)
		{
			cerr << "couldn't get " << node->nodeid_m << " - buddyid: " << node->edges_m[uniqueid] << endl;
		}

		assert(buddy);
		assert(!buddy->dead_m);

		if (cnVERBOSE)
		{
			cerr << "    " << buddy << endl;
		}

		int buniqueid = buddy->getBuddy(bdir);

		if (buniqueid == -1) 
		{ 
			if (cnVERBOSE) { cerr << "node's buddy's buddy is not me" << endl; }
			return; 
		}

		assert(buddy->edges_m[buniqueid].nodeid_m == node->nodeid_m);

		if (cnVERBOSE) { cerr << " --> " << node->edges_m[uniqueid] << endl; }

		// str
		string astr = node->str_m;
		if (dir == R)
		{ 
			astr = CanonicalMer_t::rc(astr); 
			node->revreads();
			node->revCovDistr();
		}

		string bstr = buddy->str_m;
		if (Edge_t::edgedir_dest(edir) == R)
		{ 
			bstr = CanonicalMer_t::rc(bstr); 
			buddy->revreads();
			buddy->revCovDistr();
		}

		if (cnVERBOSE)
		{
			string shift;
			shift.append(astr.length()-K+1, ' ');
			cerr << "     " << astr << endl; 
			cerr << "     " << shift << bstr << endl;
		}


		assert(astr.substr(astr.length()-K+1, K-1) == bstr.substr(0, K-1));

		string mstr = astr + bstr.substr(K-1);

		if (cnVERBOSE) { cerr << "     " << mstr << endl; }

		if (dir == R) 
		{ 
			mstr = CanonicalMer_t::rc(mstr); 
		}

		node->str_m = mstr;

		// coverage
		int   amerlen = astr.length() - K + 1;
		int   bmerlen = bstr.length() - K + 1;
		
		float ncov_tmr_fwd = node->getCov(FWD,TMR);
		float ncov_nml_fwd = node->getCov(FWD,NML);
		float ncov_tmr_rev = node->getCov(REV,TMR);
		float ncov_nml_rev = node->getCov(REV,NML);
		
		float ccov_tmr_fwd = buddy->getCov(FWD,TMR);
		float ccov_nml_fwd = buddy->getCov(FWD,NML);
		float ccov_tmr_rev = buddy->getCov(REV,TMR);
		float ccov_nml_rev = buddy->getCov(REV,NML);
		
		// add coverage info for the new base-pairs and update coverage of overlapping region	

		assert(node->cov_distr_tmr.size() == node->cov_distr_nml.size());
		
		/*
		int p = (node->cov_distr_tmr.size())-K+1;
		assert(p>0);
		for (int l = 0; l < K-1; l++) {
			// tumor
			int node_cov_tmr = node->cov_distr_tmr[p+l].fwd + node->cov_distr_tmr[p+l].rev;
			int buddy_cov_tmr = buddy->cov_distr_tmr[l].fwd + buddy->cov_distr_tmr[l].rev;
			if(node_cov_tmr < buddy_cov_tmr) { node->cov_distr_tmr[p+l] = buddy->cov_distr_tmr[l]; }
			// normal
			int node_cov_nml = node->cov_distr_nml[p+l].fwd + node->cov_distr_nml[p+l].rev;
			int buddy_cov_nml = buddy->cov_distr_nml[l].fwd + buddy->cov_distr_nml[l].rev;
			if(node_cov_nml < buddy_cov_nml) { node->cov_distr_nml[p+l] = buddy->cov_distr_nml[l]; }
		}
		
		int p = (node->cov_distr_tmr.size())-K+1;
		assert(p>0);
		for (int l = 0; l < K-1; l++) {
			node->cov_distr_tmr[p+l] = buddy->cov_distr_tmr[l];	// tumor		
			node->cov_distr_nml[p+l] = buddy->cov_distr_nml[l];	// tumor		
		}
		*/
		
		for (unsigned int j = (K-1); j < buddy->cov_distr_tmr.size(); ++j) {
			node->cov_distr_tmr.push_back(buddy->cov_distr_tmr[j]); // tumor
			node->cov_distr_nml.push_back(buddy->cov_distr_nml[j]); // normal
			node->cov_status.push_back(buddy->cov_status[j]);
		}
		

		// copy over barcodes
		node->bxset_tmr_fwd.insert(buddy->bxset_tmr_fwd.begin(), buddy->bxset_tmr_fwd.end());
		node->bxset_tmr_rev.insert(buddy->bxset_tmr_rev.begin(), buddy->bxset_tmr_rev.end());
		node->bxset_nml_fwd.insert(buddy->bxset_nml_fwd.begin(), buddy->bxset_nml_fwd.end());
		node->bxset_nml_rev.insert(buddy->bxset_nml_rev.begin(), buddy->bxset_nml_rev.end());
		
					
		node->computeMinCov(); // recompute min coverage 
		
		node->cov_tmr_m_fwd = ((ncov_tmr_fwd * amerlen) + (ccov_tmr_fwd * bmerlen)) / (amerlen + bmerlen);
		node->cov_nml_m_fwd = ((ncov_nml_fwd * amerlen) + (ccov_nml_fwd * bmerlen)) / (amerlen + bmerlen);

		node->cov_tmr_m_rev = ((ncov_tmr_rev * amerlen) + (ccov_tmr_rev * bmerlen)) / (amerlen + bmerlen);
		node->cov_nml_m_rev = ((ncov_nml_rev * amerlen) + (ccov_nml_rev * bmerlen)) / (amerlen + bmerlen);
		
		// reads
		unordered_set<ReadId_t>::const_iterator ri;
		for (ri = buddy->reads_m.begin(); ri != buddy->reads_m.end(); ++ri)
		{
			node->reads_m.insert(*ri);
		}

		// add buddy read starts
		int shift = amerlen;
		for (unsigned int i = 0; i < buddy->readstarts_m.size(); ++i)
		{
			ReadStart_t & rs = buddy->readstarts_m[i];
			node->readstarts_m.push_back(ReadStart_t(rs.readid_m, rs.nodeoffset_m+shift, rs.trim5_m, rs.ori_m));
		}

		// flip read starts
		if (dir == R)
		{
			node->revreads();
			node->revCovDistr();
		}

		node->sortReadStarts();

		// dead flag
		buddy->dead_m = true;

		// isRef
		//node->isRef_m |= buddy->isRef_m;
		if(buddy->isRef()) { node->setIsRef(); }
		if(buddy->isNormal()) { node->setIsNormal(); }
		if(buddy->isTumor()) { node->setIsTumor(); }
				
		// node edges
		node->edges_m.erase(node->edges_m.begin()+uniqueid);

		// now move over the buddy edges
		for (int i = 0; i < (int) buddy->edges_m.size(); ++i)
		{
			if (i == buniqueid) { continue; }

			Edge_t ne(buddy->edges_m[i]);

			if (edir == FR || edir == RF)
			{
				ne.dir_m = Edge_t::flipme(ne.dir_m);
			}

			Node_t * other = getNode(ne);

			if (cnVERBOSE) { cerr << "Keeping: " << ne << endl; }

			// NOTE: this comprassion can generate the error: 
			// scalpel: Path.cc:66: std::string Path_t::str(): Assertion `retval.substr(retval.length()-K+1) == nstr.substr(0, K-1)' failed.  
			if (other == buddy) 
			{
				cerr << "circle to buddy" << endl;
				ne.nodeid_m = node->nodeid_m;
				node->edges_m.push_back(ne);
			}
			else
			{
				node->edges_m.push_back(ne);
				other->updateEdge(buddy->nodeid_m, Edge_t::fliplink(buddy->edges_m[i].dir_m),
					node->nodeid_m, Edge_t::fliplink(ne.dir_m));
			}
		}		
	}
}


// compress
//////////////////////////////////////////////////////////////

void Graph_t::compress(int compid)
{
	if(verbose) { cerr << "compressing graph:"; }

	MerTable_t::iterator mi;

	for (mi = nodes_m.begin(); mi != nodes_m.end(); ++mi)
	{
		if(mi->second->component_m == compid) { //only analyze the selected connected component 
		
			if (mi->second->dead_m)  { continue; }
			//if (mi->second->isRef()) { continue; }
			if (mi->second->isSpecial()) { continue; }

			compressNode(mi->second, F);
			compressNode(mi->second, R);
		}
	}

	cleanDead();
}

// cleanDead
///////////////////////////////////////////////////////////////

void Graph_t::cleanDead()
{
	set<Mer_t> deadnodes;

	MerTable_t::iterator mi;
	for (mi = nodes_m.begin(); mi != nodes_m.end(); ++mi)
	{
		if (mi->second->dead_m) 
		{
			deadnodes.insert(mi->second->nodeid_m);
		}
	}

	if(verbose) { cerr << "  removing " << deadnodes.size() << " dead nodes" << endl; }

	set<Mer_t>::iterator di;
	for (di = deadnodes.begin(); di != deadnodes.end(); ++di)
	{
		mi = nodes_m.find(*di);
		assert(mi != nodes_m.end());
		assert(mi->second->dead_m);

		delete mi->second;
		nodes_m.erase(mi);
	}
}


// removeNode
//////////////////////////////////////////////////////////////

void Graph_t::removeNode(Node_t * node)
{
	assert(node);
	assert(!node->dead_m);

	node->dead_m = true;

	for(unsigned int i = 0; i < node->edges_m.size(); ++i)
	{
		Node_t * nn = getNode(node->edges_m[i]);

		if ((nn) && (nn != node))
		{
			nn->removeEdge(node->nodeid_m, Edge_t::fliplink(node->edges_m[i].dir_m));
		}
	}
}


// removeLowCov
//////////////////////////////////////////////////////////////

void Graph_t::removeLowCov(bool docompression, int compid)
{
	if(verbose) { cerr << endl << "removing low coverage:"; }

	int lowcovnodes = 0;
	
	double avgcov = ((double) totalreadbp_m) / ((double) ref_m->rawseq.length());
	//cerr << "avgcov: " << avgcov << endl;

	MerTable_t::iterator mi;

	for (mi = nodes_m.begin(); mi != nodes_m.end(); ++mi)
	{
		if(mi->second->component_m == compid) { //only process the selected connected component 
		
			Node_t * node = mi->second;

			//if (node->isRef())    { continue; }
			if (node->isSpecial())    { continue; }
			//if (node->touchRef_m) { continue; }

			//if ( (node->getMinCov() <= LOW_COV_THRESHOLD) || (node->getMinCov() <= (MIN_COV_RATIO*avgcov)) ||
			if ( (node->getMinCovMinQV() <= LOW_COV_THRESHOLD) || (node->getMinCovMinQV() <= (MIN_COV_RATIO*avgcov)) ||
				(node->getTotTmrCov() == 1 && node->getTotNmlCov() == 1) )
			{
				++lowcovnodes;
				removeNode(node);
			}
		}
	}

	if (verbose) { cerr << " found " << lowcovnodes; }

	cleanDead();
	if(docompression) { compress(compid); }

	if(verbose) { printStats(compid); }
}


// removeShortLinks
/////////////////////////////////////////////////////////////

void Graph_t::removeShortLinks(int compid)
{
	int links = 0;

	double avgcov = ((double) totalreadbp_m) / ((double) ref_m->rawseq.length());

	if (verbose) { cerr << endl << "remove short links: "; }

	MerTable_t::iterator mi;

	for (mi = nodes_m.begin(); mi != nodes_m.end(); ++mi)
	{
		if(mi->second->component_m == compid) { //only process the selected connected component 
		
			Node_t * cur = mi->second;

			//if (cur->isRef()) { continue; }
			if (cur->isSpecial()) { continue; }

			int deg = cur->edges_m.size();
			int len = cur->getSize();

			//if ((deg >= 2) && (len < MAX_LINK_LEN) )
			if ((deg >= 2) && (len < MAX_LINK_LEN) && (cur->getMinCov() <= floor(sqrt(avgcov))) )
			{
				int LEN=0;
				string MOTIF = "";
				//stringstream STR;
				
				// do not remove short-links within STRs: small bubbles are normal in STRs due to the DeBruijn graph represenation.   
				findTandems(cur->str_m, "shortlink", MAX_UNIT_LEN, MIN_REPORT_UNITS, MIN_REPORT_LEN, DIST_FROM_STR, K-1, LEN, MOTIF);
				//STR << LEN << MOTIF;
				//cout << cur->str_m << endl;
				//cout << "MS:" << STR.str() << endl;				
				if(LEN==0) {
					removeNode(cur);
					++links;
				}
			}
		}
	}

	if(verbose) { cerr << " removed links: " << links << endl; }

	if (links) { compress(compid); }

	if(verbose) { printStats(compid); }
}

// removeTips
//////////////////////////////////////////////////////////////

void Graph_t::removeTips(int compid)
{
	int tips = 0;
	int round = 0;

	do
	{
		++round;
		tips = 0;

		if (verbose) { cerr << endl << "remove tips round: " << round; }

		MerTable_t::iterator mi;

		for (mi = nodes_m.begin(); mi != nodes_m.end(); ++mi)
		{
			if(mi->second->component_m == compid) { //only process the selected connected component 
			
				Node_t * cur = mi->second;

				//if (cur->isRef()) { continue; }
				if (cur->isSpecial()) { continue; }

				int deg = cur->edges_m.size();
				int len = cur->strlen() - K + 1;

				if ((deg <= 1) && (len < MAX_TIP_LEN))
				{
					removeNode(cur);
					++tips;
				}
			}
		}

		if(verbose) { cerr << " removed: " << tips << endl; }

		if (tips) { compress(compid); }
	}
	while (tips);

	if(verbose) { printStats(compid); }
}

// greedyTrim
//////////////////////////////////////////////////////////////

class CovCmp
{
public:

  CovCmp(Graph_t * g) : _g(g) {}

  bool operator() (const string & a, const string & b)
  {
    MerTable_t::iterator ai = _g->nodes_m.find(a);
    MerTable_t::iterator bi = _g->nodes_m.find(b);

    return ai->second->getTotCov() > bi->second->getTotCov();
  }
  
  Graph_t * _g;
};

void Graph_t::greedyTrim()
{
	int branches = 0;

	cerr << endl << "greedy trim" << endl;

	MerTable_t::iterator mi;

  vector<string> nodelist;

	for (mi = nodes_m.begin(); mi != nodes_m.end(); ++mi)
	{
    nodelist.push_back(mi->first);
  }

  CovCmp covcmp(this);
  sort (nodelist.begin(), nodelist.end(), covcmp);

  for (unsigned int i = 0; i < nodelist.size(); ++i)
  {
    mi = nodes_m.find(nodelist[i]);

		Node_t * cur = mi->second;

    if (i < 10)
    {
      cerr << i << " " << nodelist[i] << " " << cur->getTotCov() << endl;
    }

    // figure out the best edges in the forward and reverse

    if (cur->edges_m.size() < 2) { continue; }

    Edge_t bestf(cur->edges_m[0]); float covf = -1; int degf = 0;
    Edge_t bestr(cur->edges_m[0]); float covr = -1; int degr = 0;

    for (unsigned int j = 0; j < cur->edges_m.size(); ++j)
    {
      Edge_t & edge = cur->edges_m[j];
      Node_t * other = getNode(edge);

      if (edge.isDir(F)) { ++degf; if (other->getTotCov() > covf) { bestf = edge; covf = other->getTotCov(); } }
      else               { ++degr; if (other->getTotCov() > covr) { bestr = edge; covr = other->getTotCov(); } }
    }

    // prune away all the other edges

    if (degf > 1 || degr > 1)
    {
      ++branches;

      for (unsigned int i = 0; i < cur->edges_m.size(); ++i)
      {
        Edge_t & edge = cur->edges_m[i];
        Node_t * other = getNode(edge);

        bool removeEdge = true;

        if (edge.isDir(F)) { if (edge.nodeid_m == bestf.nodeid_m && edge.dir_m == bestf.dir_m) { removeEdge = false; } }
        else               { if (edge.nodeid_m == bestr.nodeid_m && edge.dir_m == bestr.dir_m) { removeEdge = false; } }

        if (removeEdge)
        {
          other->removeEdge(cur->nodeid_m, Edge_t::fliplink(edge.dir_m));
        }
      }

      cur->edges_m.clear();

      if (covf != -1) { cur->edges_m.push_back(bestf); }
      if (covr != -1) { cur->edges_m.push_back(bestr); }
    }
	}

	cerr << " removed: " << branches << endl;

	if (branches) { compress(0); }

	if(verbose) { printStats(0); }
}



// threadReads
//////////////////////////////////////////////////////////////
/* 
**
**         1\         /3
**           \       /
**            ---N---
**           /       \
**         2/         \4
**
**  Check for individual reads that span 1-N-3, 1-N-4, 2-N-3, 2-N-4
**
*/

void Graph_t::threadReads(int compid)
{
    if (MIN_THREAD_READS == -1)
    {
      cerr << "Skipping threading reads" << endl;
      return;
    }

	MerTable_t::iterator mi;

	bool oldverbose = VERBOSE;

	//VERBOSE = true;

	int thread = 1;
	int threadround = 0;

	while (thread)
	{
		thread = 0;
		++threadround;

		cerr << endl << "threading round: " << threadround << " " << nodes_m.size() << " nodes" << endl;

		if (VERBOSE)
		{
			cerr << "================================================================" << endl;
			printGraph();

			char buffer [1024];
			sprintf(buffer, "thread_%d.dot", threadround);
			printDot(buffer, 0);

			cerr << "================================================================" << endl;
		}

		for (mi = nodes_m.begin(); mi != nodes_m.end(); ++mi)
		{
			
			if(mi->second->component_m == compid) { //only process the selected connected component 
			
				Node_t * cur = mi->second;

				//if (cur->isRef()) { continue; }
				if (cur->isSpecial()) { continue; }

				//if ((cur->degree(F) > 1) || (cur->degree(R) > 1))
				if ((cur->degree(F) > 1) && (cur->degree(R) > 1))
				{
					unsigned int numedges = cur->edges_m.size();

					vector<unsigned int> threaded;
					threaded.resize(numedges);

					for (unsigned int ti = 0; ti < numedges; ++ti)
					{
						threaded[ti] = 0;
					}

					if (VERBOSE)
					{
						cerr << endl;
						cerr << "=Checking: " << cur->nodeid_m << " " << cur->strlen() << "bp" << endl;
						cerr << "================================================================" << endl;

						for (unsigned int i = 0; i < numedges; ++i)
						{
							cerr << "  " << i << ": " << cur->nodeid_m << ":" << cur->edges_m[i] << endl;
						}
					}

					if (cur->isTandem()) 
					{ 
						if (VERBOSE) { cerr << " tandem, skipping" << endl << endl; }

						continue; 
					}

					// Find reads that span across the node

					for (unsigned int e1i = 0; e1i < numedges; ++e1i)
					{
						Edge_t & e1 = cur->edges_m[e1i];
						Node_t * n1 = getNode(e1);

						if (e1.startdir() != F) { continue; }

						set<ReadId_t> e1reads;
						for (unsigned int j = 0; j < e1.readids_m.size(); ++j)
						{
							e1reads.insert(e1.readids_m[j]);
						}

						for (unsigned int e2i = 0; e2i < numedges; ++e2i)
						{
							Edge_t & e2 = cur->edges_m[e2i];
							Node_t * n2 = getNode(e2);
						

							if (e2.startdir() != R) { continue; }

							// check if any reads span edges from e1 to e2
							//////////////////////////////////////////////

							set<ReadId_t> overlap;

							if (VERBOSE) cerr << e1i << ":" << e1.label() << " == " << e2i << ":" << e2.label() << " :";

							for (unsigned int j = 0; j < e2.readids_m.size(); ++j)
							{
								int jj = e2.readids_m[j];

								if (e1reads.find(jj) != e1reads.end())
								{
									if (VERBOSE) cerr << " " << jj;
									overlap.insert(jj);
								}
							}

							if (VERBOSE) cerr << endl;


							// check if any mates of n1 span to n2
							//////////////////////////////////////

							if (VERBOSE) cerr << "mates: " << endl;

							set<ReadId_t> mateoverlap;
							unordered_set<ReadId_t>::const_iterator s1, s2;

							for (s1 =  n1->reads_m.begin();
							s1 != n1->reads_m.end();
							++s1)
							{
								s2 = n2->reads_m.find(readid2info[*s1].mateid_m);

								if (s2 != n2->reads_m.end())
								{
									// mates are in n1 and n2
									// TODO: check orientation and spacing
									if (VERBOSE) cerr << " " << *s1 << " " << *s2 << endl;
								}
							}

							// overlap has reads that span from e1 (in F dir) to e2 (in R dir)

							if ( (int)overlap.size() >= MIN_THREAD_READS)
							{
								++(threaded[e1i]);
								++(threaded[e2i]);
							}
						}
					}

					// see what fraction of nodes were completely threaded

					unsigned int threadcnt = 0;
					for (unsigned int ti = 0; ti < numedges; ++ti)
					{
						if (threaded[ti] > 0)
						{
							++threadcnt;
						}
					}
				

					if (VERBOSE)
					{ 
						cerr << "threads through " << threadcnt << " of " << numedges << endl;
					}

					// if the node was completely resolved, go ahead and resolve it

					if (threadcnt == numedges)
					{
						if (VERBOSE) { cerr << endl << "All edges threaded, applying..." << endl; }

						int copy = 0;
						vector<Node_t*> newnodes;

						for (unsigned int e1i = 0; e1i < numedges; ++e1i)
						{
							Edge_t & e1 = cur->edges_m[e1i];

							if (e1.startdir() != F) { continue; }

							set<ReadId_t> e1reads;
							for (unsigned int j = 0; j < e1.readids_m.size(); ++j)
							{
								e1reads.insert(e1.readids_m[j]);
							}

							for (unsigned int e2i = 0; e2i < numedges; ++e2i)
							{
								Edge_t & e2 = cur->edges_m[e2i];

								if (e2.startdir() != R) { continue; }

								set<ReadId_t> overlap;

								for (unsigned int j = 0; j < e2.readids_m.size(); ++j)
								{
									int jj = e2.readids_m[j];

									if (e1reads.find(jj) != e1reads.end())
									{
										overlap.insert(jj);
									}
								}

								// overlap has reads that span from e1 (in F dir) to e2 (in R dir)

								if ((int)overlap.size() >= MIN_THREAD_READS)
								{
									++copy;

									char buffer [1024];
									sprintf(buffer, "%s_%d", cur->nodeid_m.c_str(), copy);

									if (VERBOSE)
									{
										cerr << "thread " << cur->nodeid_m << " " << cur->strlen() << "bp" << endl;
										cerr << "  1: " << e1 << endl
											<< "  2: " << e2 << endl
											<< "  r[" << overlap.size() << "]:";

										set<ReadId_t>::iterator si;
										for (si = overlap.begin(); si != overlap.end(); ++si) { cerr << " " << *si; }
										cerr << endl;

										cerr << "  making " << buffer << endl << endl;
									}

									// create new node and updated data structures accordingly 
									Node_t * copy = new Node_t(cur->nodeid_m);

									copy->nodeid_m = buffer;      // nodeid
									copy->setMinQV(MIN_QUAL_CALL);
									copy->setK(K);
									copy->str_m = cur->str_m;     // sequence
									copy->component_m = cur->component_m;
								
									copy->cov_tmr_m_fwd = overlap.size(); // tumor coverage forward
									copy->cov_tmr_m_rev = overlap.size(); // tumor coverage reverse
									copy->cov_nml_m_fwd = overlap.size(); // normal coverage forward
									copy->cov_nml_m_rev = overlap.size(); // normal coverage reverse
								
									copy->cov_distr_tmr.resize(copy->str_m.size());
									copy->cov_distr_nml.resize(copy->str_m.size());
									copy->cov_status = cur->cov_status; // T=tumot,N=normal,B=both,E=empty
									copy->cov_distr_tmr = cur->cov_distr_tmr;
									copy->cov_distr_nml = cur->cov_distr_nml;;
									//copy->updateCovDistr(overlap.size()); // coverage distribution

									if(cur->isRef()) { copy->setIsRef(); }
									if(cur->isNormal()) { copy->setIsNormal(); }
									if(cur->isTumor()) { copy->setIsTumor(); }
								
									// reads and edges
									set<ReadId_t>::iterator si;
									for (si = overlap.begin(); si != overlap.end(); ++si)
									{
										copy->addEdge(e1.nodeid_m, e1.dir_m, *si);
										copy->addEdge(e2.nodeid_m, e2.dir_m, *si);
									}

									// TODO: update read starts?

									copy->touchRef_m = cur->touchRef_m;

									newnodes.push_back(copy);
								}
							}
						}

						if (!newnodes.empty())
						{
							if (VERBOSE) 
							{ 
								cerr << "cleaning up: " << cur->nodeid_m << endl; 

								for (unsigned int e1i = 0; e1i < cur->edges_m.size(); ++e1i)
								{
									Edge_t & e1 = cur->edges_m[e1i];
									cerr << "   " << e1 << endl;
								} 
							} 

							removeNode(cur);

							++thread;

							for (unsigned int j = 0; j < newnodes.size(); ++j)
							{
								Node_t * nn = newnodes[j];
								nodes_m.insert(make_pair(nn->nodeid_m, nn));

								if (VERBOSE) { cerr << "  swapping in: " << nn->nodeid_m << endl; }

								for (unsigned int k = 0; k < nn->edges_m.size(); ++k)
								{
									Edge_t & e = nn->edges_m[k];
									Node_t * other = getNode(e.nodeid_m);

									if (VERBOSE) { cerr << "    edge: " << e << endl; }

									for (unsigned int r = 0; r < e.readids_m.size(); ++r)
									{
										other->addEdge(nn->nodeid_m, Edge_t::fliplink(e.dir_m), e.readids_m[r]);
									}

									if (VERBOSE) { cerr << "      " << other << endl; }
								}
							}

							if (VERBOSE) { cerr << endl << endl; }
						}
					}
					else
					{
						if (VERBOSE) { cerr << "Not all edges threaded, skipping" << endl << endl; }
					}
				}
			}
		}
		
		cerr << "  round: " << threadround << " threaded: " << thread << endl;

		if (thread > 0)
		{
			cleanDead();
			compress(compid);
		}
	}

	VERBOSE=oldverbose;

	if(verbose) { printStats(0); }
}

// checkReadStarts
//////////////////////////////////////////////////////////////

void Graph_t::checkReadStarts(int compid)
{
	if(verbose) { cerr << "checking read starts.... "; }

	int all = 0;
	int bad = 0;

	MerTable_t::iterator mi;

	for (mi = nodes_m.begin(); mi != nodes_m.end(); ++mi)
	{
		if(mi->second->component_m == compid) { //only analyze the selected connected component 

			Node_t * cur = mi->second;

			for (unsigned int i = 0; i < cur->readstarts_m.size(); ++i)
			{
				ReadStart_t & rstart = cur->readstarts_m[i];
				ReadId_t rid         = rstart.readid_m;
				ReadInfo_t & rinfo   = readid2info[rid];

				string ckmer;
				string rkmer = rinfo.seq_m.substr(rstart.trim5_m, K);

				++all;

				if (rstart.ori_m == R)
				{
					ckmer = cur->str_m.substr(rstart.nodeoffset_m-K+1, K);
					ckmer = CanonicalMer_t::rc(ckmer);
				}
				else
				{
					ckmer = cur->str_m.substr(rstart.nodeoffset_m, K);
				}

				if ((rkmer != ckmer)) // || VERBOSE)
				{
					cerr << "Checking " << rid << " " << rinfo.readname_m 
						<< " " << rstart.ori_m 
						<< " offset:" << rstart.nodeoffset_m 
						<< " trim5:" << rstart.trim5_m << endl;
					cerr << "  " << rkmer << endl;
					cerr << "  " << ckmer << endl;

					cur->print(cerr) << endl;

					if (rkmer != ckmer)
					{
						++bad;
						cerr << "mismatch: " << cur->str_m << endl;
					}
					else
					{
						cerr << "ok" << endl;
					}
				}
			}
		}
	}

	if(verbose) { cerr << " found " << bad << " bad starts out of " << all << endl; }
}

// updateContigReadStarts
//////////////////////////////////////////////////////////////

void Graph_t::updateContigReadStarts()
{
	MerTable_t::iterator mi;

	for (mi = nodes_m.begin(); mi != nodes_m.end(); ++mi)
	{
		Node_t * cur = mi->second;

		for (unsigned int i = 0; i < cur->readstarts_m.size(); ++i)
		{
			ReadId_t rid = cur->readstarts_m[i].readid_m;

			readid2info[rid].contigid_m = cur->nodeid_m;
			readid2info[rid].readstartidx_m = i;
		}
	}
}

// bundleMates
//////////////////////////////////////////////////////////////

void Graph_t::bundleMates()
{
	MerTable_t::iterator mi;

	// compute contig links

	int links = 0;

	for (mi = nodes_m.begin(); mi != nodes_m.end(); ++mi)
	{
		Node_t * cur = mi->second;

		for (unsigned int i = 0; i < cur->readstarts_m.size(); ++i)
		{
			ReadId_t mateid = readid2info[cur->readstarts_m[i].readid_m].mateid_m;

			if (mateid != -1)
			{
				Mer_t matecontig = readid2info[mateid].contigid_m;

				if (matecontig != "")
				{
					cur->addContigLink(matecontig, cur->readstarts_m[i].readid_m);
					++links;
				}
			}
		}
	}

	cerr << "created " << links << " contig links" << endl;


	// bundle consistent links

	for (mi = nodes_m.begin(); mi != nodes_m.end(); ++mi)
	{
		Node_t * cur = mi->second;

		if (cur->contiglinks_m.size() != 0)
		{
			// cur is linked to at least 1 other contig

			ContigLinkMap_t::iterator mi;

			for (mi =  cur->contiglinks_m.begin();
			mi != cur->contiglinks_m.end();
			++mi)
			{
				ContigLinkList_t * list = mi->second;

				Node_t * other = getNode(mi->first);

				if (0)
				{
					cerr << cur->nodeid_m << "(" << cur->strlen() <<  ") -- " 
						<< list->linkCnt() << " -- " << mi->first
						<< "(" << other->strlen() << ")" << endl;
				}

				ContigLinkList_t bundles [4];

				int internal = 0;
				int invalid = 0;

				int lastpos = -1;
				int lastmpos = -1;

				for (unsigned int i = 0; i < list->linkCnt(); ++i)
				{
					ContigLink_t & link = list->linklist_m[i];

					ReadId_t rid         = link.rid_m;
					ReadInfo_t & rinfo   = readid2info[rid];
					ReadStart_t & rstart = cur->readstarts_m[rinfo.readstartidx_m];

					ReadId_t mid         = rinfo.mateid_m;
					ReadInfo_t & minfo   = readid2info[mid];
					ReadStart_t & mstart = other->readstarts_m[minfo.readstartidx_m];

					Edgedir_t linkdir = FF;

					if      ((rstart.ori_m == F) && (mstart.ori_m == F)) { linkdir = FF; }
					else if ((rstart.ori_m == F) && (mstart.ori_m == R)) { linkdir = FR; }
					else if ((rstart.ori_m == R) && (mstart.ori_m == F)) { linkdir = RF; }
					else if ((rstart.ori_m == R) && (mstart.ori_m == R)) { linkdir = RR; }

					int linkdist = 0;

					if (cur == other)
					{
						++internal;

						if      (linkdir == FR) { linkdist = mstart.nodeoffset_m - rstart.nodeoffset_m; }
						else if (linkdir == RF) { linkdist = rstart.nodeoffset_m - mstart.nodeoffset_m; }
						else
						{
							++invalid;
						}
					}
					else
					{
						int adist = cur->strlen() - rstart.nodeoffset_m;
						if (rstart.ori_m == R) { adist = rstart.nodeoffset_m; }

						int bdist = other->strlen() - mstart.nodeoffset_m;
						if (mstart.ori_m == R) { bdist = mstart.nodeoffset_m; }

						linkdist = INSERT_SIZE - (adist + bdist);
					}

					int dup = 0;

					if ((rstart.nodeoffset_m == lastpos) && (mstart.nodeoffset_m == lastmpos))
					{
						dup = 1;
					}

					lastpos  = rstart.nodeoffset_m;
					lastmpos = mstart.nodeoffset_m;

					if (0)
					{
						int lo = linkdist - 2 * INSERT_STDEV;
						int hi = linkdist + 2 * INSERT_STDEV;

						cerr << Edge_t::toString(linkdir) << ":" << linkdist << "\t"
							<< lo << "\t" << hi << "\t"
							<< rinfo.code_m << "\t" << dup << "\t"
							<< rid << "\t"  << rinfo.readname_m <<  "\t" << rstart.nodeoffset_m << "\t" << rstart.ori_m << "\t"
							<< mid << "\t"  << minfo.readname_m <<  "\t" << mstart.nodeoffset_m << "\t" << mstart.ori_m << endl;
					}

					if (!dup)
					{
						bundles[linkdir].addLink(rid, linkdist);
					}
					else
					{
						bundles[linkdir].addDup();
					}
				}

				//cerr << "--" << endl;

				for (int i = 0; i < 4; ++i)
				{
					int cnt = bundles[i].linkCnt();
					int dups = bundles[i].dupCnt();

					if (cnt > 0)
					{
						float mean  = bundles[i].mean();
						float stdev = bundles[i].stdev(mean);

						cerr << cur->nodeid_m << "(" << cur->strlen() <<  ") "
							<< Edge_t::toString((Edgedir_t) i) << " "
							<< mi->first << "(" << other->strlen() << ")\t" 
							<< mean << "\t+/-\t" << stdev << "\t"
							<< cnt << "\t(" << dups << ")" << endl;
					}
				}

				//cerr << endl;
			}
		}
	}
}


// scaffold
//////////////////////////////////////////////////////////////

void Graph_t::scaffoldContigs()
{
	cerr << endl << "== scaffolding ==" << endl;

	checkReadStarts(0);
	updateContigReadStarts();
	bundleMates();

	cerr << endl << "==" << endl << endl;
}

// printGraph
//////////////////////////////////////////////////////////////

void Graph_t::printGraph()
{
	MerTable_t::iterator gi;
	for (gi = nodes_m.begin(); gi != nodes_m.end(); ++gi)
	{
		cout << gi->second << endl;
	}
}

// printStats
//////////////////////////////////////////////////////////////

void Graph_t::printStats(int compid)
{
	int edgecnt = 0;
	int span = 0;

	MerTable_t::iterator gi;
	for (gi = nodes_m.begin(); gi != nodes_m.end(); ++gi)
	{
		if(gi->second->component_m == compid) {
			edgecnt += gi->second->edges_m.size();
			span += gi->second->strlen();
		}
	}

	cerr << "  " << compid << ": nodes: " << nodes_m.size() 
		<< " edges: " << edgecnt
		<< " span: " << span << endl;
}

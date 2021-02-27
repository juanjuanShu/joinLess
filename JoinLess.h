#pragma once
#include "stdafx.h"
#include "Types.h"
#include "Common.h"

class JoinLess
{
public:
	JoinLess(vector<InstanceType>& instances, double min_prev, double min_cond_prob,double star_radius);
	void execute();

private:
	double _min_prev;
	double _min_cond_prob;
	double _star_radius;
	map<FeatureType,map<InstanceIdType,LocationType>> _instances;
	map<unsigned int,ColocationPackage> _prevalentColocation;
	map<FeatureType, map<InstanceIdType,NeighborhoodsObjectType>> _starNeighborhoods;
	RelationsType _relations;

	RelationsType _gen_relations(vector<InstanceType>& instances);
	void _gen_star_neighborhoods();
	ColocationSetType  _generateCandidateColocations_k(int k);
	ColocationSetType  _generateCandidateColocations_2();
	bool  _isSubsetPrevalent(ColocationType& candidates, int k);
	void  _filterStarInstances(ColocationSetType candidates, int k);
	void  _gen_star_k_instances(NeighborhoodsObjectType neighborhoodsObject, FeatureType centerFeature£¬ int k);
};


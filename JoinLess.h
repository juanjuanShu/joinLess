#pragma once
#include "stdafx.h"
#include "Types.h"
#include "Common.h"

struct Rule {
	ColocationType antecedent;
	ColocationType consequent;
	double conf;
};

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
	map<unsigned int, ColocationSetType> _prevalentColocation;
	map<FeatureType, map<InstanceIdType,NeighborhoodsObjectType>> _starNeighborhoods;
	RelationsType _relations;
	map<FeatureType, unsigned int> _numOfInstances;
	map<unsigned int,map<ColocationType,TableRowNeighborhoodsType>> _cliqueInstances;

	RelationsType _gen_relations(vector<InstanceType>& instances);
	void _gen_star_neighborhoods();
	ColocationSetType  _generateCandidateColocations_k(int k);
	ColocationSetType  _generateCandidateColocations_2();
	bool  _isSubsetPrevalent(ColocationType& candidates, int k);
	map<ColocationType, TableRowNeighborhoodsType>  _filterStarInstances(ColocationSetType candidates, int k);

	 void _gen_star_k_instances(
		  map<ColocationType, TableRowNeighborhoodsType> &ans,
		  const ColocationSetType& candidates,
		const NeighborhoodsObjectType& neighborhoodsObject, 
		FeatureType centerFeature, 
		int k);
	
	void _gen_star_k_instances_recursion(
		map<ColocationType, TableRowNeighborhoodsType>& ans,
		const NeighborhoodsObjectType& neighborhoodsObject,
		int k,
		int pos,
		int remainder,//还要填充的位置有几个
		ColocationType& tmp_colocation,
		NeighborhoodsObjectType &tmp_neighborhoods_object, 
		const ColocationSetType& candidates
	);

	void _select_coarse_prevalent_colocations(map<ColocationType, TableRowNeighborhoodsType> &colocationInstanceMap);

	bool isPrevalentByPartiIndex(const ColocationType& colocations, const TableRowNeighborhoodsType& tableRowNeighborhoods);

	void _filterCliqueInstances(map<ColocationType, TableRowNeighborhoodsType>& colocationInstanceMap, int k);

	bool isClique(NeighborhoodsObjectType neighborhoodsObjectType, int k);

	void _selectPrevalentColocations(int k);

	bool isSubCandidate(const ColocationType& tmp_colocation, const ColocationSetType& candidates);

	unsigned int  getRowInstancesOfColocationSub(const ColocationType& colocationSub);
	vector<Rule> _generateRules();
};


#include "JoinLess.h"
#include "Common.h"

JoinLess::JoinLess(vector<InstanceType>& instances,double min_prev, double min_cond_prob, double star_radius)
	:_min_prev(min_prev),
	_min_cond_prob(min_cond_prob),
	_star_radius(star_radius)
{
	for (auto it = instances.begin(); it != instances.end(); it++) {
		auto instanceId = get<InstanceIdType>(*it);
		auto feature = get<FeatureType>(*it);
		auto location = get<LocationType>(*it);

		_instances[feature][instanceId] = location;

		_prevalentColocation[1][{feature}].push_back({ instanceId });
	}

	_relations = _gen_relations(instances);
}

RelationsType JoinLess::_gen_relations(vector<InstanceType>& instances) {
	RelationsType relations;

	for (unsigned int i = 0; i < instances.size(); ++i) {
		InstanceType& instance1 = instances[i];
		auto id1 = get<InstanceIdType>(instance1);
		auto feature1 = get<FeatureType>(instance1);
		auto &location1 = get<LocationType>(instance1);
		for (unsigned int j = i + 1; j < instances.size(); ++j) {
			InstanceType& instance2 = instances[j];
			auto id2 = get<InstanceIdType>(instance2);
			auto feature2 = get<FeatureType>(instance2);
			auto location2 = get<LocationType>(instance2);
			
			//这样做的前提是输入数据的feature有序
			//if (feature1 >= feature2) continue;

			//这样做也是错误的，因为边和边之间可以是同类关系，因为最终合成的是中心和边的关系
			//if (feature1 == feature2) continue;
			
			if ((new Common())->hasRelation(location1, location2,_star_radius)) {
				// feature1.id1 should less than feature2.id2.这里未做处理，是因为默认输入数据有序
				_relations.insert({ { feature1 ,id1 }, { feature2 ,id2 } });
			}
		}
		
	}

	return relations;
}
	
void JoinLess::_gen_star_neighborhoods() {
	for (auto &relation : _relations) {
		auto &starCenter = relation.first;
		auto &starEdge = relation.second;

		// The star neighborhoods include the star center self.
		FeatureType centerFeature = starCenter.first;
		InstanceIdType centerId = starCenter.second;
		if (!_starNeighborhoods.count(centerFeature) || !_starNeighborhoods[centerFeature].count(centerId)) {
			_starNeighborhoods[centerFeature][centerId].push_back(starCenter);
		}

		_starNeighborhoods[centerFeature][centerId].push_back(starEdge);
	}
}

ColocationSetType JoinLess::_generateCandidateColocations_2() {
	vector<FeatureType> colocations;
	vector<ColocationType> candidateColocations;

	auto &prevalent = _prevalentColocation[1];
	//获取到实例类型，排序
	for (auto it_data = prevalent.begin(); it_data != prevalent.end(); it_data++) {
		colocations.push_back((*it_data).first[0]);
	}
	sort(colocations.begin(), colocations.end());
	//A B C
	for (unsigned int i = 0; i < colocations.size() - 1; i++) {
		for (unsigned int j = i + 1; j < colocations.size(); j++) {
			candidateColocations.push_back({ colocations[i],colocations[j] });
		}
	}

	return candidateColocations;
}

bool  JoinLess::_isSubsetPrevalent(ColocationType& candidates, int k) {
	if (k <= 2) return true;

	for (unsigned int i = 0; i < candidates.size(); i++) {
		ColocationType candidatesCopy(candidates);
		candidatesCopy.erase(candidatesCopy.begin() + i);
		if (!_prevalentColocation[k - 1].count(candidatesCopy)) {
			return false;
		}
	}

	return true;
}

ColocationSetType JoinLess::_generateCandidateColocations_k(int k)
{
	if (k == 2)  return _generateCandidateColocations_2();
	
	vector<ColocationType> candidateColocations;
	
	//get
	ColocationSetType C;
	ColocationPackage& colocationPackage = _prevalentColocation[k - 1];
	for (auto it = colocationPackage.begin(); it != colocationPackage.end(); it++) {
		C.push_back((*it).first);
	}
	sort(C.begin(), C.end());

	//存储
	vector<FeatureType> colocationSet;
	map < ColocationType, ColocationType> trie_colocationSet = {};
	for (unsigned int i = 0; i < C.size(); ++i) {
		colocationSet = C[i];
		FeatureType lastElement = colocationSet.back();
		colocationSet.pop_back();
		if (trie_colocationSet.find(colocationSet) == trie_colocationSet.end()) {
			trie_colocationSet.insert({ colocationSet,{lastElement} });
		}
		else {
			trie_colocationSet[colocationSet].push_back(lastElement);
		}
	}

	//连接
	for (auto& item : trie_colocationSet) {
		ColocationType candidate = item.first;
		//如果后面的k只有一个，则无法连接
		if (item.second.size() >= 2) {
			for (auto it_value = (item.second).begin(); it_value != (item.second).end() - 1; it_value++) {
				for (auto it_value1 = it_value + 1; it_value1 != (item.second).end(); it_value1++) {
					candidate.push_back(*it_value);
					candidate.push_back(*it_value1);
					if (_isSubsetPrevalent(candidate, k)) {
						candidateColocations.push_back(candidate);
					}
				}
			}
		}
	}

	return candidateColocations;
}

void JoinLess::_filterStarInstances(ColocationSetType candidates, int k) {
	//收集所有候选colocation的第一个Feature 
	set<FeatureType> centerFeatureSet;
	for (auto& candidate : candidates) {
		centerFeatureSet.insert(candidate[0]);
	}

	//候选colocation的第一个Feature一定是星型的中心
	for (auto centerFeature : centerFeatureSet) {
		auto& idStarNeighborhoods = _starNeighborhoods[centerFeature];
		for (auto idStarNeighborhood : idStarNeighborhoods) {
			InstanceIdType centerId = idStarNeighborhood.first;
			NeighborhoodsObjectType neighborhoodsObject = idStarNeighborhood.second;

			if (neighborhoodsObject.size() < k) continue;

			//auto colocationInstanceMap = _gen_star_k_instances(neighborhoodsObject,k);
			//A1 B1 C1 生成k的表实例
			_gen_star_k_instances(neighborhoodsObject, centerFeature,k);
			//得是candidates里的值

		}
	}
	
}

void _gen_star_k_instances_recursion(
	map<ColocationType, TableRowNeighborhoodsType> &ans,
	const NeighborhoodsObjectType& neighborhoodsObject,
	int k,
	int pos,
	int remainder,//还要填充的位置有几个
	ColocationType &tmp_colocation,
	NeighborhoodsObjectType tmp_neighborhoods_object
) {
	if (pos + remainder > neighborhoodsObject.size()) return;

	if (0 == remainder) {
		ans[tmp_colocation].push_back(tmp_neighborhoods_object);
	}

	tmp_neighborhoods_object[k - remainder] = neighborhoodsObject[pos];
	tmp_colocation[k - remainder] = neighborhoodsObject[pos].first;

	_gen_star_k_instances_recursion(ans, neighborhoodsObject, k, pos + 1, remainder - 1,tmp_colocation, tmp_neighborhoods_object);
}

void JoinLess::_gen_star_k_instances(const NeighborhoodsObjectType &neighborhoodsObject,FeatureType centerFeature, int k)
{
	//NeighborhoodsObjectType
	map<ColocationType, TableRowNeighborhoodsType> ans;
	ColocationType tmp_colocation;
	NeighborhoodsObjectType tmp_neighborhoods_object;

	//中心Feature是第一个元素
	tmp_colocation[0] = centerFeature;
	tmp_neighborhoods_object[0] = neighborhoodsObject[0];

	_gen_star_k_instances_recursion(ans, neighborhoodsObject,k,1,k - 1, 
		tmp_colocation, tmp_neighborhoods_object);
}

void JoinLess::execute()
{
	_gen_star_neighborhoods();
	int k = 2;

	while (_prevalentColocation.count(k - 1) && !_prevalentColocation[k - 1].empty()) {
		ColocationSetType candidates = _generateCandidateColocations_k(k);
		_filterStarInstances(candidates,k);
		k++;
	}

}

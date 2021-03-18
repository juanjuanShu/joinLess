#include "JoinLess.h"
#include "Common.h"

JoinLess::JoinLess(vector<InstanceType>& instances, double min_prev, double min_cond_prob, double star_radius)
	:_min_prev(min_prev),
	_min_cond_prob(min_cond_prob),
	_star_radius(star_radius)
{
	for (auto it = instances.begin(); it != instances.end(); it++) {
		auto instanceId = get<InstanceIdType>(*it);
		auto feature = get<FeatureType>(*it);
		auto& location = get<LocationType>(*it);

		_instances[feature][instanceId] = location;

		_prevalentColocation[1].push_back({ {feature} });

		//_colocationNum (特征，数量)
		auto ret = _numOfInstances.insert({ feature,1 });
		if (!ret.second) { ++ret.first->second; }

		_cliqueInstances[1][{feature}].push_back({ make_pair(feature, instanceId) });
	}

	_gen_relations(instances);
}

void JoinLess::_gen_relations(vector<InstanceType>& instances) {
	for (unsigned int i = 0; i < instances.size(); ++i) {
		InstanceType& instance1 = instances[i];
		auto id1 = get<InstanceIdType>(instance1);
		auto feature1 = get<FeatureType>(instance1);
		auto& location1 = get<LocationType>(instance1);
		for (unsigned int j = i + 1; j < instances.size(); ++j) {
			InstanceType& instance2 = instances[j];
			auto id2 = get<InstanceIdType>(instance2);
			auto feature2 = get<FeatureType>(instance2);
			auto& location2 = get<LocationType>(instance2);

			//这样做的前提是输入数据的feature有序
			//if (feature1 >= feature2) continue;

			//这样做也是错误的，因为边和边之间可以是同类关系，因为最终合成的是中心和边的关系
			//if (feature1 == feature2) continue;

			if ((new Common())->hasRelation(location1, location2, _star_radius)) {
				// feature1.id1 should less than feature2.id2.这里未做处理，是因为默认输入数据有序
				_relations.insert({ { feature1 ,id1 }, { feature2 ,id2 } });
			}
		}

	}
}

void JoinLess::_gen_star_neighborhoods() {
	for (auto& relation : _relations) {
		auto& starCenter = relation.first;
		auto& starEdge = relation.second;

		// The star neighborhoods include the star center self.
		FeatureType centerFeature = starCenter.first;
		InstanceIdType centerId = starCenter.second;
		if (!_starNeighborhoods.count(centerFeature) || !_starNeighborhoods[centerFeature].count(centerId)) {
			_starNeighborhoods[centerFeature][centerId].push_back(starCenter);
		}

		FeatureType edgeFeature = starEdge.first;
		if (edgeFeature != centerFeature) {
			_starNeighborhoods[centerFeature][centerId].push_back(starEdge);
		}
	}
}

ColocationSetType JoinLess::_generateCandidateColocations_2() {
	vector<ColocationType> candidateColocations;
	ColocationType colocations;

	for (auto& numOfInstance : _numOfInstances) {
		FeatureType feature = numOfInstance.first;
		colocations.push_back(feature);
	}

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

	auto& prevalent = _prevalentColocation[k - 1];
	for (unsigned int i = 0; i < candidates.size(); i++) {
		ColocationType candidatesCopy(candidates);
		candidatesCopy.erase(candidatesCopy.begin() + i);
		if (!isSubCandidate(candidatesCopy, prevalent)) {
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
	ColocationSetType  C = _prevalentColocation[k - 1];

	//存储
	vector<FeatureType> colocationSet;
	map < ColocationType, ColocationType> trie_colocationSet = {};
	for (auto it = C.begin(); it != C.end(); it++) {
		colocationSet = (*it);
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
					ColocationType tmpCandidate(candidate);
					tmpCandidate.push_back(*it_value);
					tmpCandidate.push_back(*it_value1);
					if (_isSubsetPrevalent(tmpCandidate, k)) {
						candidateColocations.push_back(tmpCandidate);
					}
				}
			}
		}
	}

	return candidateColocations;
}

map<ColocationType, TableRowNeighborhoodsType> JoinLess::_filterStarInstances(ColocationSetType candidates, int k) {
	map<ColocationType, TableRowNeighborhoodsType> colocationInstanceMap;

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

			//对于每一个中心生成的表实例
			_gen_star_k_instances(colocationInstanceMap, candidates, neighborhoodsObject, centerFeature, k);
		}
	}

	return colocationInstanceMap;
}

bool JoinLess::isSubCandidate(const ColocationType& tmp_colocation, const ColocationSetType& candidates) {
	return binary_search(candidates.begin(), candidates.end(), tmp_colocation);
}

void JoinLess::_gen_star_k_instances_recursion(
	map<ColocationType, TableRowNeighborhoodsType>& ans,
	const NeighborhoodsObjectType& neighborhoodsObject,
	int k,
	int pos,
	int remainder,//还要填充的位置有几个
	ColocationType& tmp_colocation,
	NeighborhoodsObjectType& tmp_neighborhoods_object,
	const ColocationSetType& candidates
) {
	if (pos + remainder > neighborhoodsObject.size()) return;

	if (0 == remainder) {
		//进一步判读是不是是candidates的子集
		if (isSubCandidate(tmp_colocation, candidates)) {
			ans[tmp_colocation].push_back(tmp_neighborhoods_object);
		}
		return;
	}

	//取该元素
	tmp_neighborhoods_object[k - remainder] = neighborhoodsObject[pos];
	tmp_colocation[k - remainder] = neighborhoodsObject[pos].first;
	_gen_star_k_instances_recursion(ans, neighborhoodsObject, k, pos + 1, remainder - 1, tmp_colocation, tmp_neighborhoods_object, candidates);

	//不取该元素，remainder的位置保留不变
	_gen_star_k_instances_recursion(ans, neighborhoodsObject, k, pos + 1, remainder, tmp_colocation, tmp_neighborhoods_object, candidates);
}

void JoinLess::_gen_star_k_instances(map<ColocationType, TableRowNeighborhoodsType>& ans, const ColocationSetType& candidates,
	const NeighborhoodsObjectType& neighborhoodsObject,
	FeatureType centerFeature,
	int k)
{
	NeighborhoodsObjectType newNeighborhoodsObject;

	//预分配空间
	ColocationType tmp_colocation(k);
	NeighborhoodsObjectType tmp_neighborhoods_object(k);

	//中心Feature是第一个元素
	tmp_colocation[0] = centerFeature;
	tmp_neighborhoods_object[0] = neighborhoodsObject[0];

	_gen_star_k_instances_recursion(ans, neighborhoodsObject, k, 1, k - 1, tmp_colocation, tmp_neighborhoods_object, candidates);

}

bool JoinLess::isPrevalentByPartiIndex(const ColocationType& colocations, const TableRowNeighborhoodsType& tableRowNeighborhoods) {
	bool isPrevalent = true;

	//初始化位图
	map<FeatureType, vector<bool>> bitMap;
	for (unsigned int i = 0; i < colocations.size(); i++) {
		FeatureType feature = colocations[i];
		//numOfInstances[feature]是feature的实例数
		bitMap[feature] = vector<bool>(_numOfInstances[feature], false);
	}

	//赋值 
	// A B ;对于A 
	// 0 1
	for (unsigned int i = 0; i < colocations.size(); i++) {
		FeatureType feature = colocations[i];
		for (auto rowInstance : tableRowNeighborhoods) {
			InstanceIdType id = get<InstanceIdType>(rowInstance[i]);
			bitMap[feature][id - 1] = true;
		}
	}

	//计算
	for (auto it_bit = bitMap.begin(); it_bit != bitMap.end(); it_bit++) {
		FeatureType feature = (*it_bit).first;
		vector<bool> flag = (*it_bit).second;

		int count = 0;
		for (unsigned int i = 0; i < flag.size(); i++) {
			if (flag[i]) count++;
		}

		if (count * 1.0 / flag.size() < _min_prev) {
			isPrevalent = false;
			break;
		}
	}

	return isPrevalent;
}

void JoinLess::_select_coarse_prevalent_colocations(map<ColocationType, TableRowNeighborhoodsType>& colocationInstanceMap) {
	//根据参与度进行剪枝
	for (auto it = colocationInstanceMap.begin(); it != colocationInstanceMap.end();) {
		auto& colocations = (*it).first;
		auto& tableRowNeighborhoods = (*it).second;

		if (!isPrevalentByPartiIndex(colocations, tableRowNeighborhoods)) {
			it = colocationInstanceMap.erase(it);
		}
		else {
			++it;
		}
	}
}

bool JoinLess::isClique(NeighborhoodsObjectType neighborhoodsObject, int k) {
	bool flag = false;

	//A1 B2 C1 => B2 C1
	auto it = neighborhoodsObject.begin();
	neighborhoodsObject.erase(it);

	ColocationType subColocation;
	for (auto it = neighborhoodsObject.begin(); it != neighborhoodsObject.end(); it++) {
		FeatureType feature = get<FeatureType>(*it);
		subColocation.push_back(feature);
	}

	//查找 B2 C1是不是子集
	if (!_cliqueInstances[k - 1].empty() && !_cliqueInstances[k - 1][subColocation].empty()) {
		auto& subinstance = _cliqueInstances[k - 1][subColocation];
		flag = binary_search(subinstance.begin(), subinstance.end(), neighborhoodsObject);
	}

	return flag;
}

void JoinLess::_filterCliqueInstances(map<ColocationType, TableRowNeighborhoodsType>& colocationInstanceMap, int k) {
	for (auto& colocationInstance : colocationInstanceMap) {
		auto& colocations = colocationInstance.first;
		auto& tableRowNeighborhoods = colocationInstance.second;

		for (auto it = tableRowNeighborhoods.begin(); it != tableRowNeighborhoods.end();) {
			NeighborhoodsObjectType neighborhoodsObjectType = (*it);
			if (!isClique(neighborhoodsObjectType, k)) {
				it = tableRowNeighborhoods.erase(it);
			}
			else {
				++it;
			}
		}
	}

	_cliqueInstances[k] = colocationInstanceMap;
}

void JoinLess::_selectPrevalentColocations(int k) {
	auto& colocationInstanceMap = _cliqueInstances[k];

	for (auto it = colocationInstanceMap.begin(); it != colocationInstanceMap.end();) {
		auto& colocations = (*it).first;
		auto& tableRowNeighborhoods = (*it).second;

		if (!isPrevalentByPartiIndex(colocations, tableRowNeighborhoods)) {
			it = colocationInstanceMap.erase(it);
		}
		else {
			_prevalentColocation[k].push_back(colocations);
			++it;
		}
	}
}

//判断一个集合是另一个集合的子集
bool issubset(ColocationType colocation_sub, ColocationType  colocation) {
	set<FeatureType> sub_set(colocation_sub.begin(), colocation_sub.end());
	set<FeatureType> colocatioin_set(colocation.begin(), colocation.end());
	for (auto& sub_item : sub_set) {
		if (colocatioin_set.find(sub_item) == colocatioin_set.end()) {
			return false;
		}
	}

	return true;
}

vector<unsigned int>  getFeatureIdx(const ColocationType& colocation, const ColocationType& antecedent) {
	vector<unsigned int> featureIdx;

	int pos = 0;
	//A B ;A B C 
	for (unsigned int i = 0; i < colocation.size(); i++) {
		if (colocation[i] == antecedent[pos]) {
			featureIdx.push_back(i);
			pos++;
		}
		if (pos == antecedent.size())  break;
	}

	return featureIdx;
}

unsigned int getProjectNumOfColocation(TableRowNeighborhoodsType tableInstance, vector<unsigned int> featureIdx) {
	set<NeighborhoodsObjectType> rowInstanceProjectSet;

	for (auto rowInstance : tableInstance) {
		NeighborhoodsObjectType rowInstanceObjects;
		//得到投影的模式的行实例个数
		for (unsigned int i = 0; i < featureIdx.size(); i++) {
			rowInstanceObjects.push_back(rowInstance[featureIdx[i]]);
		}
		rowInstanceProjectSet.insert(rowInstanceObjects);
	}

	return rowInstanceProjectSet.size();
}

unsigned int  JoinLess::getRowInstancesOfColocationSub(const ColocationType& colocationSub) {
	int k = colocationSub.size();
	return _cliqueInstances[k][colocationSub].size();
}

void JoinLess::_generateRules() {
	//获取colocationSubSet
	ColocationSetType colocationSubSet;
	ColocationSetType colocationOneSet = _prevalentColocation[1];
	for (auto colocationOne : colocationOneSet) {
		colocationSubSet.push_back(colocationOne);
	}

	int length = _cliqueInstances.size();
	for (unsigned int k = 2; k <= length; k++) {
		map<ColocationType, TableRowNeighborhoodsType> colocationPackages = _cliqueInstances[k];

		for (auto colocationPackage : colocationPackages) {
			ColocationType colocation = colocationPackage.first;
			TableRowNeighborhoodsType tableInstance = colocationPackage.second;

			for (auto colocationSub : colocationSubSet) {
				if (!issubset(colocationSub, colocation)) {
					continue;
				}

				//abc - bc = a, a =>bc  
				ColocationType antecedent;
				set_difference(colocation.begin(), colocation.end(), colocationSub.begin(), colocationSub.end(), back_inserter(antecedent));

				vector<unsigned int> featureIdx = getFeatureIdx(colocation, antecedent);

				//获得分子：abc在ab投影下的行实例数
				unsigned int projectNumOfColocation = getProjectNumOfColocation(tableInstance, featureIdx);

				//分母
				unsigned int rowInstancesOfColocationSub = getRowInstancesOfColocationSub(antecedent);

				//条件概率是  abc中a去重 / bc
				double conf = projectNumOfColocation * 1.0 / rowInstancesOfColocationSub;
				if (conf >= _min_cond_prob) {
					_rules.insert(move(Rule{ antecedent, colocationSub, conf }));
				}
			}

			//这一轮得到的低级放入下一轮 作为 colocationSubSet
			colocationSubSet.push_back(colocation);
		}
	}
}

set<Rule> JoinLess::execute()
{
	_gen_star_neighborhoods();

	int k = 2;
	while (_prevalentColocation.count(k - 1) && !_prevalentColocation[k - 1].empty()) {
		ColocationSetType candidates = _generateCandidateColocations_k(k);
		map<ColocationType, TableRowNeighborhoodsType> SIk = _filterStarInstances(candidates, k);

		if (k == 2) {
			_cliqueInstances[2] = SIk;
		}
		else {
			//根据参与度进行剪枝,SIk woule be changed
			_select_coarse_prevalent_colocations(SIk);
			//A1 B1 C1，判断B1 C1是不是BC下的，不是剪枝
			_filterCliqueInstances(SIk, k);
		}

		//两次剪枝后根据真正的参与度进行剪枝
		_selectPrevalentColocations(k);
		k++;
	}

	_generateRules();

	return _rules;
}

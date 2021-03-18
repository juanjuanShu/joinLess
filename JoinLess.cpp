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

		//_colocationNum (����������)
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

			//��������ǰ�����������ݵ�feature����
			//if (feature1 >= feature2) continue;

			//������Ҳ�Ǵ���ģ���Ϊ�ߺͱ�֮�������ͬ���ϵ����Ϊ���պϳɵ������ĺͱߵĹ�ϵ
			//if (feature1 == feature2) continue;

			if ((new Common())->hasRelation(location1, location2, _star_radius)) {
				// feature1.id1 should less than feature2.id2.����δ����������ΪĬ��������������
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

	//�洢
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

	//����
	for (auto& item : trie_colocationSet) {
		ColocationType candidate = item.first;
		//��������kֻ��һ�������޷�����
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

	//�ռ����к�ѡcolocation�ĵ�һ��Feature 
	set<FeatureType> centerFeatureSet;
	for (auto& candidate : candidates) {
		centerFeatureSet.insert(candidate[0]);
	}

	//��ѡcolocation�ĵ�һ��Featureһ�������͵�����
	for (auto centerFeature : centerFeatureSet) {
		auto& idStarNeighborhoods = _starNeighborhoods[centerFeature];
		for (auto idStarNeighborhood : idStarNeighborhoods) {
			InstanceIdType centerId = idStarNeighborhood.first;
			NeighborhoodsObjectType neighborhoodsObject = idStarNeighborhood.second;

			if (neighborhoodsObject.size() < k) continue;

			//����ÿһ���������ɵı�ʵ��
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
	int remainder,//��Ҫ����λ���м���
	ColocationType& tmp_colocation,
	NeighborhoodsObjectType& tmp_neighborhoods_object,
	const ColocationSetType& candidates
) {
	if (pos + remainder > neighborhoodsObject.size()) return;

	if (0 == remainder) {
		//��һ���ж��ǲ�����candidates���Ӽ�
		if (isSubCandidate(tmp_colocation, candidates)) {
			ans[tmp_colocation].push_back(tmp_neighborhoods_object);
		}
		return;
	}

	//ȡ��Ԫ��
	tmp_neighborhoods_object[k - remainder] = neighborhoodsObject[pos];
	tmp_colocation[k - remainder] = neighborhoodsObject[pos].first;
	_gen_star_k_instances_recursion(ans, neighborhoodsObject, k, pos + 1, remainder - 1, tmp_colocation, tmp_neighborhoods_object, candidates);

	//��ȡ��Ԫ�أ�remainder��λ�ñ�������
	_gen_star_k_instances_recursion(ans, neighborhoodsObject, k, pos + 1, remainder, tmp_colocation, tmp_neighborhoods_object, candidates);
}

void JoinLess::_gen_star_k_instances(map<ColocationType, TableRowNeighborhoodsType>& ans, const ColocationSetType& candidates,
	const NeighborhoodsObjectType& neighborhoodsObject,
	FeatureType centerFeature,
	int k)
{
	NeighborhoodsObjectType newNeighborhoodsObject;

	//Ԥ����ռ�
	ColocationType tmp_colocation(k);
	NeighborhoodsObjectType tmp_neighborhoods_object(k);

	//����Feature�ǵ�һ��Ԫ��
	tmp_colocation[0] = centerFeature;
	tmp_neighborhoods_object[0] = neighborhoodsObject[0];

	_gen_star_k_instances_recursion(ans, neighborhoodsObject, k, 1, k - 1, tmp_colocation, tmp_neighborhoods_object, candidates);

}

bool JoinLess::isPrevalentByPartiIndex(const ColocationType& colocations, const TableRowNeighborhoodsType& tableRowNeighborhoods) {
	bool isPrevalent = true;

	//��ʼ��λͼ
	map<FeatureType, vector<bool>> bitMap;
	for (unsigned int i = 0; i < colocations.size(); i++) {
		FeatureType feature = colocations[i];
		//numOfInstances[feature]��feature��ʵ����
		bitMap[feature] = vector<bool>(_numOfInstances[feature], false);
	}

	//��ֵ 
	// A B ;����A 
	// 0 1
	for (unsigned int i = 0; i < colocations.size(); i++) {
		FeatureType feature = colocations[i];
		for (auto rowInstance : tableRowNeighborhoods) {
			InstanceIdType id = get<InstanceIdType>(rowInstance[i]);
			bitMap[feature][id - 1] = true;
		}
	}

	//����
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
	//���ݲ���Ƚ��м�֦
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

	//���� B2 C1�ǲ����Ӽ�
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

//�ж�һ����������һ�����ϵ��Ӽ�
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
		//�õ�ͶӰ��ģʽ����ʵ������
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
	//��ȡcolocationSubSet
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

				//��÷��ӣ�abc��abͶӰ�µ���ʵ����
				unsigned int projectNumOfColocation = getProjectNumOfColocation(tableInstance, featureIdx);

				//��ĸ
				unsigned int rowInstancesOfColocationSub = getRowInstancesOfColocationSub(antecedent);

				//����������  abc��aȥ�� / bc
				double conf = projectNumOfColocation * 1.0 / rowInstancesOfColocationSub;
				if (conf >= _min_cond_prob) {
					_rules.insert(move(Rule{ antecedent, colocationSub, conf }));
				}
			}

			//��һ�ֵõ��ĵͼ�������һ�� ��Ϊ colocationSubSet
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
			//���ݲ���Ƚ��м�֦,SIk woule be changed
			_select_coarse_prevalent_colocations(SIk);
			//A1 B1 C1���ж�B1 C1�ǲ���BC�µģ����Ǽ�֦
			_filterCliqueInstances(SIk, k);
		}

		//���μ�֦����������Ĳ���Ƚ��м�֦
		_selectPrevalentColocations(k);
		k++;
	}

	_generateRules();

	return _rules;
}

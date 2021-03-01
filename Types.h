#pragma once
#include "stdafx.h"

//ʵ��<ʵ��ID,�ռ��������ͣ�λ��>
using InstanceIdType = unsigned int;
using FeatureType = unsigned char;
using LocationType = pair<double, double>;
using InstanceType = tuple<InstanceIdType, FeatureType, LocationType>;

//��ʵ���ͱ�ʵ��
using RowInstanceType = vector<InstanceIdType>;
using TableInstanceType = vector<RowInstanceType>;

//co-locationģʽ
using ColocationType = vector<FeatureType>;
//��ѡģʽ,ֻ�洢����    vector<vector<FeatureType>>  
using ColocationSetType = vector<ColocationType>;
using ColocationPackage = map<ColocationType,TableInstanceType>;

using InstanceObjectType = pair<FeatureType, InstanceIdType>;
using NeighborhoodsObjectType = vector<InstanceObjectType>;
using TableRowNeighborhoodsType = vector<NeighborhoodsObjectType>;

using RelationPairType = pair<InstanceObjectType, InstanceObjectType>;
using RelationsType = set<RelationPairType>;



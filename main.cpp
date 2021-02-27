#include "stdafx.h"
#include "Types.h"
#include "JoinLess.h"
#include <iostream>
#include <sstream>

void loadData(vector<InstanceType>& instances);

int main() {
    vector<InstanceType> instances;
    loadData(instances);

    JoinLess joinLess(instances,0.2,0,4);
    joinLess.execute();

    return 0;
}

void loadData(vector<InstanceType>& instances) {
    instances.push_back({ 1, 'A', {0, 2} });
    instances.push_back({ 2, 'A', {11, 3} });
    instances.push_back({ 3, 'A', {12,3} });
    instances.push_back({ 4, 'A', {10, 4} });

    instances.push_back({ 1, 'B', {0, 1} });
    instances.push_back({ 2, 'B', {13,1.5} });
    instances.push_back({ 3, 'B', {20, 5} });
    instances.push_back({ 4, 'B', {11,2} });
    instances.push_back({ 5, 'B', {31,11} });

    instances.push_back({ 1, 'C', {12,1.5} });
    instances.push_back({ 2, 'C', {0, 3} });
    instances.push_back({ 3, 'C', {30, 10} });
    
    //A1 B1 {0, 2} {0, 1}    1
    //A1 C2 {0, 2} {0, 3}    1
    //A3 C1 {12,3} {12,1.5} 2.25
    //A2 B4 {11,3} {11,2}    1
    //A3 B4 {12,3} {11,2}    2
    //B2 C1  {13,1.5} {12,1.5} 1
    //B4 C1  {11,2}  {12,1.5}  1.25
    //B5 C3  {31,11} {30, 10}  4

    //ifstream ifs("colocation.csv", ios::in);
    /*ifstream ifs("1.txt", ios::in);

    string line;
    while (getline(ifs, line)) {
        for (int i = 0; i < line.size(); ++i) {
            if (line[i] == ',') {
                line[i] = ' ';
            }
        }

        stringstream ss(line);
        FeatureType feature;  InstanceIdType instanceId; double x, y;
        ss >> feature >> instanceId >> x >> y;
        instances.push_back(make_tuple(instanceId, feature, make_pair(x, y)));
    }*/
}
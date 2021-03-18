#include "stdafx.h"
#include "Types.h"
#include "JoinLess.h"
#include <iostream>
#include <sstream>
#include <chrono>

using chrono::high_resolution_clock;
using chrono::milliseconds;

void loadData(vector<InstanceType>& instances, string inputPath);
void visualization(set<Rule> rules);

int main(int argc, char* argv[]) {
    if (argc != 5) {
        cout << "Argument number must be 4" << endl;
        cout << "./JoinLess :minimum_prevalence ,minimum_rule_probability, maximum_neighbourhood_distance ,inputPath" << endl;
        return 0;
    }
    double min_prev = stod(argv[1]), min_cond_prob = stod(argv[2]), star_radius = stod(argv[3]);
    string inputPath(argv[4]);

    vector<InstanceType> instances;
    loadData(instances, inputPath);

    high_resolution_clock::time_point beginTime = high_resolution_clock::now();

    JoinLess joinLess(instances, min_prev, min_cond_prob, star_radius);
    set<Rule> rules = joinLess.execute();

    high_resolution_clock::time_point endTime = high_resolution_clock::now();
    milliseconds timeInterval = chrono::duration_cast<milliseconds>(endTime - beginTime);
    cout << timeInterval.count() << "ms" << endl;

    visualization(rules);
    return 0;
}

void loadData(vector<InstanceType>& instances, string inputPath) {
    /* instances.push_back({ 1, 'A', {0, 2} });
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
     instances.push_back({ 3, 'C', {30, 10} });*/


    ifstream ifs(inputPath, ios::in);

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
    }
}

void visualization(set<Rule> rules) {
    ColocationType antecedent;
    ColocationType consequent;

    for (auto& rule : rules)
    {
        antecedent = rule.antecedent;
        for (auto& item : antecedent) {
            cout << item;
            int size = antecedent.size() - 1;
            if (item == antecedent[size])
                cout << " ";
            else
                cout << "^";
        }

        cout << " => ";

        consequent = rule.consequent;
        for (auto& item : consequent) {
            cout << item;
            int size = consequent.size() - 1;
            if (item == consequent[size])
                cout << " ";
            else
                cout << "^";
        }

        cout << "   confidence  : " << rule.conf;

        cout << endl;
    }
}
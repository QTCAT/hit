#include <Rcpp.h>
#include <map>
#include <vector>
#include <string>

using namespace Rcpp;
using std::string;
using std::vector;
using std::map;

void findNames(List x, vector<string>& name);
CharacterVector names(List x);
void dendIndex(map<string, int>& dandNameIndex,
               CharacterVector& dendNames, CharacterVector& newNames);
IntegerVector clusterIndex(List x, map<string, int>& dandNameIndex);
void runDend(List x, map< int, IntegerVector>& hierarchy, int& counter,
             int superset, int level, NumericVector& height, 
             map<string, int>& dandNameIndex);
List  dend2hier(List x, NumericVector height, CharacterVector newNames);

// [[Rcpp::export]]
List dend2hier(List x, NumericVector height, CharacterVector newNames)
{
    CharacterVector dendNames = names(x);
    map<string, int> dandNameIndex; 
    dendIndex(dandNameIndex, dendNames, newNames);
    // initial cluster
    IntegerVector cluster = wrap(dandNameIndex);
    cluster.attr("height") = height[0];
    vector<int> nullSubset;
    cluster.attr("subset") = nullSubset;
    // output 
    map< int, IntegerVector> out;
    out[0] = cluster;
    int counter = 1;
    // run through the dendrogram
    runDend(x, out, counter, 0, 0, height, dandNameIndex);
    return wrap(out);
}

void runDend(List x, map< int, IntegerVector>& hierarchy, int& counter,
             int superset, int level, NumericVector& height, 
             map<string, int>& dandNameIndex)
{
    RObject node;
    double nodeHeight = 0;
    int newLevel = 0;
    for (int i = 0; i < x.size(); ++i) {
        node = as<RObject>(x[i]);
        for (int j = 0; j < height.size(); ++j) {
            nodeHeight = node.attr("height");
            if (nodeHeight <= height[j])
                newLevel = j;
            else
                break;
        }
        if (newLevel > level) {
            IntegerVector cluster;
            vector<int>  subset;
            if (node.hasAttribute("leaf") && node.attr("leaf")) {
                // make cluster
                string name = as<string>(node.attr("label"));
                cluster = dandNameIndex[name];
                cluster.attr("height") = height[newLevel];
                cluster.attr("superset") = superset +1;
                hierarchy[counter] = cluster;
                // modify super cluster
                subset = hierarchy[superset].attr("subset");
                subset.push_back(counter +1);
                hierarchy[superset].attr("subset") = subset;
                ++ counter;
            }
            else {
                // make cluster
                cluster = clusterIndex(x[i], dandNameIndex);
                cluster.attr("height") = height[newLevel];
                cluster.attr("superset") = superset +1;
                vector<int> nullSubset;
                cluster.attr("subset") = nullSubset;
                hierarchy[counter] = cluster;
                // modify super cluster
                subset = hierarchy[superset].attr("subset");
                subset.push_back(counter +1);
                hierarchy[superset].attr("subset") = subset;
                int newSuperset = counter;
                ++ counter;
                runDend(x[i], hierarchy, counter, newSuperset, 
                        newLevel, height, dandNameIndex);
            }
        }
        else if (!node.hasAttribute("leaf")) {
            runDend(x[i], hierarchy, counter, superset, 
                    level, height, dandNameIndex);
        }
    }
}

// dendrogram names in order to the new names
void dendIndex(map<string, int>& dandNameIndex,
               CharacterVector& dendNames, CharacterVector& newNames)
{
    for (int i = 0; i < dendNames.size(); ++i) {
        for (int j = 0; j < newNames.size(); ++j) {
            if (dendNames[i] == newNames[j]) {
                dandNameIndex[as<string>(dendNames[i])] = j + 1;
                break;
            }
        }
    }
}

// index of cluster
IntegerVector clusterIndex(List x, map<string, int>& dandNameIndex)
{
    CharacterVector name = names(x);
    int n = name.size();
    IntegerVector out(n);
    for (int i = 0; i < n; ++i) {
       out[i] = dandNameIndex[as<string>(name[i])];
    }
    return out.sort();
}

// names of cluster
// [[Rcpp::export]]
CharacterVector names(List x)
{
    vector<string> out;
    findNames(x, out);
    return wrap(out);
}

void findNames(List x, vector<string>& name) 
{
    for (int i = 0; i < x.size(); ++i) {
        RObject node = as<RObject>(x[i]);
        if (node.hasAttribute("label"))
            name.push_back(as<string>(node.attr("label")));
        else
            findNames(x[i], name);
    }
}

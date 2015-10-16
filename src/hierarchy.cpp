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


//' @title Create a hierarchy from a dendrogram
//' 
//' @description A function which can be called from R. It is creating a 
//' hierarchy from a dendrogram.
//' 
//' @param x A dendrogram S3 R object.
//' @param height A vector of heights at which nodes are grouped.
//' @param newNames Labels of the variabels which should be part of the 
//' hierarchy.
//'
//'@keywords internal
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
    map<int, IntegerVector> out;
    out[0] = cluster;
    int counter = 1;
    // run through the dendrogram
    runDend(x, out, counter, 0, 0, height, dandNameIndex);
    return wrap(out);
}


// @title Creat a node of a hierarchy
// 
// @description A function which recursively is called to generate all nodes 
// in the hierarchy. Only call from within a C++ function!
// 
// @param x A dendrogram S3 R object.
// @param hierarchy A map to which the node is added.
// @param counter A interger for the position of the node in the 
// hierarchy.
// @param superset A integer giving the position of the next higher node.
// @param level A integer of the level of heights.
// @param height A vector of heights at which nodes are grouped.
// @param dandNameIndex Name index to add to node.
// 
// @keywords internal
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
            RObject superNode = as<RObject>(hierarchy[superset]);
            vector<int>  superSubset;
            if (node.hasAttribute("leaf") && node.attr("leaf")) {
                // make cluster
                string name = as<string>(node.attr("label"));
                cluster = dandNameIndex[name];
                cluster.attr("height") = height[newLevel];
                cluster.attr("superset") = superset +1;
                hierarchy[counter] = cluster;
                // modify super cluster
                if (superNode.hasAttribute("subset")) {
                    superSubset = hierarchy[superset].attr("subset");
                    superSubset.push_back(counter +1);
                    hierarchy[superset].attr("subset") = superSubset;
                }
                else
                    hierarchy[superset].attr("subset") = counter +1;
                ++ counter;
            }
            else {
                // make cluster
                cluster = clusterIndex(x[i], dandNameIndex);
                cluster.attr("height") = height[newLevel];
                cluster.attr("superset") = superset +1;
                hierarchy[counter] = cluster;
                // modify super cluster
                if (superNode.hasAttribute("subset")) {
                    superSubset = hierarchy[superset].attr("subset");
                    superSubset.push_back(counter +1);
                    hierarchy[superset].attr("subset") = superSubset;
                }
                else
                    hierarchy[superset].attr("subset") = counter +1;
                int newSuperset = counter;
                ++ counter;
                runDend(x[i], hierarchy, counter, newSuperset, 
                        newLevel, height, dandNameIndex);
            }
        }
        else if (!node.hasAttribute("leaf"))
            runDend(x[i], hierarchy, counter, superset, 
                    level, height, dandNameIndex);
    }
}


// @title dendrogram names in order to the new names
//
// @description A function which creates a nameIndex. Only call from within a 
// C++ function!
//
// @param dandNameIndex A map to which the dendrogram names index is writen.
// @param dendNames The names of variables which are part of the dendrogram.
// @param newNames Labels of the variabels which should be part of the 
// hierarchy.
// 
// @keywords internal
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


// @title Index of cluster
// 
// @description A function which creates a Index of clusters. Only call from 
// within a C++ function!
//
// @param x A dendrogram S3 R object.
// @papam dandNameIndex  A map of name indexes.
//
//@keywords internal
IntegerVector clusterIndex(List x, map<string, int>& dandNameIndex)
{
    CharacterVector name = names(x);
    int n = name.size();
    IntegerVector out(n);
    for (int i = 0; i < n; ++i)
        out[i] = dandNameIndex[as<string>(name[i])];
    return out.sort();
}


// @titile Names of cluster
// 
// @description A function which creates a vector of names
// 
// @param x A dendrogram S3 R object.
//
//@keywords internal
CharacterVector names(List x)
{
    vector<string> out;
    findNames(x, out);
    return wrap(out);
}


// @title Find names of dendrogram
// 
// @description A function which finds names by recursively calling it self.
//
// @param x A dendrogram S3 R object.
// @param name A vector in which the names are writrn.
//
//@keywords internal
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

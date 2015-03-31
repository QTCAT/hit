#include <Rcpp.h>

using namespace Rcpp;



//  make hierarchy from clustering
// [[Rcpp::export]]
List clusterInx(IntegerMatrix x) {
    int n = x.ncol();
    std::vector<bool> chacked(n);
    std::map< int, std::vector<int> > clustMemb;
    std::map< int, std::vector<int> > clustLevel;
    std::map< int, int > clustNo;
    std::vector<int> inx;
    int nInx = 0;
    int a, b;
    bool newClust = false;
    int clustCount = 1;
    int clustNow = 1;
    // run thru all rows
    for (int k = 0; k < x.nrow(); ++ k) {
        // set chack to zero
        for (int i = 0; i < chacked.size (); ++ i) {
            chacked[i] = true;
        }
        // compare element 'i' with all following elements
        for (int i = 0; i < n; ++ i) {
            if (IntegerVector::is_na(x(k, i))) { 
              // if NA nexed
              continue;
            } else {
                // if 'i' is not part already part of a cluster make it part 
                // of a new cluster
                if (chacked[i]) {
                    inx.push_back(i);
                    // if further elements belong to the same cluster add them to
                    for (int j = i +1; j < n; ++ j) {
                        if (x(k, i) == x(k, j)) {
                            inx.push_back(j);
                        }
                    }
                    nInx = inx.size();
                    a = inx[0];
                    // if 'k' is in row 2 or higher, test if similar cluster 
                    // existed at previous level as well
                    if (k > 0) {
                        for (int j = 0; j < n; ++ j) {
                            if (x(k -1, a) == x(k -1, j)) {
                                newClust = true;
                                for (int l = 0; l < nInx; ++ l) {
                                    if (j == inx[l]) {
                                        newClust = false;
                                        break;
                                    }
                                }
                            }
                            if (newClust) {
                                break;
                            }
                        }
                    } else {
                        newClust = true;
                    }
                    // if cluster is new at it to the index tree
                    if (newClust) {
                        for (int j = 0; j < inx.size (); ++ j) {
                            b = inx[j];
                            clustMemb[clustCount].push_back(b +1); // make R index 
                            chacked[b] = false;
                        }
                        clustLevel[k +1].push_back(clustCount); // push to R index
                        clustNo[a] = clustCount;
                        ++ clustCount; 
                    } else {
                        for (int j = 0; j < inx.size (); ++ j) {
                            b = inx[j];
                            chacked[b] = false;
                        }
                        clustLevel[k +1].push_back(clustNo[a]); // push to R index
                    }
                    newClust = false;
                }
                inx.clear();
            }
        }
    }
    return List::create(Named("hierarchyCluster", clustLevel),
                        Named("clusterMembers", clustMemb));
} // clusterInx

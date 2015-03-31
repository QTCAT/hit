#include <Rcpp.h>

using namespace Rcpp;



// find significant cluster from p-value matrix
// [[Rcpp::export]]
NumericMatrix sigCluster(NumericVector pVal, List clustLevel, List clustMemb,
                         const int p, const double alpha, 
                         int minDistInx, int maxDistInx) 
{
    // R Indices to C++
    minDistInx --;
    maxDistInx --;
    // Initialisation
    NumericMatrix sigClust(2, p);
    for (int l = 0; l < p; ++ l) {
        sigClust(1, l) = 1;
    }
    IntegerVector gIDs;
    std::vector<int> gInx;
    int k = 1;
    int m, h;
    bool newClust = true;
    bool isZero = false;
    bool isNotZero = false;
    // start searching between min and max distance
    for (int i = minDistInx; i > maxDistInx; -- i) {
        gIDs = clustLevel[i];
        // test each cluster at dist i
        for (int j = 0; j < gIDs.size(); ++ j) {
            m = gIDs[j] -1; // R Indix to C++
            gInx = clustMemb[m]; 
            // test if all variables of cluster have already been significant 
            // before at smaller distances
            for (int l = 0; l < gInx.size(); ++ l) {
                h = gInx[l] -1; // R Indix to C++
                if (sigClust(0, h) == 0) {
                    isZero = true;
                }
                if (sigClust(0, h) != 0) {
                    isNotZero = true;
                }
                // is variable already part of sig cluster?
                // or p-value is too large jump to next cluster
                if ((isZero && isNotZero) || (pVal[m] > alpha)) {
                    newClust = false;
                    break;
                }
            }
            // if newClust == true make sig cluster
            if (newClust) {
                for(int l = 0; l < gInx.size(); ++ l) {
                    h = gInx[l] -1; // R Indix to C++
                    sigClust(0, h) = k;
                    sigClust(1, h) = pVal[m];
                }
                ++ k;
            }
            newClust = true;
            isZero = false;
            isNotZero = false;
            gInx.clear();  
        }
    }
    // ordert cluster index 
    NumericVector cInx = sort_unique(sigClust(0, _));
    for (int i = 0; i < cInx.size(); ++ i) {
        for (int l = 0; l < p; ++ l) {
            if (sigClust(0, l) == cInx[i]) {
                sigClust(0, l) = i;
            }
        }
    }
    return sigClust;
} // end sigCluster()

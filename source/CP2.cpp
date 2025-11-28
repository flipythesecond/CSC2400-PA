/************************************************************************************************************************************
    Name: Driver.cpp                                                                                                                |
    Compile: g++ Driver.cpp -o runPA1                                                                                               |
    Authors: Bekhruz Anvarov, Ben Nunley, Lance Johnston                                                                            |
    Date: 11/26/2025                                                                                                                   |
    Credit Statement:                                                                                                               |
    Example code such as how to use <ctime> library, etc. was provided by April Crockett during CSC-1300 and CSC-1310 courses.      |
    Pseudo code reference for bubbleSort & mergeSort provided by Prantar Ghosh in lecture notes.                                    |
\***********************************************************************************************************************************/

#include <algorithm>
#include "CP2.h"


double euclideanDist(const City& a, const City& b) {
    double dx = a.xCoord - b.xCoord;
    double dy = a.yCoord - b.yCoord;
    return sqrt(dx * dx + dy * dy);
}

bool compareXCoord(const City& a, const City& b){
    return a.xCoord < b.xCoord;
}

bool compareYCoord(const City& a, const City& b){
    return a.yCoord < b.yCoord;
}

ClosestResult BFClosest(City cities[], int n){
    ClosestResult best;

    // best.dist = numeric_limits<double>::max(); // cmin = inf (from textbook)
    best.dist = euclideanDist(cities[0], cities[1]);
    best.cityID1 = cities[0].cityID;
    best.cityID2 = cities[1].cityID;

    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){

            double itDist = euclideanDist(cities[i], cities[j]);

            if (itDist < best.dist){
                best.dist = itDist;
                best.cityID1 = cities[i].cityID;
                best.cityID2 = cities[j].cityID;
            }
        }
    }
    return best;
}

ClosestResult BFRange(const vector<City> &pts, int left, int right){
    ClosestResult best;

    int count = right - left + 1;

    if (count < 2){
        cout << "[DEBUG]: ERROR count is less than 2..." << endl;
    }

    best.dist = euclideanDist(pts[left], pts[left + 1]);
    best.cityID1 = pts[left].cityID;
    best.cityID2 = pts[left + 1].cityID;


    //Check pairs
    for(int i = left; i <= right; i++){
        for (int j = i + 1; j <= right; j++){
            
            double itDist = euclideanDist(pts[i], pts[j]);
            
            if (itDist < best.dist) {
                best.dist = itDist;
                best.cityID1 = pts[i].cityID;
                best.cityID2 = pts[j].cityID;
            }
        }
    }
    return best;
}

ClosestResult closestUtil(vector<City> &ptsX, int left, int right){
    int n = right - left + 1;

    if (n <= 3){
        return BFRange(ptsX, left, right);
    }

    int mid = left + (right - left) / 2;
    double midX = ptsX[mid].xCoord;

    ClosestResult leftRes = closestUtil(ptsX, left, mid);
    ClosestResult rightRes = closestUtil(ptsX, mid + 1, right);
    ClosestResult best;

    if (leftRes.dist < rightRes.dist){
        best = leftRes;
    }
    else{
        best = rightRes;
    }

    double currDist = best.dist;
    vector<City> strip;
    strip.reserve(n);
    
    for(int i = left; i <= right; i++){
        if(fabs(ptsX[i].xCoord - midX) < currDist){
            strip.push_back(ptsX[i]);
        }
    }

    sort(strip.begin(), strip.end(), compareYCoord);

    int sSize = strip.size();

    for(int i = 0; i < sSize; i++){
        for(int j = i + 1; j < sSize; j++){

            if (strip[j].yCoord - strip[i].yCoord >= currDist){
                break;
            }

            double itDist = euclideanDist(strip[i], strip[j]);

            if (itDist < best.dist){
                best.dist = itDist;
                best.cityID1 = strip[i].cityID;
                best.cityID2 = strip[j].cityID;
                currDist = itDist;
            }
        }
    }
    return best;
}

// TODO: divdedAndConquer implementation here
ClosestResult divideAndConquer(City cities[], int n){

    vector<City> ptsX;
    ptsX.reserve(n);

    for(int i = 0; i < n; i++){
        ptsX.push_back(cities[i]);
    }

    sort(ptsX.begin(), ptsX.end(), compareXCoord);

    return closestUtil(ptsX, 0, n - 1);

}


/*-----------Function Declerations for reference--------------------*\
bool compareYCoord(const City& a, const City& b);
bool compareXCoord(const City& a, const City& b);
ClosestResult BFClosest(City cities[], int n);
ClosestResult BFRange(const vector<City> &pts, int left, int right);
ClosestResult closestUtil(vector<City> &ptsX, int left, int right);
ClosestResult divideAndConquer(City cities[], int n);
***********************************************************************/

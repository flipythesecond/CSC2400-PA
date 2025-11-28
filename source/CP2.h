/************************************************************************************************************************************
    Name: CP2.h                                                                                                                     |                                                                                            |
    Authors: Bekhruz Anvarov, Ben Nunley, Lance Johnston                                                                            |
    Date: 11/26/2025                                                                                                                |
                                                                                                                                    |
\***********************************************************************************************************************************/

#ifndef CP2_H
#define CP2_H

// Library includes
#include <iostream>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
using namespace std;

// City Struct
struct City {
    int cityID;
    double xCoord;
    double yCoord;
};

struct ClosestResult {
    double dist;
    int cityID1;
    int cityID2;
};


// Function Declarations
void runFlights();
void runClosestPair();
//
void bubbleSort(double FlightTimeHour[], double FlightCost[], int size);
void mergeSort(double arr[], int left, int right);
void merge(double arr[], int left, int mid, int right);
//
double euclideanDist(const City &a, const City &b);
bool compareYCoord(const City& a, const City& b);
bool compareXCoord(const City& a, const City& b);
ClosestResult BFClosest(City cities[], int n);
ClosestResult BFRange(const vector<City> &pts, int left, int right);
ClosestResult closestUtil(vector<City> &ptsX, int left, int right);
ClosestResult divideAndConquer(City cities[], int n);



//
#endif

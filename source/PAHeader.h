/************************************************************************************************************************************
    Name: PAHeader.h                                                                                                                |                                                                                            |
    Authors: Bek Anvarov, Ben Nunley, Lance Johnston                                                                                |
    Date: 11/26/2025                                                                                                                |
    Purpose:                                                                                                                        |
      Header file for Project Assignment (PA) Driver.cpp, it declares:                                                              |
        Data structures for representing cities and closest-pair results.                                                           |
        Driver functions for running the flights (Checkpoint 1) and closest pair (Checkpoint 2) algorithims.                        |                                                                                                             |
        Sorting functions (Bubble Sort and Merge Sort) used for flight data.                                                        |
        Brute force and divide and conquer functions for the closest pair of cities problem.                                        |                                                                                                                                                       |
\***********************************************************************************************************************************/

#ifndef PAHEADER_H
#define PAHEADER_H

// Library includes
#include <iostream>
#include <iomanip>
#include <ctime>        // clock(), CLOCK_PER_SEC for runtime mesurements
#include <cmath>        // sqrt, fabs, pow for distance calculations
#include <fstream>      // ifstream/ofstream for reading/writing files
#include <string>
#include <sstream>      // string stream for parsing
#include <vector>       // vector used in D&C closest pair algorithim
#include <algorithm>
using namespace std;

// City Struct
//Represents a single city with an ID and (x,y) coordiantes, cities.txt parsed (i.e) line 1: 1 2XX.XX 3XX.XX
struct City {
    int cityID;    
    double xCoord; 
    double yCoord;
};

//ClosestResult struct
// Stores result of a closest pair computation of the distance and the two city IDs
struct ClosestResult {
    double dist;
    int cityID1;
    int cityID2;
};


//---------------------- Function Declarations -----------------------------//

//Sectioned off the checkpoints for organization purposes
void runBubbleAndMerge(); // Checkpoint 1 function call, uses flights.txt 
void runClosestPair();   // Checkpoint 2 function call, uses cities.txt

// Checkpoint 1 algorithims and helper functions
void bubbleSort(double FlightTimeHour[], double FlightCost[], int size);
void mergeSort(double arr[], int left, int right);
void merge(double arr[], int left, int mid, int right);

// Checkpoint 2 algorithims and helper functions
double euclideanDist(const City &a, const City &b); // Computes Euclidean distance between two cities using their (x, y) coordinates
bool compareYCoord(const City& a, const City& b);   // Sort cities by y-coordinates
bool compareXCoord(const City& a, const City& b);   // Sort cities by x-coordinates
ClosestResult BFClosest(City cities[], int n);      // Brute-force closest pair over the first n cities in a plain array
ClosestResult BFRange(const vector<City> &pts, int left, int right); // Brute-force closest pair over a subrange [left, right] of a vector
ClosestResult closestUtil(vector<City> &ptsX, int left, int right);  // Recursive helper for divide-and-conquer closest pair
ClosestResult divideAndConquer(City cities[], int n); // // Copies cities[] into a vector, sorts by x, and calls closestUtil



//
#endif

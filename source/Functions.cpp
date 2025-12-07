/************************************************************************************************************************************
    Name: Functions.cpp                                                                                                             |
    Authors: Bek Anvarov, Ben Nunley, Lance Johnston                                                                            |
    Date: 11/26/2025                                                                                                                |
    Purpose:                                                                                                                        |
      Implements the core algorithms used in PA1:                                                                                   |
        Bubble Sort and Merge Sort for processing flight times and costs               (Checkpoint 1).                              |
        Euclidean distance computation for 2D city coordinates.                        (Checkpoint 2)                               |
        Brute force closest-pair algorithm over an array or subrange of cities.        (Checkpoint 2)                               |
        Divide and conquer closest pair algorithm using sorting by (x, y) coordinates. (Checkpoint 2)                               |
\***********************************************************************************************************************************/

#include "PAHeader.h"


//--- Bubble Sort modified from April Crockett and Prantar Ghosh's Pseudo Code ---//
void bubbleSort(double FlightTimeHour[], double FlightCost[], int size){
    
    for(int i = 0; i < size - 1; i++){
        for(int j = 0; j < size - i - 1; j++){
            if(FlightTimeHour[j] > FlightTimeHour[j + 1]){
                double tempT = FlightTimeHour[j];
                FlightTimeHour[j] = FlightTimeHour[j + 1];
                FlightTimeHour[j + 1] = tempT;
                
                //cout << "Time Swapped " << FlightTimeHour[j + 1] << " and " << FlightTimeHour[j] << endl;
                //break;
            }
            
            if(FlightCost[j] > FlightCost[j + 1]){
                double tempC = FlightCost[j];
                FlightCost[j] = FlightCost[j + 1];
                FlightCost[j + 1] = tempC;
                
                //cout << "Cost Swapped " << FlightCost[j + 1] << " and " << FlightCost[j] << endl;
            }
        }
    }

    // Open and write Bubble Sort output to file and close
    ofstream outputFile1("output/FtimeBubSort.txt", ios::app);
    ofstream outputFile2("output/FcostbubSort.txt", ios::app);
    for(int i = 0; i < size; i++){
        outputFile1 << FlightTimeHour[i] << " ";
    }
    outputFile1 << endl;
    for(int i = 0; i < size; i++){
        outputFile2 << FlightCost[i] << " ";
    }
    outputFile2 << endl;
    outputFile1.close();
    outputFile2.close();

}

//--- MergeSort modified from April Crockett ---//
void mergeSort(double arr[], int left, int right){

    //Base case
    if (left >= right ){
        return;
    }
    
    // Initlize mid point
    int mid = left + (right - left) / 2;    
    
    // Recursively sort first and second halves
    mergeSort(arr, left, mid);
    mergeSort(arr, mid + 1, right);
    merge(arr, left, mid, right);

}

//--- Merge algorithim modified from April Crockett ---//
void merge(double arr[], int left, int mid, int right){

    //Initilize variables
    int mergedSize = right - left + 1;
    int mergePos = 0;
    int leftPos = left;
    int rightPos = mid + 1;

    // Create temporary arrays
    double* temp = new double[mergedSize];

    // Merge the two subarrays into temp[]
    while (leftPos <= mid && rightPos <= right) {
        if(arr[leftPos] <= arr[rightPos]) {
            temp[mergePos] = arr[leftPos];
            leftPos++;
            mergePos++;
         
        } 
        else {
            temp[mergePos] = arr[rightPos];
            rightPos++;
            mergePos++;
        }
    }

    // Copy the remaining elements of left subarray
    while (leftPos <= mid) {
        temp[mergePos] = arr[leftPos];
        leftPos++;
        mergePos++;
    }
    
    // Copy the remaining elements of right subarray
    while (rightPos <= right) {
        temp[mergePos] = arr[rightPos];
        rightPos++;
        mergePos++;
    }

    // Copy the merged elements back into original array
    for (int i = 0; i < mergedSize; i++) {
        arr[left + i] = temp[i];
    }

    // Delete temporary array
    delete[] temp;
}


// Computes Euclidean distance between two cities using their (x, y) coordinates
double euclideanDist(const City& a, const City& b) {
    double dx = a.xCoord - b.xCoord;
    double dy = a.yCoord - b.yCoord;
    return sqrt(dx * dx + dy * dy);
}

// Comparator for sorting cities by x coordinate (used in divide and conquer)
bool compareXCoord(const City& a, const City& b){
    return a.xCoord < b.xCoord;
}

//// Comparator for sorting cities by y coordinate (used when processing the strip)
bool compareYCoord(const City& a, const City& b){
    return a.yCoord < b.yCoord;
}

//Brute force closest pair over the first n cities in the array
//Checks all O(n^2) pairs and returns the smallest distance and the city IDs
ClosestResult BFClosest(City cities[], int n){
    ClosestResult best;

    // Initilized best distance between the first two cities, either that or numeric limits initilized to the highest value
    // best.dist = numeric_limits<double>::max(); // cmin = inf (from textbook)
    best.dist = euclideanDist(cities[0], cities[1]); 
    best.cityID1 = cities[0].cityID;
    best.cityID2 = cities[1].cityID;

    // Try every pair (i, j) with i < j
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){

            double itDist = euclideanDist(cities[i], cities[j]);

            // If this pair is closer, update the best result
            if (itDist < best.dist){
                best.dist = itDist;
                best.cityID1 = cities[i].cityID;
                best.cityID2 = cities[j].cityID;
            }
        }
    }
    return best;
}

// Brute force closest pair in a subrange [left, right] of a vector of cities
// Used as a helper for the divide-and-conquer algorithm when n is small
ClosestResult BFRange(const vector<City> &pts, int left, int right){
    ClosestResult best;

    int count = right - left + 1;

    // Should not happen in normal use; just a debug message
    if (count < 2){
        cout << "[DEBUG]: ERROR count is less than 2..." << endl;
    }

    // Initialize best using the first two cities in the range
    best.dist = euclideanDist(pts[left], pts[left + 1]);
    best.cityID1 = pts[left].cityID;
    best.cityID2 = pts[left + 1].cityID;


    // Check all pairs in [left, right]
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

// Recursive divide and conquer helper for closest pair
// Assumes ptsX is sorted by x coordinate and works on subarray [left, right]
ClosestResult closestUtil(vector<City> &ptsX, int left, int right){
    int n = right - left + 1;

    // Base case: small number of points → just brute force them
    if (n <= 3){
        return BFRange(ptsX, left, right);
    }

    // Find middle index and its x coordinate
    int mid = left + (right - left) / 2;
    double midX = ptsX[mid].xCoord;

    // Recursively find closest pairs in left and right halves
    ClosestResult leftRes = closestUtil(ptsX, left, mid);
    ClosestResult rightRes = closestUtil(ptsX, mid + 1, right);
    ClosestResult best;
    // Choose the closer of the two results as the current best
    if (leftRes.dist < rightRes.dist){
        best = leftRes;
    }
    else{
        best = rightRes;
    }

    double currDist = best.dist;

    // Build the strip: points within currDist of the vertical line x = midX
    vector<City> strip;
    strip.reserve(n);
    
    for(int i = left; i <= right; i++){
        if(fabs(ptsX[i].xCoord - midX) < currDist){
            strip.push_back(ptsX[i]);
        }
    }

    // Sort the strip by y coordinate to apply the strip scanning step
    sort(strip.begin(), strip.end(), compareYCoord);

    int sSize = strip.size();

    // Check each point in the strip against the next few points by y
    for(int i = 0; i < sSize; i++){
        for(int j = i + 1; j < sSize; j++){

            // If y distance is already >= currDist, no need to check further
            if (strip[j].yCoord - strip[i].yCoord >= currDist){
                break;
            }

            double itDist = euclideanDist(strip[i], strip[j]);

            // Update best if we find a closer pair that crosses the midline
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

// Divide and conquer driver for closest pair
// Copies the array of cities into a vector, sorts by x, and calls closestUtil
// Time complexity should be O(n log n)
ClosestResult divideAndConquer(City cities[], int n){

    vector<City> ptsX;
    ptsX.reserve(n);

    // Copy cities into vector so we can sort safely
    for(int i = 0; i < n; i++){
        ptsX.push_back(cities[i]);
    }

    // Sort by x-coordinate once at the start
    sort(ptsX.begin(), ptsX.end(), compareXCoord);

    // Recursive call to compute the closest pair
    return closestUtil(ptsX, 0, n - 1);

}
//
int knapMax(const vector<int> &weights, int W){
    
    int n = static_cast<int>(weights.size());
    int wi, vi, skip, take;

    vector<vector<int>> A(n + 1, vector<int>(W + 1, 0));

    for(int i = 1; i <= n; i++){

        wi = weights[i - 1];
        vi = 1;
        
        for(int C = 1; C <= W; C++){
            if(wi > C){

                A[i][C] = A[i - 1][C];
            }
            else{

                skip = A[i-1][C];
                take = A[i-1][C - wi] + vi;
                A[i][C] = max(skip, take);
            }


        }

    }
    return A[n][W];
}

/*-----------Function Declerations for reference--------------------*\
bool compareYCoord(const City& a, const City& b);
bool compareXCoord(const City& a, const City& b);
ClosestResult BFClosest(City cities[], int n);
ClosestResult BFRange(const vector<City> &pts, int left, int right);
ClosestResult closestUtil(vector<City> &ptsX, int left, int right);
ClosestResult divideAndConquer(City cities[], int n);
***********************************************************************/

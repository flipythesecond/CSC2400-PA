/************************************************************************************************************************************
    Name: Driver.cpp                                                                                                                |
    Compile: g++ Driver.cpp -o runProg2                                                                        |
    Authors: Bek Anvarov, Ben Nunley, Lance Johnston                                                                                |
    Date: 11/7/20                                                                                                                   |
    Purpose:                                                                                                                        |
        This Program runs the functions for the Project Assignment for each checkpoint                                              |
                                                                                                                                    |
        Checkpoint 1 - Reads flight data, sorts it using Bubble Sort and Merge Sort, and                                            |
        records the sorted results and runtimes to output files                                                                     |
                                                                                                                                    |
        Checkpoint 2 - Reads city coordinate data, computes the closest pair of cities using                                        |
        both a brute-force algorithm and a divide and conquer algorithm, and records the                                            |
        closest pairs and runtimes for comparison.                                                                                  |
                                                                                                                                    |
    Credit Statement:                                                                                                               |
        Example code such as how to use <ctime> library, etc. was provided by April Crockett during CSC-1300 and CSC-1310 courses.  |
        Pseudo code reference for bubbleSort & mergeSort provided by Prantar Ghosh in lecture notes.                                |
\***********************************************************************************************************************************/

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

//----------- Main Function -----------//
int main() {
    
    // Program Start Message
    cout << "Starting Program..." << endl;

    /***************************************\
    | Uncomment to run either function,     |
    | haven't ran both at the same time yet |
    \***************************************/

    //runBubbleAndMerge(); // (Checkpoint 1)


    runClosestPair();   // (Checkpoint 2)
    
    return 0;
}

//------------- Checkpoint 1 Flights Function ----------------//
void runBubbleAndMerge(){
    cout << "\n[DEBUG]: Running flights bubble/merge sort..." << endl;
    // Variable Declarations
    string lineTest;
    fstream file("flights.txt");
    
    /* check if file opened successfully
    if (!file.is_open()) {
        cerr << "ERROR: Could not open flights.txt" << endl;
        return 1;
    }
    cout << "Successfully opened flights.txt" << endl;
    */

  

    clock_t tStartBub, tStopBub;
    clock_t tStartMer, tStopMer;
    double compute_time_Bub = 0, compute_time_Mer = 0;
    

    
    /*vector<double> FlightTimeHour;
    vector<double> FlightCost;
    vector<int> CityNum; */

    
    double FlightTimeHour[100];
    double FlightCost[100];
    int CityNum[100];
    
    // truncate files at start of program (prevents duplicates)
    ofstream outputFile1("FtimeBubSort.txt", ios::trunc);
    ofstream outputFile2("FcostbubSort.txt", ios::trunc);
    ofstream outputFile3("FtimeMerSort.txt", ios::trunc);
    ofstream outputFile4("FcostMerSort.txt", ios::trunc);
    ofstream runTimeFile("runtimes.txt", ios::trunc);
    runTimeFile.close();
    
    //ofstream runTimeFile("runtimes.txt", ios::trunc);
    // initlize n and delim for parsing
    int n;
    char delim;
    
    // While loop to read each line from input file
    while(getline(file, lineTest)) {
        if(lineTest.empty()) {
            continue;
        }
     
        n = 0;
        
        stringstream ss(lineTest);
        
        while(ss >> delim) {
             
                if(delim == '('){
                        ss >> CityNum[n]  >> delim >> FlightTimeHour[n] >> delim >> FlightCost[n] >> delim;
                        //cout << "City: " << CityNum[n]  << " flightTimeHour: " << FlightTimeHour[n] << " flightCost: " << FlightCost[n] << endl;
                        n++;
                }
        }

        // Copying original arrays for merge sort
        double mergeTime[100];
        double mergeCost[100];
        for(int i = 0; i < n; i++){ // For loop to copy original arrays
            mergeTime[i] = FlightTimeHour[i];
            mergeCost[i] = FlightCost[i];
        }
        
        /*Debugging Original
        cout << "\n[DEBUG] Orignial Time Array: ";
        for(int i = 0; i < n; i++){
            cout << FlightTimeHour[i] << " ";
        }
        cout << endl;

        cout << "\n[DEBUG] Orignial Cost Array: ";
        for(int i = 0; i < n; i++){
            cout << FlightCost[i] << " ";
        }
        cout << endl;
        */
        // Bubble Sort Timing
        tStartBub = clock();
        bubbleSort(FlightTimeHour, FlightCost, n);
        tStopBub = clock();
        compute_time_Bub += (double)(tStopBub - tStartBub) / CLOCKS_PER_SEC;
        
        /*Debugging Bubble
        cout << "\n[DEBUG] Bubble Time Array: ";
        for(int i = 0; i < n; i++){
            cout << FlightTimeHour[i] << " ";
        }
        cout << endl;

        cout << "\n[DEBUG] Bubble Cost Array: ";
        for(int i = 0; i < n; i++){
            cout << FlightCost[i] << " ";
        }
        */

        // Merge Sort Timing
        tStartMer = clock();
        mergeSort(mergeTime, 0 , n - 1);
        
        // Write merge sorted output to file
        for(int i = 0; i < n; i++){
            outputFile3 << mergeTime[i] << " ";

        }
        outputFile3 << endl;
        
        //
        mergeSort(mergeCost, 0 , n - 1);

        for(int i = 0; i < n; i++){
            outputFile4 << mergeCost[i] << " ";

        }
        outputFile4 << endl;
        tStopMer = clock();
        compute_time_Mer += (double)(tStopMer - tStartMer) / CLOCKS_PER_SEC;
        
        /* Debugging Merge
        cout << "\n[DEBUG] Merge Time Array: ";
        for(int i = 0; i < n; i++){
            cout << FlightTimeHour[i] << " ";
        }
        cout << endl;

        cout << "\n[DEBUG] Merge Cost Array: ";
        for(int i = 0; i < n; i++){
            cout << FlightCost[i] << " ";
        }
        */

        //Runtime output and conversion
        ofstream runtimeFile("runtimes.txt", ios::app);
        //Conversion to nanoseconds
        long long BubTime = compute_time_Bub*1000000; //  Long Long for removing smaller 
        long long MerTime = compute_time_Mer*1000000; //  values i.e. (3.714.....e-14^10)
        
        
        runtimeFile << "(" << BubTime << ", " << MerTime << ")" << endl;   
        runtimeFile.close(); 
        
    }
    
    outputFile3.close();
    outputFile4.close();
    file.close();


    cout << "[DEBUG]: Flights function finished running..." <<  endl;

    

}

//---------- Checkpoint 2 Brute Force & Divide and Conquer --------------//
void runClosestPair(){
    

    cout << "[DEBUG]: Running closest pair BF & DC..." << endl;

    ifstream inFile("cities.txt");

    if(!inFile.is_open()){
        cout << "[DEBUG]: (ERROR) Couldn't open cities.txt" << endl;
        return;
    }

    clock_t tStartBF, tStopBF;
    clock_t tStartDC, tStopDC;
    double compute_time_BF = 0, compute_time_DC = 0;
    

    if(!inFile.is_open()){
        cout << "[DEBUG]: (ERROR) Coudln't open cities.txt" << endl;
    }

    City cities[100];
    int numCities = 0;
    string lineTest;

   while(getline(inFile, lineTest)){
        if(lineTest.empty()){
            continue;
        }

        stringstream ss(lineTest);

        ss >> cities[numCities].cityID
           >> cities[numCities].xCoord
           >> cities[numCities].yCoord;

        numCities++;

        if(numCities >= 100){
            cout << "[DEBUG]: (WARNING!) reached 100 cities..." << endl;
            break;
        }
   }
    inFile.close();

    ofstream bfOut("BF-Closest.txt", ios::trunc);
    ofstream dcOut("DC-Closest.txt", ios::trunc);
    ofstream runtime("runtimes.txt", ios::trunc);

    
   if(!bfOut.is_open() || !dcOut.is_open() || !runtime.is_open()){
    cout << "[DEBUG]: (ERROR) Coudln't open one of the city output files..." << endl;
    return;
   }

    for(int i = 50; i <= numCities; i++){

        tStartBF = clock();
        ClosestResult bfRes = BFClosest(cities, i);
        tStopBF = clock();
        compute_time_BF += (double)(tStopBF - tStartBF) / CLOCKS_PER_SEC;

        tStartDC = clock();
        ClosestResult dcRes = divideAndConquer(cities, i);
        tStopDC = clock();
        compute_time_DC = (double)(tStopDC - tStartDC) / CLOCKS_PER_SEC;


        bfOut << bfRes.cityID1 << " " << bfRes.cityID2 << " " << bfRes.dist << endl;
        dcOut << dcRes.cityID1 << " " << dcRes.cityID2 << " " << dcRes.dist << endl;

        long long BFTime = compute_time_BF*1000000; //  Long Long for removing smaller 
        long long DCTime = compute_time_DC*1000000; //  values i.e. (3.714.....e-14^10)
        
        runtime << "(" << BFTime << ", " << DCTime << ")" << endl;
    }

    bfOut.close();
    dcOut.close();
    runtime.close();

    cout << "[DEBUG]: Closest pair funcation finished running..." << endl;

}


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
    ofstream outputFile1("FtimeBubSort.txt", ios::app);
    ofstream outputFile2("FcostbubSort.txt", ios::app);
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


// // Computes Euclidean distance between two cities using their (x, y) coordinates
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

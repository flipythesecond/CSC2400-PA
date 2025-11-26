/*
    Name: Driver.cpp
    Compile: g++ Driver.cpp -o runPA1
    Authors: Bekhruz Anvarov, Ben Nunley, Lance Johnston
    Date: 11/7/20
    Credit Statement: 
    Example code such as how to use <ctime> library, etc. was provided by April Crockett during CSC-1300 and CSC-1310 courses.
    Pseudo code reference for bubbleSort & mergeSort provided by Prantar Ghosh in lecture notes.
*/

// Library includes
#include <iomanip>
#include <iostream>
#include <ctime>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
using namespace std;

// Function Prototypes
void bubbleSort(double FlightTimeHour[], double FlightCost[], int size);
void mergeSort(double arr[], int left, int right);
void merge(double arr[], int left, int mid, int right);

int main() {
    
    // Program Start Message
    cout << "Starting Program..." << endl;


    // Variable Declarations
    string lineTest;
    fstream file("given/flights.txt");
    
    /* check if file opened successfully
    if (!file.is_open()) {
        cerr << "ERROR: Could not open ../given/flights.txt" << endl;
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
    ofstream outputFile1("output/FtimeBubSort.txt", ios::trunc);
    ofstream outputFile2("output/FcostbubSort.txt", ios::trunc);
    ofstream outputFile3("output/FtimeMerSort.txt", ios::trunc);
    ofstream outputFile4("output/FcostMerSort.txt", ios::trunc);
    ofstream runTimeFile("output/runtimes.txt", ios::trunc);
    
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
        ofstream runtimeFile("output/runtimes.txt", ios::app);
        //Conversion to nanoseconds
        long long BubTime = compute_time_Bub*1000000; //  Long Long for removing smaller 
        long long MerTime = compute_time_Mer*1000000; //  values i.e. (3.714.....e-14^10)
        
        runtimeFile << "(" << BubTime << ", " << MerTime << ")" << endl;   
        runtimeFile.close(); 
        
    }
    
    outputFile3.close();
    outputFile4.close();
    return 0;
}

// Bubble Sort modified from April Crockett and Prantar Ghosh's Pseudo Code
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

//MergeSort modified from April Crockett 
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

// Merge algorithim modified from April Crockett
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

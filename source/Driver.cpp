/************************************************************************************************************************************
    Name: Driver.cpp                                                                                                                |
    Compile: g++ source/Driver.cpp source/Functions.cpp                                                                             |
    Authors: Bek Anvarov, Ben Nunley, Lance Johnston                                                                                |
    Date: 11/7/20                                                                                                                   |
    Purpose:                                                                                                                        |
        This Program runs the functions for the Project Assignment for each checkpoint                                              |
                                                                                                                                    |
        Checkpoint 1 – Reads flight data, sorts it using Bubble Sort and Merge Sort, and                                            |
        records the sorted results and runtimes to output files                                                                     |
                                                                                                                                    |
        Checkpoint 2 – Reads city coordinate data, computes the closest pair of cities using                                        |
        both a brute-force algorithm and a divide and conquer algorithm, and records the                                            |
        closest pairs and runtimes for comparison.                                                                                  |
                                                                                                                                    |
    Credit Statement:                                                                                                               |
        Example code such as how to use <ctime> library, etc. was provided by April Crockett during CSC-1300 and CSC-1310 courses.  |
        Pseudo code reference for bubbleSort & mergeSort provided by Prantar Ghosh in lecture notes.                                |
\***********************************************************************************************************************************/

// Library includes
#include "PAHeader.h" // Contains Functions and Libs


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
    runTimeFile.close();
    
    //ofstream runTimeFile("output/runtimes.txt", ios::trunc);
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
    file.close();


    cout << "[DEBUG]: Flights function finished running..." <<  endl;

    

}

//---------- Checkpoint 2 Brute Force & Divide and Conquer --------------//
void runClosestPair(){
    cout << "[DEBUG]: Running closest pair BF & DC..." << endl;

    ifstream inFile("given/cities.txt");

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

        /*
        if(numCities >= 100){
            cout << "[DEBUG]: (WARNING!) reached 100 cities..." << endl;
            break;
        }
        */
   }
    inFile.close();

    ofstream bfOut("output/BF-Closest.txt", ios::trunc);
    ofstream dcOut("output/DC-Closest.txt", ios::trunc);
    ofstream runtime("output/runtimes.txt", ios::trunc);

    
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
        compute_time_DC += (double)(tStopDC - tStartDC) / CLOCKS_PER_SEC;


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

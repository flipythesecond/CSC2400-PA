# CSC 2400 – Programming Assignment 1

## Overview



This program implements two main parts for PA1:

**BEFORE COMPILING UNCOMMENT FUNCTION CALLS IN DRIVER.cpp**

**To Compile using GCC/g++ run:** g++ source/Driver.cpp source/Function.cpp
1. **Checkpoint 1 – Flight Sorting**
   - Reads flight data from `given/flights.txt`.
   - Uses **Bubble Sort** and **Merge Sort** to sort:
     - Flight times (hours)
     - Flight costs
   - Writes the sorted results and runtimes to:
     - `output/FtimeBubSort.txt`
     - `output/FcostBubSort.txt`
     - `output/FtimeMerSort.txt`
     - `output/FcostMerSort.txt`
     - `output/runtimes.txt` (Bubble vs Merge runtime per line)

2. **Checkpoint 2 – Closest Pair of Cities**
   - Reads city coordinates from `given/cities.txt`.
   - Represents each city with an ID, x-coordinate, and y-coordinate.
   - Computes the closest pair of cities using:
     - A **brute-force** algorithm `BFClosest`
     - A **divide-and-conquer** algorithm `divideAndConquer`
   - For increasing values of *n* (from 50 up to the total number of cities), it:
     - Finds the closest pair using each algorithm.
     - Writes results to:
       - `output/BF-Closest.txt` (brute force)
       - `output/DC-Closest.txt` (divide & conquer)
       - `output/runtimes.txt` (BF vs DC runtime per line)
      
3. **Checkpoint 3 – Knapsack-Based Round-Trip Optimization**
   * This checkpoint uses **dynamic programming (0/1 Knapsack)** to maximize the number of cities that can be visited under a fixed **$5000 round-trip budget**, based on round-trip ticket cost data.

   ### Overview
   - Reads round-trip ticket data from `given/roundtrip_costs.txt`, where each line represents a starting city and contains `(city, cost)` pairs.
   - Parses each line to extract ticket costs and ignores any tickets that individually exceed the budget.
   - Converts ticket costs from `double` to `int` so they can be used as knapsack weights.

   ### Knapsack Model
   * Each round-trip ticket is treated as a knapsack item:
      - **Weight** = round-trip cost  
      - **Value** = 1 (one city visited)
   
   * The objective is to select a subset of tickets that **maximizes the number of cities visited** without exceeding the budget.
   
   ### Dynamic Programming Algorithm
   * The function `knapMax(weights, W)` builds a DP table:
---

## File Structure


**[NOTE]**: Both runtimes are outputted to runtimes.txt
        either call runBubbleAndMerge() or runClosestPair()

**[LOG]**: Log files are git ignored, these are outputs
       generated from the program itself
```text
.
├── source/
│   ├── Driver.cpp       # main() – calls runBubbleAndMerge() / runClosestPair() / runRoundTrip()
│   ├── Functions.cpp    # algorithms for PA1 (sorting + closest pair + DP Knapsack)
│   └── PAHeader.h       # structs (City, ClosestResult) + function prototypes
├── given/
│   ├── flights.txt         # input flight data
│   ├── cities.txt          # input city coordinates
|   └── roundtrip_costs.txt # input round trip costs
├── output/
│   ├── FtimeBubSort.txt  # [LOG]
│   ├── FcostBubSort.txt  # [LOG]
│   ├── FtimeMerSort.txt  # [LOG]
│   ├── FcostMerSort.txt  # [LOG]
│   ├── BF-Closest.txt    # [LOG]
│   ├── DC-Closest.txt    # [LOG]
│   ├── trip_nums.txt     # [LOG] Knapsack Algorithim
│   ├── runtimes.txt      # [LOG] Bubble vs Merge   
│   └── runtimes.txt      # [LOG] BF vs DC        
├── .gitignore
├── .gitattributes
└── README.md

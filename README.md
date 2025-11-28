# CSC 2400 – Programming Assignment 1

## Overview



This program implements two main parts for PA1:

**BEFORE COMPILING UNCOMMENT FUNCTION CALLS IN DRIVER.cpp**
**To Compile run:** g++ source/Driver.cpp source/Function.cpp
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

---

## File Structure


**[NOTE]**: Both runtimes are outputted to runtimes.txt
        either call runBubbleAndMerge() or runClosestPair()

**[LOG]**: Log files are git ignored, these are outputs
       generated from the program itself
```text
.
├── source/
│   ├── Driver.cpp       # main() – calls runBubbleAndMerge() / runClosestPair()
│   ├── Functions.cpp    # algorithms for PA1 (sorting + closest pair)
│   └── PA1.h            # structs (City, ClosestResult) + function prototypes
├── given/
│   ├── flights.txt      # input flight data
│   └── cities.txt       # input city coordinates
├── output/
│   ├── FtimeBubSort.txt  # [LOG]
│   ├── FcostBubSort.txt  # [LOG]
│   ├── FtimeMerSort.txt  # [LOG]
│   ├── FcostMerSort.txt  # [LOG]
│   ├── BFClosest.txt     # [LOG]
│   ├── DCClosest.txt     # [LOG]
│   ├── runtimes.txt      # [LOG] Bubble vs Merge   
│   └── runtimes.txt      # [LOG] BF vs DC        
├── .gitignore
├── .gitattributes
└── README.md

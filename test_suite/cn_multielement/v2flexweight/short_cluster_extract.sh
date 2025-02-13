#!/bin/bash

g++ -O3 -o extract_clusters accelerated_extract_clusters.cpp -pthread

./extract_clusters setup.in
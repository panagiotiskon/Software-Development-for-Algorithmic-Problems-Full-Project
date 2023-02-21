## Getting Started:

This project was implemented for the course Developing Software for Algorithmic Problems (K23), during the winter semester 2022-2023 CS NKUA. 

The completion of the project came through three main phases. 

Initially, the Incremental and Convex Hull algorithms were implemented using a given set of points. Next, to improve the effectiveness of 
the previously developed algorithms Local Search and Simulated Annealing(Locally - Globally) were developed. Lastly, a preprocessing stage was created where different parameters of the algorithms were already chosen based on the size of the set and the time required for each algorithm to run. 

The program performs all different algorithm combinations to a folder containing test sets and prints the results in PROJECT_output.txt for further observation.

## Prerequisites:

 - CGAL v. 5.0.2-3 and above.

## How to run:

1. cmake -DCGAL_DIR=$CMAKE_INSTALLED_PREFIX/lib/CGAL -DCMAKE_BUILD_TYPE=Release .

2. make

3. ./PROJECT_3 -i data/uniform/ -o PROJECT_output.txt


## More about the Project:

[project_description.pdf](https://github.com/panagiotiskon/Software-Development-for-Algorithmic-Problems-Full-Project/files/10794330/project_description.pdf)

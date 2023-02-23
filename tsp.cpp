#include <iostream>
#include <fstream>
#include <cstring>

using namespace std;

void readInputFile(string inputFile) {
    ifstream myFile;
    myFile.open(inputFile, ios::in);
    if(!myFile) {
        cout << "NO SOLUTION";
        exit(1);
    }
    int width, height, coord1, coord2;
    double dist;
    myFile >> width >> height;
    int matrix[width][height];

    while (myFile >> coord1 >> coord2 >> dist) {
        matrix[coord1][coord2] = dist;
    }
    myFile.close();
}

int main(int argc, char *argv[]) {
    readInputFile("ex2.in");
    return 0;
}
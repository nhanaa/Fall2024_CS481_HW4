/*
  Name: Pax Nguyen
  Email: ntnguyen11@crimson.ua.edu
  Course Section: CS 481
  Homework #: 4
  Instructions to compile the program: mpic++ -g -Wall -o hw4 hw4.cpp
  Instructions to run the program: mpiexec -n <mumber_of_processes> ./hw4 <board_size> <number_of_generations> <number_of_processes> <output_directory>
*/

#include <vector>
#include <mpi.h>
#include <random>
#include <iostream>
#include <chrono>
#include <fstream>

using namespace std;

inline int idx(int row, int col, int width) {
  return row * width + col;
}

void checkForError(bool localOk, string fname, string message, MPI_Comm comm) {
  if (!localOk) {
    int myRank;
    MPI_Comm_rank(comm, &myRank);
    if (myRank == 0) {
      fprintf(stderr, "Proc %d > In %s, %s\n", myRank, fname.c_str(), message.c_str());
      fflush(stderr);
    }
    MPI_Finalize();
    exit(-1);
  }
}

void initLocal(vector<int>& localArray, vector<int>& counts, vector<int>& displs,
              int boardSize, int myRank, MPI_Comm comm) {
  int localOk = 1;
  string fname = "initLocal";
  vector<int> board;
  int localSize = counts[myRank];

  if (myRank == 0) { // Root process, so create the board
    board = vector<int>(boardSize * boardSize);
    if (board.size() == 0) {
      localOk = 0;
    }
    checkForError(localOk, fname, "Can't allocate temporary vector", comm);

    // Init the game of life board

    for (int i = 0; i < boardSize * boardSize; i++) {
      srand(123|i);
      board[i] = rand() % 2;
    }

    #ifdef DEBUG
      cout << "Process " << myRank << " > Initial board:" << endl;
      for (int i = 0; i < boardSize; i++) {
        for (int j = 0; j < boardSize; j++) {
          cout << board[idx(i, j, boardSize)] << " ";
        }
        cout << endl;
      }
    #endif
  }

  // Scatter the board
  MPI_Scatterv(board.data(), counts.data(), displs.data(), MPI_INT,
               localArray.data(), localSize, MPI_INT, 0, comm);
}

void gatherBoard(vector<int>& gatheredBoard, vector<int>& localArray, int localSize, int boardSize, int myRank, vector<int>& counts, vector<int>& displs, MPI_Comm comm) {
  MPI_Gatherv(localArray.data(), localSize, MPI_INT, gatheredBoard.data(), counts.data(), displs.data(), MPI_INT, 0, comm);
}

int countNeighbors(vector<vector<int>> &board, int i, int j) {
  int aliveCount = 0;

  for (int x = i - 1; x <= i + 1; x++) {
    for (int y = j - 1; y <= j + 1; y++) {
      if (x == i && y == j) {
        continue;
      }

      if (board[x][y] == 1) {
        aliveCount++;
      }
    }
  }

  return aliveCount;
}

void writeToFile(const vector<int>& board, int boardSize,
                string outputDir, int numGenerations, int numProcesses, int myRank) {
  string filename = outputDir + "/output" + to_string(boardSize) + "." +
                    to_string(numGenerations) + "." + to_string(numProcesses) + "." +
                    to_string(myRank);

  ofstream outputFile(filename);
  if (!outputFile.is_open()) {
    cout << "Error opening file " << filename << endl;
    return;
  }

  for (int i = 0; i < boardSize; i++) {
    for (int j = 0; j < boardSize; j++) {
      outputFile << board[idx(i, j, boardSize)] << " ";
    }
    outputFile << endl;
  }

  outputFile.close();
}

void localGameOfLife(vector<int>& localArray, int boardSize, vector<int>& counts, vector<int>& displs, int myRank, int commSize, int numGenerations, MPI_Comm comm) {
  int localSize = counts[myRank];
  int localCols = boardSize;
  int localRows = localSize / boardSize;
  int upperRank, lowerRank;

  MPI_Status status;

  // Create 2D array from 1D vector with 2 ghost rows and columns
  vector<vector<int>> localPortion(localRows + 2, vector<int>(localCols + 2));
  vector<vector<int>> newLocalPortion(localRows + 2, vector<int>(localCols + 2));

  // Map data from 1D to 2D array
  for (int i = 0; i < localRows + 2; i++) {
    for (int j = 0; j < localCols + 2; j++) {
      if (i == 0 || j == 0 || i == localRows + 1 || j == localCols + 1) { // Set ghost cells to 0
        localPortion[i][j] = 0;
        newLocalPortion[i][j] = 0;
      } else {
        localPortion[i][j] = localArray[idx(i - 1, j - 1, boardSize)];
        newLocalPortion[i][j] = 0;
      }
    }
  }

  upperRank = myRank + 1 >= commSize ? MPI_PROC_NULL : myRank + 1;
  lowerRank = myRank - 1 < 0 ? MPI_PROC_NULL : myRank - 1;

  // Now process the game of life for given number of generations
  for (int gen = 0; gen < numGenerations; gen++) {
    bool localChanged = false;

    // Exchange ghost rows with neighbors
    if (myRank % 2 == 0) {
      // Exchange up
      MPI_Sendrecv(&localPortion[localRows][0], localCols + 2, MPI_INT, upperRank, 0,
                   &localPortion[localRows + 1][0], localCols + 2, MPI_INT, upperRank, 0, comm, &status);
    } else {
      // Exchange down
      MPI_Sendrecv(&localPortion[1][0], localCols + 2, MPI_INT, lowerRank, 0,
                   &localPortion[0][0], localCols + 2, MPI_INT, lowerRank, 0, comm, &status);
    }

    // Do the second set of exchanges
    if (myRank % 2 == 1) {
      // Exchange up
      MPI_Sendrecv(&localPortion[localRows][0], localCols + 2, MPI_INT, upperRank, 1,
                   &localPortion[localRows + 1][0], localCols + 2, MPI_INT, upperRank, 1, comm, &status);
    } else {
      // Exchange down
      MPI_Sendrecv(&localPortion[1][0], localCols + 2, MPI_INT, lowerRank, 1,
                   &localPortion[0][0], localCols + 2, MPI_INT, lowerRank, 1, comm, &status);
    }

    // Compute the new local portion
    for (int i = 1; i <= localRows; i++) {
      for (int j = 1; j <= localCols; j++) {
        int aliveCount = countNeighbors(localPortion, i, j);

        if (localPortion[i][j] == 1) {
          newLocalPortion[i][j] = (aliveCount == 2 || aliveCount == 3) ? 1 : 0;
        } else {
          newLocalPortion[i][j] = (aliveCount == 3) ? 1 : 0;
        }

        if (newLocalPortion[i][j] != localPortion[i][j]) {
          localChanged = true;
        }
      }
    }

    // Each MPI process sends its rank to reduction, root MPI process collects the result
    bool globalChanged = false;
    MPI_Allreduce(&localChanged, &globalChanged, 1, MPI_C_BOOL, MPI_LOR, comm);
    if (!globalChanged) {
      cout << "The sum of all flag is " << globalChanged << " after k=" << gen << endl;
      break;
    }

    // Put a barrier here to make sure all processes have computed the new local portion
    MPI_Barrier(comm);

    // Swap the new local portion with the old local portion
    swap(localPortion, newLocalPortion);

    #ifdef DEBUG
      vector<int> gatheredBoard(boardSize * boardSize);
      for (int i = 1; i <= localRows; i++) {
        for (int j = 1; j <= localCols; j++) {
          localArray[idx(i - 1, j - 1, boardSize)] = localPortion[i][j];
        }
      }
      gatherBoard(gatheredBoard, localArray, localSize, boardSize, myRank, counts, displs, comm);
      if (myRank == 0) {
        cout << "Process " << myRank << " > Board at generation " << gen + 1 << ":" << endl;
        for (int i = 0; i < boardSize; i++) {
          for (int j = 0; j < boardSize; j++) {
            cout << gatheredBoard[idx(i, j, boardSize)] << " ";
          }
          cout << endl;
        }
      }
    #endif
  }

  // Copy the new local portion back to the local array
  for (int i = 1; i <= localRows; i++) {
    for (int j = 1; j <= localCols; j++) {
      localArray[idx(i - 1, j - 1, boardSize)] = localPortion[i][j];
    }
  }
}

int main(int argc, char **argv) {
  if(argc != 5) {
    cout << "Usage: ./prog <board_size> <num_generations> <num_threads> <output_dir>" << endl;
    return 1;
  }

  int boardSize = atoi(argv[1]);
  int numGenerations = atoi(argv[2]);
  int numProcesses = atoi(argv[3]);
  string outputDir = argv[4];

  int myRank, commSize;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &commSize);

  vector<int> counts(commSize);
  vector<int> displs(commSize);
  int remain = boardSize % commSize;

  // Compute counts and displacements
  for (int i = 0; i < commSize; i++) {
    counts[i] = boardSize / commSize + ((i < remain) ? 1 : 0);
    counts[i] = counts[i] * boardSize;
  }

  displs[0] = 0;
  for (int i = 1; i < commSize; i++) {
    displs[i] = displs[i - 1] + counts[i - 1];
  }
  int localBoardSize = counts[myRank];

  #ifdef DEBUG
    cout << "Process " << myRank << " > boardSize = " << boardSize << ", localBoardSize = " << localBoardSize << endl;
  #endif

  vector<int> localArray(counts[myRank]);
  initLocal(localArray, counts, displs, boardSize, myRank, MPI_COMM_WORLD);
  auto start = chrono::high_resolution_clock::now();

  // Run game of life
  localGameOfLife(localArray, boardSize, counts, displs, myRank, commSize, numGenerations, MPI_COMM_WORLD);
  vector<int> finalBoard;
  if (myRank == 0) {
    finalBoard = vector<int>(boardSize * boardSize);
  }

  // Gather the final board
  gatherBoard(finalBoard, localArray, localBoardSize, boardSize, myRank, counts, displs, MPI_COMM_WORLD);
  auto end = chrono::high_resolution_clock::now();

  if (myRank == 0) {
    cout << "Time taken: " << chrono::duration<double>(end - start).count() << " seconds" << endl;
  }

  if (myRank == 0) {
    #ifdef DEBUG
      cout << "Process " << myRank << " > Final board:" << endl;
      for (int i = 0; i < boardSize; i++) {
        for (int j = 0; j < boardSize; j++) {
          cout << finalBoard[idx(i, j, boardSize)] << " ";
        }
        cout << endl;
      }
    #endif
    writeToFile(finalBoard, boardSize, outputDir, numGenerations, numProcesses, myRank);
  }

  MPI_Finalize();
  return 0;
}

#include <poly/common.h>
#include <poly/image.h>

void inClassExample() {
    Eigen::VectorXd a;
    Eigen::VectorXd b;
    Eigen::VectorXd v;

    a.resize(8);
    b.resize(8);
    v.resize(7);

    a(0) = 4.0;
    a(1) = 3.0;
    a(2) = 4.0;
    a(3) = 3.0;
    a(4) = 5.0;
    a(5) = 4.0;
    a(6) = 3.0;
    a(7) = 2.0;

    b(0) = 3.0;
    b(1) = 6.0;
    b(2) = 0.0;
    b(3) = 0.0;
    b(4) = 0.0;
    b(5) = 0.0;
    b(6) = 1.0;
    b(7) = 2.0;

    v(0) = a(1) - a(0);
    v(1) = a(2) - a(1);
    v(2) = a(3) - a(2);
    v(3) = a(4) - a(3);
    v(4) = a(5) - a(4);
    v(5) = a(6) - a(5);
    v(6) = a(7) - a(6);

    std::cout << "A: " << endl << a << std::endl;

    Eigen::MatrixXd A;
    A.resize(4, 4);
    A.setZero();

    Eigen::VectorXd B;
    B.resize(4);

    A(0, 0) = 4;
    A(1, 0) = -2;

    A(0, 1) = -2;
    A(1, 1) = 4;
    A(2, 1) = -2;

    A(1, 2) = -2;
    A(2, 2) = 4;
    A(3, 2) = -2;

    A(2, 3) = -2;
    A(3, 3) = 4;

    cout << "Q: " << endl << A << endl;

    B(0) = -(2.0 * (-b(1) - v(1) + v(2)));
    B(1) = -(2.0 * (-v(2) + v(3)));
    B(2) = -(2.0 * (-v(3) + v(4)));
    B(3) = -(2.0 * (-b(6) - v(4) + v(5)));

    cout << "B: " << endl << B << endl;

    Eigen::BiCGSTAB<Eigen::MatrixXd > solver;
    solver.compute(A);
    Eigen::VectorXd x = solver.solve(B);
    std::cout << "#iterations:     " << solver.iterations() << std::endl;
    std::cout << "estimated error: " << solver.error()      << std::endl;

    cout << "X: " << endl << x << endl;
    // /* ... update b ... */
    // x = solver.solve(b); // solve again
}

double smartAccess(int row, int col, const Eigen::MatrixXd& map) {
    if (row >= map.rows() || row < 0) return 0.0;
    if (col >= map.cols() || col < 0) return 0.0;
    return map(row, col);
}

Eigen::MatrixXd calculateBoundaryMap(const Eigen::MatrixXd& map) {
    Eigen::MatrixXd boundary;
    boundary.resize(map.rows(), map.cols());
    boundary.setZero();

    for (int i = 0; i < map.rows(); i++) {
        for (int j = 0; j < map.cols(); j++) {
            if (map(i, j) == 0.0) {
                if (smartAccess(i+1, j, map) == 1.0) boundary(i, j) = 1.0;
                if (smartAccess(i-1, j, map) == 1.0) boundary(i, j) = 1.0;
                if (smartAccess(i, j+1, map) == 1.0) boundary(i, j) = 1.0;
                if (smartAccess(i, j-1, map) == 1.0) boundary(i, j) = 1.0;
            }
        }
    }

    return boundary;
}

// Eigen::SparseMatrix<double> calculateBoundaryMap(const& Eigen::SparseMatrix<double> map) {
//     Eigen::SparseMatrix<double> boundary;
//     // TODO
//     return boundary;
// }

// Eigen::MatrixXd calculateGradiantMatrix(const Eigen::MatrixXd& map) {
//     Eigen::MatrixXd grad;
//     grad.resize(map.rows() - 1, map.cols() - 1);
//     grad.setZero();
//
//     for (int i = 0; i < grad.rows(); i++) {
//         for (int j = 0; j < grad.cols(); j++) {
//             grad(i, j) = map()
//         }
//     }
// }

double calculateSumOfGradiant(int i, int j, const Eigen::MatrixXd& oldImg, const Eigen::MatrixXd& newImg, const Eigen::MatrixXd& boundary, const Eigen::MatrixXd& omega) {
    double sum = 0.0;

    // cout << "OMEGA: " << endl << omega << endl;

    if (smartAccess(i, j, omega) != 1.0) cout << "MAJOR ERROR: WHAT" << endl;
    // cout << "POS: " << newImg(i, j) << endl;

    // North Box
    if (smartAccess(i - 1, j, boundary) == 1.0) {
        // cout << "N: " << oldImg(i - 1, j) << endl;
        // cout << "N OLD" << endl;
        sum -= 2.0 * (newImg(i, j) - oldImg(i - 1, j));
    } else if (smartAccess(i - 1, j, omega) == 1.0) {
        // cout << "N: " << newImg(i - 1, j) << endl;
        // cout << "N NEW" << endl;
        sum -= 2.0 * (newImg(i, j) - newImg(i - 1, j));
    } else {
        cout << "MAJOR ERROR: OUT OF EXPECTED BOUNDS NORTH" << endl;
    }

    // cout << "SUM POST N: " << sum << endl;

    // South Box
    if (smartAccess(i + 1, j, boundary) == 1.0) {
        // cout << "S: " << oldImg(i + 1, j) << endl;
        // cout << "S OLD" << endl;
        sum += 2.0 * (oldImg(i + 1, j) - newImg(i, j));
    } else if (smartAccess(i + 1, j, omega) == 1.0) {
        // cout << "S: " << newImg(i + 1, j) << endl;
        // cout << "S NEW" << endl;
        sum += 2.0 * (newImg(i + 1, j) - newImg(i, j));
    } else {
        cout << "MAJOR ERROR: OUT OF EXPECTED BOUNDS SOUTH" << endl;
    }

    // cout << "SUM POST S: " << sum << endl;

    // East Box
    if (smartAccess(i, j + 1, boundary) == 1.0) {
        // cout << "E: " << oldImg(i, j + 1) << endl;
        // cout << "E OLD" << endl;
        sum += 2.0 * (oldImg(i, j + 1) - newImg(i, j));
    } else if (smartAccess(i, j + 1, omega) == 1.0) {
        // cout << "E: " << newImg(i, j + 1) << endl;
        // cout << "E NEW" << endl;
        sum += 2.0 * (newImg(i, j + 1) - newImg(i, j));
    } else {
        // cout << "BOUNDARY: " << smartAccess(i, j + 1, boundary) << endl;
        // cout << "OMEGA: " << smartAccess(i, j + 1, omega) << endl;
        cout << "MAJOR ERROR: OUT OF EXPECTED BOUNDS EAST" << endl;
    }

    // cout << "SUM POST E: " << sum << endl;

    // West Box
    if (smartAccess(i, j - 1, boundary) == 1.0) {
        // cout << "W: " << oldImg(i, j - 1) << endl;
        // cout << "W OLD" << endl;
        sum -= 2.0 * (newImg(i, j) - oldImg(i, j - 1));
    } else if (smartAccess(i, j - 1, omega) == 1.0) {
        // cout << "W: " << newImg(i, j - 1) << endl;
        // cout << "W NEW" << endl;
        sum -= 2.0 * (newImg(i, j) - newImg(i, j - 1));
    } else {
        cout << "MAJOR ERROR: OUT OF EXPECTED BOUNDS WEST" << endl;
    }

    // cout << "SUM POST W: " << sum << endl;

    // cout << "POST: " << endl << omega << endl;

    return sum;
}

// double calculateSumOfGradiant(int i, int j, const Eigen::SparseMatrix<double>& map) {
//     // TODO
//     return 0.0;
// }

int countVars(const Eigen::MatrixXd& omega) {
    int count = 0;

    for (int i = 0; i < omega.rows(); i++) {
        for (int j = 0; j < omega.cols(); j++) {
            if (omega(i, j) > 0.5) count++;
        }
    }

    return count;
}

// int countVars(const Eigen::SparseMatrix<double>& omega) {
//     int count = 0;
//
//     for (int i = 0; i < omega.rows(); i++) {
//         for (int j = 0; j < omega.cols(); j++) {
//             if (omega(i, j) > 0.5) count++;
//         }
//     }
//
//     return count;
// }

// MAJOR SLOW DOWN
int findIndex(const Eigen::MatrixXi& indexMap, int row, int col) {
    for (int i = 0; i < indexMap.rows(); i++) {
        if (indexMap(i, 0) == row && indexMap(i, 1) == col) return i;
    }

    cout << "ERROR: INDEX NOT FOUND" << endl;

    return -1;
}

Eigen::VectorXd solveMissingPoints(const Eigen::MatrixXd& oldImg, const Eigen::MatrixXd& newImg, const Eigen::MatrixXd& boundary, const Eigen::MatrixXd& omega) {
    int vars = countVars(omega);

    cout << endl;

    cout << "NUMBER OF VARS: " << vars << endl;

    Eigen::VectorXd x;
    x.resize(vars);
    x.setZero();

    Eigen::VectorXd b;
    b.resize(vars);
    b.setZero();

    Eigen::MatrixXi indexMap;
    indexMap.resize(vars, 2);
    indexMap.setZero();

    int varIndex = 0;
    for (int i = 0; i < omega.rows(); i++) {
        for (int j = 0; j < omega.cols(); j++) {
            if (omega(i, j) > 0.5) {
                indexMap(varIndex, 0) = i;
                indexMap(varIndex, 1) = j;
                x(i) = newImg(i, j);
                varIndex++;
            }
        }
    }

    cout << "INDEX MAP: " << endl << indexMap << endl;

    Eigen::MatrixXd A;
    A.resize(vars, vars);
    A.setZero();

    for (int i = 0; i < indexMap.rows(); i++) {
        int row = indexMap(i, 0);
        int col = indexMap(i, 1);

        A(i, i) = 8.0;

        if (omega(row, col + 1) > 0.5) {
            int index = findIndex(indexMap, row, col + 1);
            A(i, index) += -2.0;
        } else {
            if (boundary(row, col + 1) > 0.5) {
                b(i) -= 2.0 * oldImg(row, col + 1);
            } else {
                cout << "MAJOR ERROR: NOT OMEGA OR BOUNDARY" << endl;
            }
        }

        if (omega(row, col - 1) > 0.5) {
            int index = findIndex(indexMap, row, col - 1);
            A(i, index) += -2.0;
        } else {
            if (boundary(row, col - 1) > 0.5) {
                b(i) -= 2.0 * oldImg(row, col - 1);
            } else {
                cout << "MAJOR ERROR: NOT OMEGA OR BOUNDARY" << endl;
            }
        }

        if (omega(row - 1, col) > 0.5) {
            int index = findIndex(indexMap, row - 1, col);
            A(i, index) += -2.0;
        } else {
            if (boundary(row - 1, col) > 0.5) {
                b(i) -= 2.0 * oldImg(row - 1, col);
            } else {
                cout << "MAJOR ERROR: NOT OMEGA OR BOUNDARY" << endl;
            }
        }

        if (omega(row + 1, col) > 0.5) {
            int index = findIndex(indexMap, row + 1, col);
            A(i, index) += -2.0;
        } else {
            if (boundary(row + 1, col) > 0.5) {
                b(i) -= 2.0 * oldImg(row + 1, col);
            } else {
                cout << "MAJOR ERROR: NOT OMEGA OR BOUNDARY" << endl;
            }
        }

        b(i) += calculateSumOfGradiant(row, col, oldImg, newImg, boundary, omega);
    }

    cout << "A: " << endl << A << endl;

    cout << "b: " << endl << b << endl;
    cout << "-b: " << endl << -b << endl;


    Eigen::BiCGSTAB<Eigen::MatrixXd > solver;
    solver.compute(A);
    x = solver.solve(-b);
    std::cout << "#iterations:     " << solver.iterations() << std::endl;
    std::cout << "estimated error: " << solver.error()      << std::endl;

    cout << "X: " << endl << x << endl;

    // double test = calculateSumOfGradiant(2, 3, oldImg, newImg, boundary, omega);

    // cout << "OMEGA:" << endl << omega << endl;
    // cout << "BOUNDARY:" << endl << boundary << endl;
    //
    // cout << "GRAD (2, 3): " << test << endl;


    // varIndex = 0;
    // for (int i = 0; i < omega.rows(); i++) {
    //     for (int j = 0; j < omega.cols(); j++) {
    //         if (omega(i, j) > 0.5) {
    //             A(varIndex, varIndex) = 8.0;
    //             if (omega())
    //         }
    //     }
    // }

    // TODO

    return x;
}

void twodExample() {
    Eigen::MatrixXd imageA;
    Eigen::MatrixXd imageB;
    Eigen::MatrixXd v;
    Eigen::MatrixXd replaceMap;
    Eigen::MatrixXd boundaryMap;

    imageA.resize(7, 7);
    imageB.resize(7, 7);
    v.resize(6, 6);
    replaceMap.resize(7, 7);
    boundaryMap.resize(7, 7);

    imageA.setZero();
    imageB.setZero();
    v.setZero();
    replaceMap.setZero();
    boundaryMap.setZero();

    // imageA:
    // 5 4 1 3 1 5 3
    // 8 9 6 7 3 2 4
    // 3 5 1 3 1 5 6
    // 2 3 1 5 1 7 0
    // 3 4 2 4 5 6 9
    // 1 2 4 2 6 3 7
    // 9 8 5 7 4 2 3

    imageA << 0.5, 0.4, 0.1, 0.3, 0.1, 0.5, 0.3,
              0.8, 0.9, 0.6, 0.7, 0.3, 0.2, 0.4,
              0.3, 0.5, 0.1, 0.3, 0.1, 0.5, 0.6,
              0.2, 0.3, 0.1, 0.5, 0.1, 0.7, 0.0,
              0.3, 0.4, 0.2, 0.4, 0.5, 0.6, 0.9,
              0.1, 0.2, 0.4, 0.2, 0.6, 0.3, 0.7,
              0.9, 0.8, 0.5, 0.7, 0.4, 0.2, 0.3;

    // imageB:
    // 6 3 2 2 4 4 1
    // 4 6 3 1 3 5 1
    // 3 4 7 2 5 8 2
    // 7 7 5 3 3 3 2
    // 8 2 8 0 6 2 9
    // 1 4 3 8 3 3 4
    // 2 5 9 3 2 1 1

    imageB << 0.6, 0.3, 0.2, 0.2, 0.4, 0.4, 0.1,
              0.4, 0.6, 0.3, 0.1, 0.3, 0.5, 0.1,
              0.3, 0.4, 0.7, 0.2, 0.5, 0.8, 0.2,
              0.7, 0.7, 0.5, 0.3, 0.3, 0.3, 0.2,
              0.8, 0.2, 0.8, 0.0, 0.6, 0.2, 0.9,
              0.1, 0.4, 0.3, 0.8, 0.3, 0.3, 0.4,
              0.2, 0.5, 0.9, 0.3, 0.2, 0.1, 0.1;

    // replace map:
    // 0 0 0 0 0 0 0
    // 0 0 0 0 0 0 0
    // 0 0 1 1 1 0 0
    // 0 0 1 1 1 0 0
    // 0 0 1 1 1 0 0
    // 0 0 0 0 0 0 0
    // 0 0 0 0 0 0 0

    replaceMap << 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0,
                  0, 0, 1, 1, 1, 0, 0,
                  0, 0, 1, 1, 1, 0, 0,
                  0, 0, 1, 1, 1, 0, 0,
                  0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0;

    // cout << "REPLACE MAP: " << endl << replaceMap << endl;

    // boundary map:
    // 0 0 0 0 0 0 0
    // 0 0 1 1 1 0 0
    // 0 1 0 0 0 1 0
    // 0 1 0 0 0 1 0
    // 0 1 0 0 0 1 0
    // 0 0 1 1 1 0 0
    // 0 0 0 0 0 0 0

    boundaryMap = calculateBoundaryMap(replaceMap);

    cout << endl;

    cout << "Boundary: " << endl << boundaryMap << endl;

    cout << "REPLACE MAP: " << endl << replaceMap << endl;

    cout << "Gradiant sum at (2,2): " << calculateSumOfGradiant(2, 2, imageB, imageA, boundaryMap, replaceMap) << endl;

    // cout << "OMEGA (2, 3): " << replaceMap(0, 1) << endl;

    // cout << endl << replaceMap << endl;

    cout << "BOUND:" << endl;

    cout << boundaryMap << endl;

    cout << endl;

    cout << "REPL:" << endl;

    cout << replaceMap << endl;

    cout << endl;

    Eigen::VectorXd points = solveMissingPoints(imageB, imageA, boundaryMap, replaceMap);
}

void countVars(const Image& omega, int& red, int& green, int& blue) {
    red = 0;
    green = 0;
    blue = 0;

    for (int i = 0; i < omega.rows(); i++) {
        for (int j = 0; j < omega.cols(); j++) {
            if (omega.smartAccess(j, i, 0) > 0.5) red++;
            if (omega.smartAccess(j, i, 1) > 0.5) green++;
            if (omega.smartAccess(j, i, 2) > 0.5) blue++;
        }
    }
}

Image calculateBoundaryMap(const Image& omega) {
    Image boundary(omega.cols(), omega.rows(), 3);
    // boundary.resize(map.rows(), map.cols());
    boundary.setZero();

    for (int c = 0; c < 3; c++) {
        for (int i = 0; i < omega.rows(); i++) {
            for (int j = 0; j < omega.cols(); j++) {
                if (omega.smartAccess(j, i, c) < 0.1) {
                    if (omega.smartAccess(j, i+1, c) > 0.5) boundary.set(j, i, c, 1.0);
                    if (omega.smartAccess(j, i-1, c) > 0.5) boundary.set(j, i, c, 1.0);
                    if (omega.smartAccess(j+1, i, c) > 0.5) boundary.set(j, i, c, 1.0);
                    if (omega.smartAccess(j-1, i, c) > 0.5) boundary.set(j, i, c, 1.0);
                }
            }
        }
    }

    return boundary;
}

void calculateIndexMaps(const Image& omega, Eigen::MatrixXi& redIndexMap, Eigen::MatrixXi& greenIndexMap, Eigen::MatrixXi& blueIndexMap) {
    redIndexMap.setZero();
    greenIndexMap.setZero();
    blueIndexMap.setZero();

    // int varIndex = 0;
    // for (int i = 0; i < omega.rows(); i++) {
    //     for (int j = 0; j < omega.cols(); j++) {
    //         if (omega(i, j) > 0.5) {
    //             indexMap(varIndex, 0) = i;
    //             indexMap(varIndex, 1) = j;
    //             x(i) = newImg(i, j);
    //             varIndex++;
    //         }
    //     }
    // }

    int redVarIndex = 0;
    int greenVarIndex = 0;
    int blueVarIndex = 0;

    for (int i = 0; i < omega.rows(); i++) {
        for (int j = 0; j < omega.cols(); j++) {
            if (omega.smartAccess(j, i, 0) > 0.5) {
                redIndexMap(redVarIndex, 0) = i;
                redIndexMap(redVarIndex, 1) = j;
                redVarIndex++;
            }
            if (omega.smartAccess(j, i, 1) > 0.5) {
                greenIndexMap(greenVarIndex, 0) = i;
                greenIndexMap(greenVarIndex, 1) = j;
                greenVarIndex++;
            }
            if (omega.smartAccess(j, i, 2) > 0.5) {
                blueIndexMap(blueVarIndex, 0) = i;
                blueIndexMap(blueVarIndex, 1) = j;
                blueVarIndex++;
            }
        }
    }
}

void calculateIndexMapsFast(const Image& omega, Eigen::MatrixXi& redIndexMap, Eigen::MatrixXi& greenIndexMap, Eigen::MatrixXi& blueIndexMap) {
    redIndexMap.setZero();
    greenIndexMap.setZero();
    blueIndexMap.setZero();

    // int varIndex = 0;
    // for (int i = 0; i < omega.rows(); i++) {
    //     for (int j = 0; j < omega.cols(); j++) {
    //         if (omega(i, j) > 0.5) {
    //             indexMap(varIndex, 0) = i;
    //             indexMap(varIndex, 1) = j;
    //             x(i) = newImg(i, j);
    //             varIndex++;
    //         }
    //     }
    // }

    int redVarIndex = 0;
    int greenVarIndex = 0;
    int blueVarIndex = 0;

    for (int i = 0; i < omega.rows(); i++) {
        for (int j = 0; j < omega.cols(); j++) {
            if (omega.smartAccess(j, i, 0) > 0.5) {
                redIndexMap(i, j) = redVarIndex++;
            }
            if (omega.smartAccess(j, i, 1) > 0.5) {
                greenIndexMap(i, j) = greenVarIndex++;
            }
            if (omega.smartAccess(j, i, 2) > 0.5) {
                blueIndexMap(i, j) = blueVarIndex++;
            }
        }
    }
}

double calculateSumOfGradiant(int i, int j, const Image& imageB, const Image& imageA, const Image& boundary, const Image& omega, int ch) {
    double sum = 0.0;

    // if (smartAccess(i, j, omega) != 1.0) cout << "MAJOR ERROR: WHAT" << endl;

    if (omega.smartAccess(j, i, ch) != 1.0) cout << "MAJOR ERROR: WHAT" << endl;

    // North Box
    // if (smartAccess(i - 1, j, boundary) == 1.0) {
    //     sum -= 2.0 * (newImg(i, j) - oldImg(i - 1, j));
    // } else if (smartAccess(i - 1, j, omega) == 1.0) {
    //     sum -= 2.0 * (newImg(i, j) - newImg(i - 1, j));
    // } else {
    //     cout << "MAJOR ERROR: OUT OF EXPECTED BOUNDS NORTH" << endl;
    // }


    // // OLD
    //
    // if (boundary.smartAccess(j, i - 1, ch) == 1.0) {
    //     sum -= 2.0 * (imageA.smartAccess(j, i, ch) - imageB.smartAccess(j, i - 1, ch));
    // } else if (omega.smartAccess(j, i - 1, ch) == 1.0) {
    //     sum -= 2.0 * (imageA.smartAccess(j, i, ch) - imageA.smartAccess(j, i - 1, ch));
    // } else {
    //     cout << "MAJOR ERROR: OUT OF EXPECTED BOUNDS NORTH" << endl;
    // }
    //
    // // South Box
    // if (boundary.smartAccess(j, i + 1, ch) == 1.0) {
    //     sum += 2.0 * (imageB.smartAccess(j, i + 1, ch) - imageA.smartAccess(j, i, ch));
    // } else if (omega.smartAccess(j, i + 1, ch) == 1.0) {
    //     sum += 2.0 * (imageA.smartAccess(j, i + 1, ch) - imageA.smartAccess(j, i, ch));
    // } else {
    //     cout << "MAJOR ERROR: OUT OF EXPECTED BOUNDS SOUTH" << endl;
    // }
    //
    // // East Box
    // if (boundary.smartAccess(j + 1, i, ch) == 1.0) {
    //     sum += 2.0 * (imageB.smartAccess(j + 1, i, ch) - imageA.smartAccess(j, i, ch));
    // } else if (omega.smartAccess(j + 1, i, ch) == 1.0) {
    //     sum += 2.0 * (imageA.smartAccess(j + 1, i, ch) - imageA.smartAccess(j, i, ch));
    // } else {
    //     cout << "MAJOR ERROR: OUT OF EXPECTED BOUNDS EAST" << endl;
    // }
    //
    // // West Box
    // if (boundary.smartAccess(j - 1, i, ch) == 1.0) {
    //     sum -= 2.0 * (imageA.smartAccess(j, i, ch) - imageB.smartAccess(j - 1, i, ch));
    // } else if (omega.smartAccess(j - 1, i, ch) == 1.0) {
    //     sum -= 2.0 * (imageA.smartAccess(j, i, ch) - imageA.smartAccess(j - 1, i, ch));
    // } else {
    //     cout << "MAJOR ERROR: OUT OF EXPECTED BOUNDS WEST" << endl;
    // }


    // North Box
    // if (boundary.smartAccess(j, i - 1, ch) == 1.0) {
    //     sum -= 2.0 * (imageA.smartAccess(j, i, ch) - imageA.smartAccess(j, i - 1, ch));
    // } else if (omega.smartAccess(j, i - 1, ch) == 1.0) {
    //     sum -= 2.0 * (imageA.smartAccess(j, i, ch) - imageA.smartAccess(j, i - 1, ch));
    // } else {
    //     cout << "MAJOR ERROR: OUT OF EXPECTED BOUNDS NORTH" << endl;
    // }
    //
    // // South Box
    // if (boundary.smartAccess(j, i + 1, ch) == 1.0) {
    //     sum -= 2.0 * (imageA.smartAccess(j, i, ch) - imageA.smartAccess(j, i + 1, ch));
    // } else if (omega.smartAccess(j, i + 1, ch) == 1.0) {
    //     sum -= 2.0 * (imageA.smartAccess(j, i, ch) - imageA.smartAccess(j, i + 1, ch));
    // } else {
    //     cout << "MAJOR ERROR: OUT OF EXPECTED BOUNDS SOUTH" << endl;
    // }
    //
    // // East Box
    // if (boundary.smartAccess(j + 1, i, ch) == 1.0) {
    //     sum -= 2.0 * (imageA.smartAccess(j, i, ch) - imageA.smartAccess(j + 1, i, ch));
    // } else if (omega.smartAccess(j + 1, i, ch) == 1.0) {
    //     sum -= 2.0 * (imageA.smartAccess(j, i, ch) - imageA.smartAccess(j + 1, i, ch));
    // } else {
    //     cout << "MAJOR ERROR: OUT OF EXPECTED BOUNDS EAST" << endl;
    // }
    //
    // // West Box
    // if (boundary.smartAccess(j - 1, i, ch) == 1.0) {
    //     sum -= 2.0 * (imageA.smartAccess(j, i, ch) - imageA.smartAccess(j - 1, i, ch));
    // } else if (omega.smartAccess(j - 1, i, ch) == 1.0) {
    //     sum -= 2.0 * (imageA.smartAccess(j, i, ch) - imageA.smartAccess(j - 1, i, ch));
    // } else {
    //     cout << "MAJOR ERROR: OUT OF EXPECTED BOUNDS WEST" << endl;
    // }

    sum -= 2.0 * (imageA.smartAccess(j, i, ch) - imageA.smartAccess(j, i - 1, ch));
    sum -= 2.0 * (imageA.smartAccess(j, i, ch) - imageA.smartAccess(j, i + 1, ch));
    sum -= 2.0 * (imageA.smartAccess(j, i, ch) - imageA.smartAccess(j + 1, i, ch));
    sum -= 2.0 * (imageA.smartAccess(j, i, ch) - imageA.smartAccess(j - 1, i, ch));

    // cout << "SUMMM: " << sum << endl;
    return sum;
}

void calculateHessianAndB(Eigen::MatrixXd& A, Eigen::VectorXd& b, const Eigen::MatrixXi& indexMap, const Image& imageB, const Image& imageA, const Image& boundary, const Image& omega, int ch) {
    A.setZero();
    b.setZero();

    for (int i = 0; i < indexMap.rows(); i++) {
        int row = indexMap(i, 0);
        int col = indexMap(i, 1);

        A(i, i) = 8.0;

        if (omega.smartAccess(col + 1, row, ch) > 0.5) {
            int index = findIndex(indexMap, row, col + 1);
            A(i, index) += -2.0;
        } else {
            if (boundary.smartAccess(col + 1, row, ch) > 0.5) {
                b(i) -= 2.0 * imageB.smartAccess(col + 1, row, ch);
            } else {
                cout << "MAJOR ERROR: NOT OMEGA OR BOUNDARY" << endl;
            }
        }

        if (omega.smartAccess(col - 1, row, ch) > 0.5) {
            int index = findIndex(indexMap, row, col - 1);
            A(i, index) += -2.0;
        } else {
            if (boundary.smartAccess(col - 1, row, ch) > 0.5) {
                b(i) -= 2.0 * imageB.smartAccess(col - 1, row, ch);
            } else {
                cout << "MAJOR ERROR: NOT OMEGA OR BOUNDARY" << endl;
            }
        }

        if (omega.smartAccess(col, row - 1, ch) > 0.5) {
            int index = findIndex(indexMap, row - 1, col);
            A(i, index) += -2.0;
        } else {
            if (boundary.smartAccess(col, row - 1, ch) > 0.5) {
                b(i) -= 2.0 * imageB.smartAccess(col, row - 1, ch);
            } else {
                cout << "MAJOR ERROR: NOT OMEGA OR BOUNDARY" << endl;
            }
        }

        if (omega.smartAccess(col, row + 1, ch) > 0.5) {
            int index = findIndex(indexMap, row + 1, col);
            A(i, index) += -2.0;
        } else {
            if (boundary.smartAccess(col, row + 1, ch) > 0.5) {
                b(i) -= 2.0 * imageB.smartAccess(col, row + 1, ch);
            } else {
                cout << "MAJOR ERROR: NOT OMEGA OR BOUNDARY" << endl;
            }
        }

        b(i) += calculateSumOfGradiant(row, col, imageB, imageA, boundary, omega, ch);
    }
}

void calculateHessianAndB(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b, const Eigen::MatrixXi& indexMap, const Image& imageB, const Image& imageA, const Image& boundary, const Image& omega, int ch) {
    b.setZero();

    for (int i = 0; i < indexMap.rows(); i++) {
        int row = indexMap(i, 0);
        int col = indexMap(i, 1);

        A.insert(i, i) = 8.0;

        if (omega.smartAccess(col + 1, row, ch) > 0.5) {
            int index = findIndex(indexMap, row, col + 1);
            A.insert(i, index) = -2.0;
        } else {
            if (boundary.smartAccess(col + 1, row, ch) > 0.5) {
                b(i) -= 2.0 * imageB.smartAccess(col + 1, row, ch);
            } else {
                cout << "MAJOR ERROR: NOT OMEGA OR BOUNDARY" << endl;
            }
        }

        if (omega.smartAccess(col - 1, row, ch) > 0.5) {
            int index = findIndex(indexMap, row, col - 1);
            A.insert(i, index) = -2.0;
        } else {
            if (boundary.smartAccess(col - 1, row, ch) > 0.5) {
                b(i) -= 2.0 * imageB.smartAccess(col - 1, row, ch);
            } else {
                cout << "MAJOR ERROR: NOT OMEGA OR BOUNDARY" << endl;
            }
        }

        if (omega.smartAccess(col, row - 1, ch) > 0.5) {
            int index = findIndex(indexMap, row - 1, col);
            A.insert(i, index) = -2.0;
        } else {
            if (boundary.smartAccess(col, row - 1, ch) > 0.5) {
                b(i) -= 2.0 * imageB.smartAccess(col, row - 1, ch);
            } else {
                cout << "MAJOR ERROR: NOT OMEGA OR BOUNDARY" << endl;
            }
        }

        if (omega.smartAccess(col, row + 1, ch) > 0.5) {
            int index = findIndex(indexMap, row + 1, col);
            A.insert(i, index) = -2.0;
        } else {
            if (boundary.smartAccess(col, row + 1, ch) > 0.5) {
                b(i) -= 2.0 * imageB.smartAccess(col, row + 1, ch);
            } else {
                cout << "MAJOR ERROR: NOT OMEGA OR BOUNDARY" << endl;
            }
        }

        b(i) -= 2.0 * (imageA.smartAccess(col, row, ch) - imageA.smartAccess(col, row - 1, ch));
        b(i) -= 2.0 * (imageA.smartAccess(col, row, ch) - imageA.smartAccess(col, row + 1, ch));
        b(i) -= 2.0 * (imageA.smartAccess(col, row, ch) - imageA.smartAccess(col + 1, row, ch));
        b(i) -= 2.0 * (imageA.smartAccess(col, row, ch) - imageA.smartAccess(col - 1, row, ch));
    }
}

void calculateHessianAndB(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b, const Eigen::MatrixXi& indexGrid, const Eigen::MatrixXi& indexMap, const Image& imageB, const Image& imageA, const Image& boundary, const Image& omega, const Eigen::Vector2i offset, int ch) {
    b.setZero();

    cout << "VARS: " << indexMap.rows() << endl;

    for (int i = 0; i < indexMap.rows(); i++) {
        int row = indexMap(i, 0);
        int col = indexMap(i, 1);

        A.insert(i, i) = 8.0;

        // cout << "HERE" << endl;

        if (omega.smartAccess(col + 1, row, ch) > 0.5) {
            // int index = findIndexFast(indexMap, row, col + 1);
            int index = indexGrid(row, col + 1);
            // cout << index << endl;
            A.insert(i, index) = -2.0;
        } else {
            if (boundary.smartAccess(col + 1, row, ch) > 0.5) {
                b(i) -= 2.0 * imageB.smartAccess(col + 1, row, offset, ch);
            } else {
                cout << "MAJOR ERROR: NOT OMEGA OR BOUNDARY" << endl;
            }
        }

        if (omega.smartAccess(col - 1, row, ch) > 0.5) {
            // int index = findIndexFast(indexMap, row, col - 1);
            int index = indexGrid(row, col - 1);
            A.insert(i, index) = -2.0;
        } else {
            if (boundary.smartAccess(col - 1, row, ch) > 0.5) {
                b(i) -= 2.0 * imageB.smartAccess(col - 1, row, offset, ch);
            } else {
                cout << "MAJOR ERROR: NOT OMEGA OR BOUNDARY" << endl;
            }
        }

        if (omega.smartAccess(col, row - 1, ch) > 0.5) {
            // int index = findIndexFast(indexMap, row - 1, col);
            int index = indexGrid(row - 1, col);
            A.insert(i, index) = -2.0;
        } else {
            if (boundary.smartAccess(col, row - 1, ch) > 0.5) {
                b(i) -= 2.0 * imageB.smartAccess(col, row - 1, offset, ch);
            } else {
                cout << "MAJOR ERROR: NOT OMEGA OR BOUNDARY" << endl;
            }
        }

        if (omega.smartAccess(col, row + 1, ch) > 0.5) {
            // int index = findIndexFast(indexMap, row + 1, col);
            int index = indexGrid(row + 1, col);
            A.insert(i, index) = -2.0;
        } else {
            if (boundary.smartAccess(col, row + 1, ch) > 0.5) {
                b(i) -= 2.0 * imageB.smartAccess(col, row + 1, offset, ch);
            } else {
                cout << "MAJOR ERROR: NOT OMEGA OR BOUNDARY" << endl;
            }
        }

        b(i) -= 2.0 * (imageA.smartAccess(col, row, ch) - imageA.smartAccess(col, row - 1, ch));
        b(i) -= 2.0 * (imageA.smartAccess(col, row, ch) - imageA.smartAccess(col, row + 1, ch));
        b(i) -= 2.0 * (imageA.smartAccess(col, row, ch) - imageA.smartAccess(col + 1, row, ch));
        b(i) -= 2.0 * (imageA.smartAccess(col, row, ch) - imageA.smartAccess(col - 1, row, ch));
    }
}

void calculateHessianAndB(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b, const Eigen::MatrixXi& indexGrid, const Eigen::MatrixXi& indexMap, const Image& imageB, const Image& imageA, const Image& boundary, const Image& omega, const Image& edge, const Eigen::Vector2i offset, int ch) {
    b.setZero();

    cout << "VARS: " << indexMap.rows() << endl;

    for (int i = 0; i < indexMap.rows(); i++) {
        int row = indexMap(i, 0);
        int col = indexMap(i, 1);

        A.insert(i, i) = 8.0;

        // cout << "HERE" << endl;

        if (omega.smartAccess(col + 1, row, ch) > 0.5) {
            // int index = findIndexFast(indexMap, row, col + 1);
            int index = indexGrid(row, col + 1);
            // cout << index << endl;
            A.insert(i, index) = -2.0;
        } else {
            if (boundary.smartAccess(col + 1, row, ch) > 0.5) {
                b(i) -= 2.0 * imageB.smartAccess(col + 1, row, offset, ch);
            } else {
                if (row == 0 || col == 0 || row == imageB.rows() - 1 || col == imageB.cols() - 1) { }
                else {
                    cout << "MAJOR ERROR: NOT OMEGA OR BOUNDARY" << endl;
                }
            }
        }

        if (omega.smartAccess(col - 1, row, ch) > 0.5) {
            // int index = findIndexFast(indexMap, row, col - 1);
            int index = indexGrid(row, col - 1);
            A.insert(i, index) = -2.0;
        } else {
            if (boundary.smartAccess(col - 1, row, ch) > 0.5) {
                b(i) -= 2.0 * imageB.smartAccess(col - 1, row, offset, ch);
            } else {
                if (row == 0 || col == 0 || row == imageB.rows() - 1 || col == imageB.cols() - 1) { }
                else {
                    cout << "MAJOR ERROR: NOT OMEGA OR BOUNDARY" << endl;
                }
            }
        }

        if (omega.smartAccess(col, row - 1, ch) > 0.5) {
            // int index = findIndexFast(indexMap, row - 1, col);
            int index = indexGrid(row - 1, col);
            A.insert(i, index) = -2.0;
        } else {
            if (boundary.smartAccess(col, row - 1, ch) > 0.5) {
                b(i) -= 2.0 * imageB.smartAccess(col, row - 1, offset, ch);
            } else {
                if (row == 0 || col == 0 || row == imageB.rows() - 1 || col == imageB.cols() - 1) { }
                else {
                    cout << "MAJOR ERROR: NOT OMEGA OR BOUNDARY" << endl;
                }
            }
        }

        if (omega.smartAccess(col, row + 1, ch) > 0.5) {
            // int index = findIndexFast(indexMap, row + 1, col);
            int index = indexGrid(row + 1, col);
            A.insert(i, index) = -2.0;
        } else {
            if (boundary.smartAccess(col, row + 1, ch) > 0.5) {
                b(i) -= 2.0 * imageB.smartAccess(col, row + 1, offset, ch);
            } else {
                if (row == 0 || col == 0 || row == imageB.rows() - 1 || col == imageB.cols() - 1) { }
                else {
                    cout << "MAJOR ERROR: NOT OMEGA OR BOUNDARY" << endl;
                }
            }
        }

        if (edge.smartAccess(col, row, ch) > 0.1) {
            b(i) -= 2.0 * (imageA.smartAccess(col, row, ch) - imageA.smartAccess(col, row - 1, ch));
            b(i) -= 2.0 * (imageA.smartAccess(col, row, ch) - imageA.smartAccess(col, row + 1, ch));
            b(i) -= 2.0 * (imageA.smartAccess(col, row, ch) - imageA.smartAccess(col + 1, row, ch));
            b(i) -= 2.0 * (imageA.smartAccess(col, row, ch) - imageA.smartAccess(col - 1, row, ch));
        }
    }
}

double absf(double val) {
    if (val < 0) return -val;
    return val;
}

void calculateHessianAndB(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b, const Eigen::MatrixXi& indexGrid, const Eigen::MatrixXi& indexMap, const Image& imageB, const Image& imageA, const Image& boundary, const Image& omega, double edgeTh, const Eigen::Vector2i offset, int ch) {
    b.setZero();

    cout << "VARS: " << indexMap.rows() << endl;

    for (int i = 0; i < indexMap.rows(); i++) {
        int row = indexMap(i, 0);
        int col = indexMap(i, 1);

        A.insert(i, i) = 8.0;

        // cout << "HERE" << endl;

        if (omega.smartAccess(col + 1, row, ch) > 0.5) {
            // int index = findIndexFast(indexMap, row, col + 1);
            int index = indexGrid(row, col + 1);
            // cout << index << endl;
            A.insert(i, index) = -2.0;
        } else {
            if (boundary.smartAccess(col + 1, row, ch) > 0.5) {
                b(i) -= 2.0 * imageB.smartAccess(col + 1, row, offset, ch);
            } else {
                if (row == 0 || col == 0 || row == imageB.rows() - 1 || col == imageB.cols() - 1) { }
                else {
                    cout << "MAJOR ERROR: NOT OMEGA OR BOUNDARY" << endl;
                }
            }
        }

        if (omega.smartAccess(col - 1, row, ch) > 0.5) {
            // int index = findIndexFast(indexMap, row, col - 1);
            int index = indexGrid(row, col - 1);
            A.insert(i, index) = -2.0;
        } else {
            if (boundary.smartAccess(col - 1, row, ch) > 0.5) {
                b(i) -= 2.0 * imageB.smartAccess(col - 1, row, offset, ch);
            } else {
                if (row == 0 || col == 0 || row == imageB.rows() - 1 || col == imageB.cols() - 1) { }
                else {
                    cout << "MAJOR ERROR: NOT OMEGA OR BOUNDARY" << endl;
                }
            }
        }

        if (omega.smartAccess(col, row - 1, ch) > 0.5) {
            // int index = findIndexFast(indexMap, row - 1, col);
            int index = indexGrid(row - 1, col);
            A.insert(i, index) = -2.0;
        } else {
            if (boundary.smartAccess(col, row - 1, ch) > 0.5) {
                b(i) -= 2.0 * imageB.smartAccess(col, row - 1, offset, ch);
            } else {
                if (row == 0 || col == 0 || row == imageB.rows() - 1 || col == imageB.cols() - 1) { }
                else {
                    cout << "MAJOR ERROR: NOT OMEGA OR BOUNDARY" << endl;
                }
            }
        }

        if (omega.smartAccess(col, row + 1, ch) > 0.5) {
            // int index = findIndexFast(indexMap, row + 1, col);
            int index = indexGrid(row + 1, col);
            A.insert(i, index) = -2.0;
        } else {
            if (boundary.smartAccess(col, row + 1, ch) > 0.5) {
                b(i) -= 2.0 * imageB.smartAccess(col, row + 1, offset, ch);
            } else {
                if (row == 0 || col == 0 || row == imageB.rows() - 1 || col == imageB.cols() - 1) { }
                else {
                    cout << "MAJOR ERROR: NOT OMEGA OR BOUNDARY" << endl;
                }
            }
        }

        double one = 2.0 * (imageA.smartAccess(col, row, ch) - imageA.smartAccess(col, row - 1, ch));
        double two = 2.0 * (imageA.smartAccess(col, row, ch) - imageA.smartAccess(col, row + 1, ch));
        double three = 2.0 * (imageA.smartAccess(col, row, ch) - imageA.smartAccess(col + 1, row, ch));
        double four = 2.0 * (imageA.smartAccess(col, row, ch) - imageA.smartAccess(col - 1, row, ch));

        if (fabs(one) > edgeTh) b(i) -= one;
        if (fabs(two) > edgeTh) b(i) -= two;
        if (fabs(three) > edgeTh) b(i) -= three;
        if (fabs(four) > edgeTh) b(i) -= four;

        // if (edge.smartAccess(col, row, ch) > 0.1) {
        //     b(i) -= 2.0 * (imageA.smartAccess(col, row, ch) - imageA.smartAccess(col, row - 1, ch));
        //     b(i) -= 2.0 * (imageA.smartAccess(col, row, ch) - imageA.smartAccess(col, row + 1, ch));
        //     b(i) -= 2.0 * (imageA.smartAccess(col, row, ch) - imageA.smartAccess(col + 1, row, ch));
        //     b(i) -= 2.0 * (imageA.smartAccess(col, row, ch) - imageA.smartAccess(col - 1, row, ch));
        // }
    }
}

void calculateHessianAndBLog(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b, const Eigen::MatrixXi& indexGrid, const Eigen::MatrixXi& indexMap, const Image& imageB, const Image& imageA, const Image& boundary, const Image& omega, const Eigen::Vector2i offset, int ch) {
    b.setZero();

    cout << "VARS: " << indexMap.rows() << endl;

    for (int i = 0; i < indexMap.rows(); i++) {
        int row = indexMap(i, 0);
        int col = indexMap(i, 1);

        A.insert(i, i) = log10(8.0);

        // cout << "HERE" << endl;

        if (omega.smartAccess(col + 1, row, ch) > 0.5) {
            // int index = findIndexFast(indexMap, row, col + 1);
            int index = indexGrid(row, col + 1);
            // cout << index << endl;
            A.insert(i, index) = log10(-2.0);
        } else {
            if (boundary.smartAccess(col + 1, row, ch) > 0.5) {
                b(i) -= 2.0 * imageB.smartAccess(col + 1, row, offset, ch);
            } else {
                cout << "MAJOR ERROR: NOT OMEGA OR BOUNDARY" << endl;
            }
        }

        if (omega.smartAccess(col - 1, row, ch) > 0.5) {
            // int index = findIndexFast(indexMap, row, col - 1);
            int index = indexGrid(row, col - 1);
            A.insert(i, index) = log10(-2.0);
        } else {
            if (boundary.smartAccess(col - 1, row, ch) > 0.5) {
                b(i) -= 2.0 * imageB.smartAccess(col - 1, row, offset, ch);
            } else {
                cout << "MAJOR ERROR: NOT OMEGA OR BOUNDARY" << endl;
            }
        }

        if (omega.smartAccess(col, row - 1, ch) > 0.5) {
            // int index = findIndexFast(indexMap, row - 1, col);
            int index = indexGrid(row - 1, col);
            A.insert(i, index) = log10(-2.0);
        } else {
            if (boundary.smartAccess(col, row - 1, ch) > 0.5) {
                b(i) -= 2.0 * imageB.smartAccess(col, row - 1, offset, ch);
            } else {
                cout << "MAJOR ERROR: NOT OMEGA OR BOUNDARY" << endl;
            }
        }

        if (omega.smartAccess(col, row + 1, ch) > 0.5) {
            // int index = findIndexFast(indexMap, row + 1, col);
            int index = indexGrid(row + 1, col);
            A.insert(i, index) = log10(-2.0);
        } else {
            if (boundary.smartAccess(col, row + 1, ch) > 0.5) {
                b(i) -= 2.0 * imageB.smartAccess(col, row + 1, offset, ch);
            } else {
                cout << "MAJOR ERROR: NOT OMEGA OR BOUNDARY" << endl;
            }
        }

        b(i) -= 2.0 * (imageA.smartAccess(col, row, ch) - imageA.smartAccess(col, row - 1, ch));
        b(i) -= 2.0 * (imageA.smartAccess(col, row, ch) - imageA.smartAccess(col, row + 1, ch));
        b(i) -= 2.0 * (imageA.smartAccess(col, row, ch) - imageA.smartAccess(col + 1, row, ch));
        b(i) -= 2.0 * (imageA.smartAccess(col, row, ch) - imageA.smartAccess(col - 1, row, ch));

        b(i) = log10(b(i));
    }
}

void imageExample() {
    cout << "Image Example" << endl;

    // Image imageA(DATA_DIR "/input/GradiantX.png");
    // Image imageB(DATA_DIR "/input/GradiantY.png");
    Image imageB(DATA_DIR "/input/GradiantX.png");
    Image imageA(DATA_DIR "/input/GradiantY.png");
    Image omega(DATA_DIR "/input/GradiantReplace1.png");
    Image boundary(imageB.cols(), imageB.rows(), 3);
    Image results(imageB.cols(), imageB.rows(), 3);

    boundary = calculateBoundaryMap(omega);
    results.setZero();

    cout << "OMEGA RED:" << endl;
    omega.printRed();

    cout << "BOUNDARY RED:" << endl;
    boundary.printRed();

    cout << "OMEGA Green:" << endl;
    omega.printGreen();

    cout << "BOUNDARY Green:" << endl;
    boundary.printGreen();

    cout << "OMEGA Blue:" << endl;
    omega.printBlue();

    cout << "BOUNDARY Blue:" << endl;
    boundary.printBlue();

    int redVars;
    int greenVars;
    int blueVars;

    countVars(omega, redVars, greenVars, blueVars);

    cout << "Red Count: " << redVars << endl;
    cout << "Green Count: " << greenVars << endl;
    cout << "Blue Count: " << blueVars << endl;

    Eigen::MatrixXi redIndexMap;
    Eigen::MatrixXi greenIndexMap;
    Eigen::MatrixXi blueIndexMap;

    redIndexMap.resize(redVars, 2);
    greenIndexMap.resize(greenVars, 2);
    blueIndexMap.resize(blueVars, 2);

    calculateIndexMaps(omega, redIndexMap, greenIndexMap, blueIndexMap);

    cout << "Red Index Map: " << endl << redIndexMap << endl;

    Eigen::VectorXd redX;
    Eigen::VectorXd greenX;
    Eigen::VectorXd blueX;

    redX.resize(redVars);
    greenX.resize(greenVars);
    blueX.resize(blueVars);

    redX.setZero();
    greenX.setZero();
    blueX.setZero();

    for (int i = 0; i < redVars; i++) {
        redX(i) = omega.smartAccess(redIndexMap(i, 1), redIndexMap(i, 0), 0);
    }

    for (int i = 0; i < greenVars; i++) {
        greenX(i) = omega.smartAccess(greenIndexMap(i, 1), greenIndexMap(i, 0), 1);
    }

    for (int i = 0; i < blueVars; i++) {
        blueX(i) = omega.smartAccess(blueIndexMap(i, 1), blueIndexMap(i, 0), 2);
    }

    Eigen::MatrixXd redA;
    Eigen::MatrixXd greenA;
    Eigen::MatrixXd blueA;

    Eigen::VectorXd redB;
    Eigen::VectorXd greenB;
    Eigen::VectorXd blueB;

    redA.resize(redVars, redVars);
    greenA.resize(greenVars, greenVars);
    blueA.resize(blueVars, blueVars);

    redB.resize(redVars);
    greenB.resize(greenVars);
    blueB.resize(blueVars);

    // const Eigen::MatrixXd& oldImg, const Eigen::MatrixXd& newImg, const Eigen::MatrixXd& boundary, const Eigen::MatrixXd& omega

    calculateHessianAndB(redA, redB, redIndexMap, imageB, imageA, boundary, omega, 0);
    calculateHessianAndB(greenA, greenB, greenIndexMap, imageB, imageA, boundary, omega, 1);
    calculateHessianAndB(blueA, blueB, blueIndexMap, imageB, imageA, boundary, omega, 2);

    // cout << "Red A:" << endl;
    //
    // cout << redA << endl;

    Eigen::BiCGSTAB<Eigen::MatrixXd> redSolver;
    Eigen::BiCGSTAB<Eigen::MatrixXd> greenSolver;
    Eigen::BiCGSTAB<Eigen::MatrixXd> blueSolver;

    redSolver.compute(redA);
    greenSolver.compute(greenA);
    blueSolver.compute(blueA);

    redX.setZero();
    greenX.setZero();
    blueX.setZero();

    redX = redSolver.solve(-redB);
    greenX = greenSolver.solve(-greenB);
    blueX = blueSolver.solve(-blueB);


    cout << "X: " << endl << redX << endl;

    std::cout << "#iterations:     " << redSolver.iterations() << std::endl;
    std::cout << "estimated error: " << redSolver.error()      << std::endl;

    for (int c = 0; c < 3; c++) {
        for (int i = 0; i < results.rows(); i++) {
            for (int j = 0; j < results.cols(); j++) {
                results.set(j, i, c, imageB.smartAccess(j, i, c));
            }
        }
    }

    for (int i = 0; i < redVars; i++) {
        results.set(redIndexMap(i, 1), redIndexMap(i, 0), 0, redX(i));
    }

    for (int i = 0; i < greenVars; i++) {
        results.set(greenIndexMap(i, 1), greenIndexMap(i, 0), 1, greenX(i));
    }

    for (int i = 0; i < blueVars; i++) {
        results.set(blueIndexMap(i, 1), blueIndexMap(i, 0), 2, blueX(i));
    }

    results.write(DATA_DIR "/output/gradient1Test.png");

    cout << "Grad Sum: " << calculateSumOfGradiant(10, 11, imageB, imageA, boundary, omega, 0) << endl;
    cout << "P: " << imageA.smartAccess(11, 10, 0) << endl;
    cout << "N: " << imageA.smartAccess(11, 9, 0) << endl;
    cout << "S: " << imageA.smartAccess(11, 11, 0) << endl;
    cout << "E: " << imageA.smartAccess(12, 10, 0) << endl;
    cout << "W: " << imageA.smartAccess(10, 10, 0) << endl;

    double P = imageA.smartAccess(11, 10, 0);
    double N = imageA.smartAccess(11, 9, 0);
    double S = imageA.smartAccess(11, 11, 0);
    double E = imageA.smartAccess(12, 10, 0);
    double W = imageA.smartAccess(10, 10, 0);

    double BN = imageB.smartAccess(11, 9, 0);
    double BW = imageB.smartAccess(10, 10, 0);

    cout << "ACTUAL GRAD: " << -2.0 * (P - N) - 2.0 * (P - W) + 2.0 * (E - P) + 2.0 * (S - P) << endl;

    cout << "B: " << redB(0) << endl;

    cout << "ACTUAL B: " << -2.0 * (P - N) - 2.0 * (P - W) + 2.0 * (E - P) + 2.0 * (S - P) - 2.0 * BN - 2.0 * BW << endl;
}

void dogDemo() {
    cout << "Dog Example" << endl;

    // Image imageA(DATA_DIR "/input/GradiantX.png");
    // Image imageB(DATA_DIR "/input/GradiantY.png");
    Image imageB(DATA_DIR "/input/water.jpg");
    Image imageA(DATA_DIR "/input/doggo.jpg");
    Image omega(DATA_DIR "/input/dogweights.png");
    Image boundary(omega.cols(), omega.rows(), 3);
    Image results(imageB.cols(), imageB.rows(), 3);

    cout << "Creating Boundary" << endl;

    boundary = calculateBoundaryMap(omega);
    results.setZero();

    // cout << "OMEGA RED:" << endl;
    // omega.printRed();
    //
    // cout << "BOUNDARY RED:" << endl;
    // boundary.printRed();
    //
    // cout << "OMEGA Green:" << endl;
    // omega.printGreen();
    //
    // cout << "BOUNDARY Green:" << endl;
    // boundary.printGreen();
    //
    // cout << "OMEGA Blue:" << endl;
    // omega.printBlue();
    //
    // cout << "BOUNDARY Blue:" << endl;
    // boundary.printBlue();

    int redVars;
    int greenVars;
    int blueVars;

    cout << "Counting Vars" << endl;

    countVars(omega, redVars, greenVars, blueVars);

    // cout << "Red Count: " << redVars << endl;
    // cout << "Green Count: " << greenVars << endl;
    // cout << "Blue Count: " << blueVars << endl;

    Eigen::MatrixXi redIndexMap;
    Eigen::MatrixXi greenIndexMap;
    Eigen::MatrixXi blueIndexMap;

    redIndexMap.resize(redVars, 2);
    greenIndexMap.resize(greenVars, 2);
    blueIndexMap.resize(blueVars, 2);

    cout << "Calculating Index Maps" << endl;

    calculateIndexMaps(omega, redIndexMap, greenIndexMap, blueIndexMap);

    // cout << "Red Index Map: " << endl << redIndexMap << endl;

    Eigen::VectorXd redX;
    Eigen::VectorXd greenX;
    Eigen::VectorXd blueX;

    redX.resize(redVars);
    greenX.resize(greenVars);
    blueX.resize(blueVars);

    redX.setZero();
    greenX.setZero();
    blueX.setZero();

    cout << "Initializing X" << endl;

    for (int i = 0; i < redVars; i++) {
        redX(i) = omega.smartAccess(redIndexMap(i, 1), redIndexMap(i, 0), 0);
    }

    for (int i = 0; i < greenVars; i++) {
        greenX(i) = omega.smartAccess(greenIndexMap(i, 1), greenIndexMap(i, 0), 1);
    }

    for (int i = 0; i < blueVars; i++) {
        blueX(i) = omega.smartAccess(blueIndexMap(i, 1), blueIndexMap(i, 0), 2);
    }

    Eigen::MatrixXd redA;
    Eigen::MatrixXd greenA;
    Eigen::MatrixXd blueA;

    Eigen::VectorXd redB;
    Eigen::VectorXd greenB;
    Eigen::VectorXd blueB;

    redA.resize(redVars, redVars);
    greenA.resize(greenVars, greenVars);
    blueA.resize(blueVars, blueVars);

    redB.resize(redVars);
    greenB.resize(greenVars);
    blueB.resize(blueVars);

    cout << "Calculating Hessians" << endl;

    // const Eigen::MatrixXd& oldImg, const Eigen::MatrixXd& newImg, const Eigen::MatrixXd& boundary, const Eigen::MatrixXd& omega

    calculateHessianAndB(redA, redB, redIndexMap, imageB, imageA, boundary, omega, 0);
    calculateHessianAndB(greenA, greenB, greenIndexMap, imageB, imageA, boundary, omega, 1);
    calculateHessianAndB(blueA, blueB, blueIndexMap, imageB, imageA, boundary, omega, 2);

    // cout << "Red A:" << endl;
    //
    // cout << redA << endl;

    Eigen::BiCGSTAB<Eigen::MatrixXd> redSolver;
    Eigen::BiCGSTAB<Eigen::MatrixXd> greenSolver;
    Eigen::BiCGSTAB<Eigen::MatrixXd> blueSolver;

    redSolver.compute(redA);
    greenSolver.compute(greenA);
    blueSolver.compute(blueA);

    // redX.setZero();
    // greenX.setZero();
    // blueX.setZero();

    cout << "Solving" << endl;

    redX = redSolver.solve(-redB);
    greenX = greenSolver.solve(-greenB);
    blueX = blueSolver.solve(-blueB);

    // cout << "X: " << endl << redX << endl;

    std::cout << "#iterations:     " << redSolver.iterations() << std::endl;
    std::cout << "estimated error: " << redSolver.error()      << std::endl;

    for (int c = 0; c < 3; c++) {
        for (int i = 0; i < results.rows(); i++) {
            for (int j = 0; j < results.cols(); j++) {
                results.set(j, i, c, imageB.smartAccess(j, i, c));
            }
        }
    }

    cout << "Writing Image" << endl;

    for (int i = 0; i < redVars; i++) {
        results.set(redIndexMap(i, 1), redIndexMap(i, 0), 0, redX(i));
    }

    for (int i = 0; i < greenVars; i++) {
        results.set(greenIndexMap(i, 1), greenIndexMap(i, 0), 1, greenX(i));
    }

    for (int i = 0; i < blueVars; i++) {
        results.set(blueIndexMap(i, 1), blueIndexMap(i, 0), 2, blueX(i));
    }

    results.write(DATA_DIR "/output/doggoTest.png");

    // cout << "Grad Sum: " << calculateSumOfGradiant(10, 11, imageB, imageA, boundary, omega, 0) << endl;
    // cout << "P: " << imageA.smartAccess(11, 10, 0) << endl;
    // cout << "N: " << imageA.smartAccess(11, 9, 0) << endl;
    // cout << "S: " << imageA.smartAccess(11, 11, 0) << endl;
    // cout << "E: " << imageA.smartAccess(12, 10, 0) << endl;
    // cout << "W: " << imageA.smartAccess(10, 10, 0) << endl;
    //
    // double P = imageA.smartAccess(11, 10, 0);
    // double N = imageA.smartAccess(11, 9, 0);
    // double S = imageA.smartAccess(11, 11, 0);
    // double E = imageA.smartAccess(12, 10, 0);
    // double W = imageA.smartAccess(10, 10, 0);
    //
    // double BN = imageB.smartAccess(11, 9, 0);
    // double BW = imageB.smartAccess(10, 10, 0);
    //
    // cout << "ACTUAL GRAD: " << -2.0 * (P - N) - 2.0 * (P - W) + 2.0 * (E - P) + 2.0 * (S - P) << endl;
    //
    // cout << "B: " << redB(0) << endl;
    //
    // cout << "ACTUAL B: " << -2.0 * (P - N) - 2.0 * (P - W) + 2.0 * (E - P) + 2.0 * (S - P) - 2.0 * BN - 2.0 * BW << endl;
}

// imageA is the source image
// imageB is the dest image
Image seamlessPoissonCloning(const Image& imageA, const Image& imageB, const Image& omega, Eigen::Vector2i omegaOffset) {
    Image results(imageB.cols(), imageB.rows(), 3);

    for (int c = 0; c < 3; c++) {
        for (int i = 0; i < results.rows(); i++) {
            for (int j = 0; j < results.cols(); j++) {
                results.set(j, i, c, imageB.smartAccess(j, i, c));
            }
        }
    }

    Image boundary(omega.cols(), omega.rows(), 3);

    // cout << "Creating Boundary" << endl;
    boundary = calculateBoundaryMap(omega);

    int redVars;
    int greenVars;
    int blueVars;

    // cout << "Counting Vars" << endl;
    countVars(omega, redVars, greenVars, blueVars);

    Eigen::MatrixXi redIndexGrid;
    Eigen::MatrixXi greenIndexGrid;
    Eigen::MatrixXi blueIndexGrid;

    Eigen::MatrixXi redIndexMap;
    Eigen::MatrixXi greenIndexMap;
    Eigen::MatrixXi blueIndexMap;

    redIndexGrid.resize(omega.cols(), omega.rows());
    greenIndexGrid.resize(omega.cols(), omega.rows());
    blueIndexGrid.resize(omega.cols(), omega.rows());

    redIndexMap.resize(redVars, 2);
    greenIndexMap.resize(greenVars, 2);
    blueIndexMap.resize(blueVars, 2);

    // cout << "Calculating Index Maps" << endl;
    calculateIndexMaps(omega, redIndexMap, greenIndexMap, blueIndexMap);
    calculateIndexMapsFast(omega, redIndexGrid, greenIndexGrid, blueIndexGrid);

    Eigen::VectorXd redX;
    Eigen::VectorXd greenX;
    Eigen::VectorXd blueX;

    redX.resize(redVars);
    greenX.resize(greenVars);
    blueX.resize(blueVars);

    redX.setZero();
    greenX.setZero();
    blueX.setZero();

    // cout << "Initializing X" << endl;
    for (int i = 0; i < redVars; i++) {
        redX(i) = omega.smartAccess(redIndexMap(i, 1), redIndexMap(i, 0), 0);
    }
    for (int i = 0; i < greenVars; i++) {
        greenX(i) = omega.smartAccess(greenIndexMap(i, 1), greenIndexMap(i, 0), 1);
    }
    for (int i = 0; i < blueVars; i++) {
        blueX(i) = omega.smartAccess(blueIndexMap(i, 1), blueIndexMap(i, 0), 2);
    }

    Eigen::SparseMatrix<double> redA;
    Eigen::SparseMatrix<double> greenA;
    Eigen::SparseMatrix<double> blueA;

    Eigen::VectorXd redB;
    Eigen::VectorXd greenB;
    Eigen::VectorXd blueB;

    redA.resize(redVars, redVars);
    greenA.resize(greenVars, greenVars);
    blueA.resize(blueVars, blueVars);

    redB.resize(redVars);
    greenB.resize(greenVars);
    blueB.resize(blueVars);

    // cout << "Calculating Hessians" << endl;
    calculateHessianAndB(redA, redB, redIndexGrid, redIndexMap, imageB, imageA, boundary, omega, omegaOffset, 0);
    calculateHessianAndB(greenA, greenB, greenIndexGrid, greenIndexMap, imageB, imageA, boundary, omega, omegaOffset, 1);
    calculateHessianAndB(blueA, blueB, blueIndexGrid, blueIndexMap, imageB, imageA, boundary, omega, omegaOffset, 2);

    // Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > redSolver;
    // Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > greenSolver;
    // Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > blueSolver;

    Eigen::ConjugateGradient< Eigen::SparseMatrix<double> > redSolver;
    Eigen::ConjugateGradient< Eigen::SparseMatrix<double> > greenSolver;
    Eigen::ConjugateGradient< Eigen::SparseMatrix<double> > blueSolver;

    redSolver.compute(redA);
    greenSolver.compute(greenA);
    blueSolver.compute(blueA);

    redX = redSolver.solve(-redB);
    greenX = greenSolver.solve(-greenB);
    blueX = blueSolver.solve(-blueB);

    cout << "ERROR: " << (redA * redX + redB).norm() << endl;

    for (int i = 0; i < redVars; i++) {
        results.set(redIndexMap(i, 1) + omegaOffset[0], redIndexMap(i, 0) + omegaOffset[1], 0, redX(i));
    }

    for (int i = 0; i < greenVars; i++) {
        results.set(greenIndexMap(i, 1) + omegaOffset[0], greenIndexMap(i, 0) + omegaOffset[1], 1, greenX(i));
    }

    for (int i = 0; i < blueVars; i++) {
        results.set(blueIndexMap(i, 1) + omegaOffset[0], blueIndexMap(i, 0) + omegaOffset[1], 2, blueX(i));
    }

    return results;
}

Image seamlessPoissonCloningToon(const Image& imageA, const Image& imageB, const Image& omega, const Image& edge, Eigen::Vector2i omegaOffset) {
    Image results(imageB.cols(), imageB.rows(), 3);

    for (int c = 0; c < 3; c++) {
        for (int i = 0; i < results.rows(); i++) {
            for (int j = 0; j < results.cols(); j++) {
                results.set(j, i, c, imageB.smartAccess(j, i, c));
            }
        }
    }

    Image boundary(omega.cols(), omega.rows(), 3);

    // cout << "Creating Boundary" << endl;
    boundary = calculateBoundaryMap(omega);

    int redVars;
    int greenVars;
    int blueVars;

    // cout << "Counting Vars" << endl;
    countVars(omega, redVars, greenVars, blueVars);

    Eigen::MatrixXi redIndexGrid;
    Eigen::MatrixXi greenIndexGrid;
    Eigen::MatrixXi blueIndexGrid;

    Eigen::MatrixXi redIndexMap;
    Eigen::MatrixXi greenIndexMap;
    Eigen::MatrixXi blueIndexMap;

    redIndexGrid.resize(omega.cols(), omega.rows());
    greenIndexGrid.resize(omega.cols(), omega.rows());
    blueIndexGrid.resize(omega.cols(), omega.rows());

    redIndexMap.resize(redVars, 2);
    greenIndexMap.resize(greenVars, 2);
    blueIndexMap.resize(blueVars, 2);

    // cout << "Calculating Index Maps" << endl;
    calculateIndexMaps(omega, redIndexMap, greenIndexMap, blueIndexMap);
    calculateIndexMapsFast(omega, redIndexGrid, greenIndexGrid, blueIndexGrid);

    Eigen::VectorXd redX;
    Eigen::VectorXd greenX;
    Eigen::VectorXd blueX;

    redX.resize(redVars);
    greenX.resize(greenVars);
    blueX.resize(blueVars);

    redX.setZero();
    greenX.setZero();
    blueX.setZero();

    // cout << "Initializing X" << endl;
    for (int i = 0; i < redVars; i++) {
        redX(i) = omega.smartAccess(redIndexMap(i, 1), redIndexMap(i, 0), 0);
    }
    for (int i = 0; i < greenVars; i++) {
        greenX(i) = omega.smartAccess(greenIndexMap(i, 1), greenIndexMap(i, 0), 1);
    }
    for (int i = 0; i < blueVars; i++) {
        blueX(i) = omega.smartAccess(blueIndexMap(i, 1), blueIndexMap(i, 0), 2);
    }

    Eigen::SparseMatrix<double> redA;
    Eigen::SparseMatrix<double> greenA;
    Eigen::SparseMatrix<double> blueA;

    Eigen::VectorXd redB;
    Eigen::VectorXd greenB;
    Eigen::VectorXd blueB;

    redA.resize(redVars, redVars);
    greenA.resize(greenVars, greenVars);
    blueA.resize(blueVars, blueVars);

    redB.resize(redVars);
    greenB.resize(greenVars);
    blueB.resize(blueVars);

    // cout << "Calculating Hessians" << endl;
    calculateHessianAndB(redA, redB, redIndexGrid, redIndexMap, imageB, imageA, boundary, omega, edge, omegaOffset, 0);
    calculateHessianAndB(greenA, greenB, greenIndexGrid, greenIndexMap, imageB, imageA, boundary, omega, edge, omegaOffset, 1);
    calculateHessianAndB(blueA, blueB, blueIndexGrid, blueIndexMap, imageB, imageA, boundary, omega, edge, omegaOffset, 2);

    // Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > redSolver;
    // Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > greenSolver;
    // Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > blueSolver;

    Eigen::ConjugateGradient< Eigen::SparseMatrix<double> > redSolver;
    Eigen::ConjugateGradient< Eigen::SparseMatrix<double> > greenSolver;
    Eigen::ConjugateGradient< Eigen::SparseMatrix<double> > blueSolver;

    redSolver.compute(redA);
    greenSolver.compute(greenA);
    blueSolver.compute(blueA);

    redX = redSolver.solve(-redB);
    greenX = greenSolver.solve(-greenB);
    blueX = blueSolver.solve(-blueB);

    cout << "ERROR: " << (redA * redX + redB).norm() << endl;

    for (int i = 0; i < redVars; i++) {
        results.set(redIndexMap(i, 1) + omegaOffset[0], redIndexMap(i, 0) + omegaOffset[1], 0, redX(i));
    }

    for (int i = 0; i < greenVars; i++) {
        results.set(greenIndexMap(i, 1) + omegaOffset[0], greenIndexMap(i, 0) + omegaOffset[1], 1, greenX(i));
    }

    for (int i = 0; i < blueVars; i++) {
        results.set(blueIndexMap(i, 1) + omegaOffset[0], blueIndexMap(i, 0) + omegaOffset[1], 2, blueX(i));
    }

    return results;
}

Image seamlessPoissonCloningToon(const Image& imageA, const Image& imageB, const Image& omega, double edgeTh, Eigen::Vector2i omegaOffset) {
    Image results(imageB.cols(), imageB.rows(), 3);

    for (int c = 0; c < 3; c++) {
        for (int i = 0; i < results.rows(); i++) {
            for (int j = 0; j < results.cols(); j++) {
                results.set(j, i, c, imageB.smartAccess(j, i, c));
            }
        }
    }

    Image boundary(omega.cols(), omega.rows(), 3);

    // cout << "Creating Boundary" << endl;
    boundary = calculateBoundaryMap(omega);

    int redVars;
    int greenVars;
    int blueVars;

    // cout << "Counting Vars" << endl;
    countVars(omega, redVars, greenVars, blueVars);

    Eigen::MatrixXi redIndexGrid;
    Eigen::MatrixXi greenIndexGrid;
    Eigen::MatrixXi blueIndexGrid;

    Eigen::MatrixXi redIndexMap;
    Eigen::MatrixXi greenIndexMap;
    Eigen::MatrixXi blueIndexMap;

    redIndexGrid.resize(omega.cols(), omega.rows());
    greenIndexGrid.resize(omega.cols(), omega.rows());
    blueIndexGrid.resize(omega.cols(), omega.rows());

    redIndexMap.resize(redVars, 2);
    greenIndexMap.resize(greenVars, 2);
    blueIndexMap.resize(blueVars, 2);

    // cout << "Calculating Index Maps" << endl;
    calculateIndexMaps(omega, redIndexMap, greenIndexMap, blueIndexMap);
    calculateIndexMapsFast(omega, redIndexGrid, greenIndexGrid, blueIndexGrid);

    Eigen::VectorXd redX;
    Eigen::VectorXd greenX;
    Eigen::VectorXd blueX;

    redX.resize(redVars);
    greenX.resize(greenVars);
    blueX.resize(blueVars);

    redX.setZero();
    greenX.setZero();
    blueX.setZero();

    // cout << "Initializing X" << endl;
    for (int i = 0; i < redVars; i++) {
        redX(i) = omega.smartAccess(redIndexMap(i, 1), redIndexMap(i, 0), 0);
    }
    for (int i = 0; i < greenVars; i++) {
        greenX(i) = omega.smartAccess(greenIndexMap(i, 1), greenIndexMap(i, 0), 1);
    }
    for (int i = 0; i < blueVars; i++) {
        blueX(i) = omega.smartAccess(blueIndexMap(i, 1), blueIndexMap(i, 0), 2);
    }

    Eigen::SparseMatrix<double> redA;
    Eigen::SparseMatrix<double> greenA;
    Eigen::SparseMatrix<double> blueA;

    Eigen::VectorXd redB;
    Eigen::VectorXd greenB;
    Eigen::VectorXd blueB;

    redA.resize(redVars, redVars);
    greenA.resize(greenVars, greenVars);
    blueA.resize(blueVars, blueVars);

    redB.resize(redVars);
    greenB.resize(greenVars);
    blueB.resize(blueVars);

    // cout << "Calculating Hessians" << endl;
    calculateHessianAndB(redA, redB, redIndexGrid, redIndexMap, imageB, imageA, boundary, omega, edgeTh, omegaOffset, 0);
    calculateHessianAndB(greenA, greenB, greenIndexGrid, greenIndexMap, imageB, imageA, boundary, omega, edgeTh, omegaOffset, 1);
    calculateHessianAndB(blueA, blueB, blueIndexGrid, blueIndexMap, imageB, imageA, boundary, omega, edgeTh, omegaOffset, 2);

    // Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > redSolver;
    // Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > greenSolver;
    // Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > blueSolver;

    Eigen::ConjugateGradient< Eigen::SparseMatrix<double> > redSolver;
    Eigen::ConjugateGradient< Eigen::SparseMatrix<double> > greenSolver;
    Eigen::ConjugateGradient< Eigen::SparseMatrix<double> > blueSolver;

    redSolver.compute(redA);
    greenSolver.compute(greenA);
    blueSolver.compute(blueA);

    redX = redSolver.solve(-redB);
    greenX = greenSolver.solve(-greenB);
    blueX = blueSolver.solve(-blueB);

    cout << "ERROR: " << (redA * redX + redB).norm() << endl;

    for (int i = 0; i < redVars; i++) {
        results.set(redIndexMap(i, 1) + omegaOffset[0], redIndexMap(i, 0) + omegaOffset[1], 0, redX(i));
    }

    for (int i = 0; i < greenVars; i++) {
        results.set(greenIndexMap(i, 1) + omegaOffset[0], greenIndexMap(i, 0) + omegaOffset[1], 1, greenX(i));
    }

    for (int i = 0; i < blueVars; i++) {
        results.set(blueIndexMap(i, 1) + omegaOffset[0], blueIndexMap(i, 0) + omegaOffset[1], 2, blueX(i));
    }

    return results;
}

// imageA is the source image
// imageB is the dest image
Image seamlessPoissonCloningLog(const Image& imageA, const Image& imageB, const Image& omega, Eigen::Vector2i omegaOffset) {
    Image results(imageB.cols(), imageB.rows(), 3);

    for (int c = 0; c < 3; c++) {
        for (int i = 0; i < results.rows(); i++) {
            for (int j = 0; j < results.cols(); j++) {
                results.set(j, i, c, imageB.smartAccess(j, i, c));
            }
        }
    }

    Image boundary(omega.cols(), omega.rows(), 3);

    // cout << "Creating Boundary" << endl;
    boundary = calculateBoundaryMap(omega);

    int redVars;
    int greenVars;
    int blueVars;

    // cout << "Counting Vars" << endl;
    countVars(omega, redVars, greenVars, blueVars);

    Eigen::MatrixXi redIndexGrid;
    Eigen::MatrixXi greenIndexGrid;
    Eigen::MatrixXi blueIndexGrid;

    Eigen::MatrixXi redIndexMap;
    Eigen::MatrixXi greenIndexMap;
    Eigen::MatrixXi blueIndexMap;

    redIndexGrid.resize(omega.cols(), omega.rows());
    greenIndexGrid.resize(omega.cols(), omega.rows());
    blueIndexGrid.resize(omega.cols(), omega.rows());

    redIndexMap.resize(redVars, 2);
    greenIndexMap.resize(greenVars, 2);
    blueIndexMap.resize(blueVars, 2);

    // cout << "Calculating Index Maps" << endl;
    calculateIndexMaps(omega, redIndexMap, greenIndexMap, blueIndexMap);
    calculateIndexMapsFast(omega, redIndexGrid, greenIndexGrid, blueIndexGrid);

    Eigen::VectorXd redX;
    Eigen::VectorXd greenX;
    Eigen::VectorXd blueX;

    redX.resize(redVars);
    greenX.resize(greenVars);
    blueX.resize(blueVars);

    redX.setZero();
    greenX.setZero();
    blueX.setZero();

    // cout << "Initializing X" << endl;
    for (int i = 0; i < redVars; i++) {
        redX(i) = omega.smartAccess(redIndexMap(i, 1), redIndexMap(i, 0), 0);
    }
    for (int i = 0; i < greenVars; i++) {
        greenX(i) = omega.smartAccess(greenIndexMap(i, 1), greenIndexMap(i, 0), 1);
    }
    for (int i = 0; i < blueVars; i++) {
        blueX(i) = omega.smartAccess(blueIndexMap(i, 1), blueIndexMap(i, 0), 2);
    }

    Eigen::SparseMatrix<double> redA;
    Eigen::SparseMatrix<double> greenA;
    Eigen::SparseMatrix<double> blueA;

    Eigen::VectorXd redB;
    Eigen::VectorXd greenB;
    Eigen::VectorXd blueB;

    redA.resize(redVars, redVars);
    greenA.resize(greenVars, greenVars);
    blueA.resize(blueVars, blueVars);

    redB.resize(redVars);
    greenB.resize(greenVars);
    blueB.resize(blueVars);

    // cout << "Calculating Hessians" << endl;
    calculateHessianAndBLog(redA, redB, redIndexGrid, redIndexMap, imageB, imageA, boundary, omega, omegaOffset, 0);
    calculateHessianAndBLog(greenA, greenB, greenIndexGrid, greenIndexMap, imageB, imageA, boundary, omega, omegaOffset, 1);
    calculateHessianAndBLog(blueA, blueB, blueIndexGrid, blueIndexMap, imageB, imageA, boundary, omega, omegaOffset, 2);

    Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > redSolver;
    Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > greenSolver;
    Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > blueSolver;

    redSolver.compute(redA);
    greenSolver.compute(greenA);
    blueSolver.compute(blueA);

    // redX.setZero();
    // greenX.setZero();
    // blueX.setZero();

    // cout << "Solving" << endl;

    redX = redSolver.solve(-redB);
    greenX = greenSolver.solve(-greenB);
    blueX = blueSolver.solve(-blueB);

    // cout << "X: " << endl << redX << endl;

    // std::cout << "#iterations:     " << redSolver.iterations() << std::endl;
    // std::cout << "estimated error: " << redSolver.error()      << std::endl;

    // for (int c = 0; c < 3; c++) {
    //     for (int i = 0; i < results.rows(); i++) {
    //         for (int j = 0; j < results.cols(); j++) {
    //             results.set(j, i, c, imageB.smartAccess(j, i, c));
    //         }
    //     }
    // }

    // cout << "Writing Image" << endl;

    for (int i = 0; i < redVars; i++) {
        results.set(redIndexMap(i, 1) + omegaOffset[0], redIndexMap(i, 0) + omegaOffset[1], 0, redX(i));
    }

    for (int i = 0; i < greenVars; i++) {
        results.set(greenIndexMap(i, 1) + omegaOffset[0], greenIndexMap(i, 0) + omegaOffset[1], 1, greenX(i));
    }

    for (int i = 0; i < blueVars; i++) {
        results.set(blueIndexMap(i, 1) + omegaOffset[0], blueIndexMap(i, 0) + omegaOffset[1], 2, blueX(i));
    }

    return results;
}

void sparseDogDemo() {
    Image imageB(DATA_DIR "/input/water.jpg");
    Image imageA(DATA_DIR "/input/doggo.jpg");
    Image omega(DATA_DIR "/input/dogweights.png");
    Eigen::Vector2i offset;
    offset[0] = 100;
    offset[1] = 0;

    Image results = seamlessPoissonCloning(imageA, imageB, omega, offset);
    results.write(DATA_DIR "/output/doggoTestSparseOffset.png");
}

void sparseDogLogDemo() {
    Image imageB(DATA_DIR "/input/water.jpg");
    Image imageA(DATA_DIR "/input/doggo.jpg");
    Image omega(DATA_DIR "/input/dogweights.png");
    Eigen::Vector2i offset;
    offset[0] = 100;
    offset[1] = 0;

    imageA.toLog();
    imageB.toLog();

    Image results = seamlessPoissonCloning(imageA, imageB, omega, offset);
    results.toLin();
    results.write(DATA_DIR "/output/doggoTestSparseOffsetLog.png");
}

void sparseFoxDemo() {
    Image imageB(DATA_DIR "/input/water.jpg");
    Image imageA(DATA_DIR "/input/foxInWaterLowExposure.png");
    Image omega(DATA_DIR "/input/foxInWaterWeights.png");
    Eigen::Vector2i offset;
    offset[0] = 0;
    offset[1] = 0;

    Image results = seamlessPoissonCloning(imageA, imageB, omega, offset);
    results.write(DATA_DIR "/output/foxDemoLinLow.png");
}

void sparseFoxLogDemo() {
    Image imageB(DATA_DIR "/input/water.jpg");
    Image imageA(DATA_DIR "/input/foxInWaterLowExposure.png");
    Image omega(DATA_DIR "/input/foxInWaterWeights.png");
    Eigen::Vector2i offset;
    offset[0] = 0;
    offset[1] = 0;

    imageA.toLog();
    imageB.toLog();

    Image results = seamlessPoissonCloning(imageA, imageB, omega, offset);
    results.toLin();
    results.write(DATA_DIR "/output/foxDemoLogLow.png");
}

void printOmegaLocation(const Image& imageB, const Image& omega, Eigen::Vector2i offset, string filename) {
    Image image(imageB.cols(), imageB.rows(), 3);
    // Image image2(imageB.cols(), imageB.rows(), 3);

    for (int i = 0; i < imageB.rows(); i++) {
        for (int j = 0; j < imageB.cols(); j++) {
            image.set(j, i, 0, imageB.smartAccess(j, i, 0));
            image.set(j, i, 1, imageB.smartAccess(j, i, 1));
            image.set(j, i, 2, imageB.smartAccess(j, i, 2));

            // image2.set(j, i, 0, imageB.smartAccess(j, i, 0));
            // image2.set(j, i, 1, imageB.smartAccess(j, i, 1));
            // image2.set(j, i, 2, imageB.smartAccess(j, i, 2));
        }
    }

    for (int i = 0; i < omega.rows(); i++) {
        for (int j = 0; j < omega.cols(); j++) {
            if (omega.smartAccess(j, i, 0) > 0.3) {
                image.set(j + offset[0], i + offset[1], 0, 1.0);
                image.set(j + offset[0], i + offset[1], 1, 0.0);
                image.set(j + offset[0], i + offset[1], 2, 0.0);
                image.set(j, i, 2, 1.0);
                image.set(j, i, 2, 0.0);
                image.set(j, i, 2, 0.0);
            }
            if (omega.smartAccess(j, i, 1) > 0.3) {
                image.set(j + offset[0], i + offset[1], 0, 1.0);
                image.set(j + offset[0], i + offset[1], 1, 0.0);
                image.set(j + offset[0], i + offset[1], 2, 0.0);
                image.set(j, i, 2, 1.0);
            }
            if (omega.smartAccess(j, i, 2) > 0.3) {
                image.set(j + offset[0], i + offset[1], 0, 1.0);
                image.set(j + offset[0], i + offset[1], 1, 0.0);
                image.set(j + offset[0], i + offset[1], 2, 0.0);
                image.set(j, i, 2, 1.0);
                image.set(j, i, 2, 0.0);
                image.set(j, i, 2, 0.0);
            }
        }
    }

    image.write(filename);
}

void printGradient(const Image& image, string fileName, double offset) {
    Image gradientIm(image.cols(), image.rows(), 3);

    for (int i = 1; i < image.rows() - 1; i++) {
        for (int j = 1; j < image.cols() - 1; j++) {
            for (int c = 0; c < 3; c++) {
                double gradient = -1.0 * image.smartAccess(j-1, i, c) + 1.0 * image.smartAccess(j+1, i, c);
                gradient += -1.0 * image.smartAccess(j, i-1, c) + 1.0 * image.smartAccess(j, i+1, c);
                gradientIm.set(j, i, c, gradient + offset);
            }
        }
    }

    gradientIm.write(fileName);
}

Image edgeDetection(const Image& image, double threshold) {
    Image edgeIm(image.cols(), image.rows(), 3);
    Image gradX(image.cols(), image.rows(), 3);
    Image gradY(image.cols(), image.rows(), 3);

    edgeIm.setZero();
    gradX.setZero();
    gradY.setZero();

    for (int i = 1; i < image.rows() - 1; i++) {
        for (int j = 1; j < image.cols() - 1; j++) {
            for (int c = 0; c < 3; c++) {
                // double grad = -1.0 * image.smartAccess(j-1, i+1, c) + 1.0 * image.smartAccess(j+1, i+1, c);
                // grad += -1.0 * image.smartAccess(j-1, i-1, c) + 1.0 * image.smartAccess(j+1, i-1, c);
                // grad += -2.0 * image.smartAccess(j-1, i, c) + 2.0 * image.smartAccess(j+1, i, c);

                double grad = -image.smartAccess(j-1, i, c) + image.smartAccess(j+1, i, c);
                gradX.set(j, i, c, grad);
            }
        }
    }

    for (int i = 1; i < image.rows() - 1; i++) {
        for (int j = 1; j < image.cols() - 1; j++) {
            for (int c = 0; c < 3; c++) {
                // double grad = -1.0 * image.smartAccess(j+1, i-1, c) + 1.0 * image.smartAccess(j+1, i+1, c);
                // grad += -1.0 * image.smartAccess(j-1, i-1, c) + 1.0 * image.smartAccess(j-1, i+1, c);
                // grad += -2.0 * image.smartAccess(j, i-1, c) + 2.0 * image.smartAccess(j, i+1, c);

                double grad = -1.0 * image.smartAccess(j, i-1, c) + 1.0 * image.smartAccess(j, i+1, c);
                gradY.set(j, i, c, grad);
            }
        }
    }

    for (int i = 0; i < image.rows(); i++) {
        for (int j = 0; j < image.cols(); j++) {
            for (int c = 0; c < 3; c++) {
                double mag = sqrt(gradX.smartAccess(j, i, c) * gradX.smartAccess(j, i, c) + gradY.smartAccess(j, i, c) * gradY.smartAccess(j, i, c));
                if (mag > threshold)
                    edgeIm.set(j, i, c, 1.0);
            }
        }
    }

    return edgeIm;
}

Image edgeDetectionX(const Image& image, double threshold) {
    Image edgeIm(image.cols(), image.rows(), 3);

    Image gradX(image.cols(), image.rows(), 3);
    Image gradY(image.cols(), image.rows(), 3);

    for (int i = 1; i < image.rows() - 1; i++) {
        for (int j = 1; j < image.cols() - 1; j++) {
            for (int c = 0; c < 3; c++) {
                double grad = -1.0 * image.smartAccess(j-1, i+1, c) + 1.0 * image.smartAccess(j+1, i+1, c);
                grad += -1.0 * image.smartAccess(j-1, i-1, c) + 1.0 * image.smartAccess(j+1, i-1, c);
                grad += -2.0 * image.smartAccess(j-1, i, c) + 2.0 * image.smartAccess(j+1, i, c);

                grad = sqrt(grad * grad);
                if (grad > threshold) edgeIm.set(j, i, c, 1.0);
            }
        }
    }

    return edgeIm;
}

Image edgeDetectionY(const Image& image, double threshold) {
    Image edgeIm(image.cols(), image.rows(), 3);

    Image gradX(image.cols(), image.rows(), 3);
    Image gradY(image.cols(), image.rows(), 3);

    for (int i = 1; i < image.rows() - 1; i++) {
        for (int j = 1; j < image.cols() - 1; j++) {
            for (int c = 0; c < 3; c++) {
                double grad = -1.0 * image.smartAccess(j+1, i-1, c) + 1.0 * image.smartAccess(j+1, i+1, c);
                grad += -1.0 * image.smartAccess(j-1, i-1, c) + 1.0 * image.smartAccess(j-1, i+1, c);
                grad += -2.0 * image.smartAccess(j, i-1, c) + 2.0 * image.smartAccess(j, i+1, c);

                grad = sqrt(grad * grad);
                if (grad > threshold) edgeIm.set(j, i, c, 1.0);
            }
        }
    }

    return edgeIm;
}

double edgeGradY(const Image& image, int row, int col) {
    // int i = row;
    // int j = col;
    //
    // double grad = -1.0 * image.smartAccess(j+1, i-1, c) + 1.0 * image.smartAccess(j+1, i+1, c);
    // grad += -1.0 * image.smartAccess(j-1, i-1, c) + 1.0 * image.smartAccess(j-1, i+1, c);
    // grad += -2.0 * image.smartAccess(j, i-1, c) + 2.0 * image.smartAccess(j, i+1, c);
    //
    // return grad;

    return 0.0;
}

void sparseFaceSwapDemo() {
    Image imageB(DATA_DIR "/input/girlOne.jpg");
    Image imageA(DATA_DIR "/input/girlTwo.jpg");
    Image omega(DATA_DIR "/input/girlTwoWeights2.png");
    Eigen::Vector2i offset;
    offset[0] = 4;
    offset[1] = 25;

    Image swap = seamlessPoissonCloning(imageA, imageB, omega, offset);
    swap.write(DATA_DIR "/output/faceSwapDemo.png");

    Image imageA2(DATA_DIR "/output/faceSwapDemo.png");
    Image imageB2(DATA_DIR "/output/faceSwapDemo.png");
    Image omega2(DATA_DIR "/input/HealingWeights11.png");

    Eigen::Vector2i offset2;
    offset2[0] = -60;
    offset2[1] = 40;

    Image fix1 = seamlessPoissonCloning(imageA2, imageB2, omega2, offset2);
    fix1.write(DATA_DIR "/output/faceSwapHealedLeft.png");

    Image imageA3(DATA_DIR "/output/faceSwapHealedLeft.png");
    Image imageB3(DATA_DIR "/output/faceSwapHealedLeft.png");
    Image omega3(DATA_DIR "/input/HealingWeights22.png");

    Eigen::Vector2i offset3;
    offset3[0] = 50;
    offset3[1] = 40;

    Image fix2 = seamlessPoissonCloning(imageA3, imageB3, omega3, offset3);
    fix2.write(DATA_DIR "/output/faceSwapHealedRight.png");

    printOmegaLocation(imageB3, omega3, offset3, DATA_DIR "/output/faceSwapHealedRightOffset.png");
    printGradient(imageB, DATA_DIR "/output/girlGradOne.png", 0.5);
    printGradient(imageA, DATA_DIR "/output/girlGradTwo.png", 0.5);
    printGradient(fix2, DATA_DIR "/output/faceSwapResults.png", 0.5);

    edgeDetection(fix2, 0.0).write(DATA_DIR "/output/faceSwapED-0_0.png");
    edgeDetection(fix2, 0.1).write(DATA_DIR "/output/faceSwapED-0_1.png");
    edgeDetection(fix2, 0.2).write(DATA_DIR "/output/faceSwapED-0_2.png");
    edgeDetection(fix2, 0.3).write(DATA_DIR "/output/faceSwapED-0_3.png");
    edgeDetection(fix2, 0.4).write(DATA_DIR "/output/faceSwapED-0_4.png");
    edgeDetection(fix2, 0.5).write(DATA_DIR "/output/faceSwapED-0_5.png");
    edgeDetection(fix2, 0.6).write(DATA_DIR "/output/faceSwapED-0_6.png");
    edgeDetection(fix2, 0.7).write(DATA_DIR "/output/faceSwapED-0_7.png");
    edgeDetection(fix2, 0.8).write(DATA_DIR "/output/faceSwapED-0_8.png");
    edgeDetection(fix2, 0.9).write(DATA_DIR "/output/faceSwapED-0_9.png");
    // edgeDetection(fix2, 0.95).write(DATA_DIR "/output/faceSwapED-0_95.png");
    // edgeDetection(fix2, 0.5).write(DATA_DIR "/output/faceSwapED-0_5.png");
    // edgeDetection(fix2, 0.6).write(DATA_DIR "/output/faceSwapED-0_6.png");

    edgeDetectionX(fix2, 0.0).write(DATA_DIR "/output/faceSwapEDX-0_0.png");
    edgeDetectionX(fix2, 0.1).write(DATA_DIR "/output/faceSwapEDX-0_1.png");
    edgeDetectionX(fix2, 0.3).write(DATA_DIR "/output/faceSwapEDX-0_3.png");
    edgeDetectionX(fix2, 0.6).write(DATA_DIR "/output/faceSwapEDX-0_6.png");
    edgeDetectionX(fix2, 0.9).write(DATA_DIR "/output/faceSwapEDX-0_9.png");
    edgeDetectionX(fix2, 0.95).write(DATA_DIR "/output/faceSwapEDX-0_95.png");
    edgeDetectionX(fix2, 0.5).write(DATA_DIR "/output/faceSwapEDX-0_5.png");
    edgeDetectionX(fix2, 0.6).write(DATA_DIR "/output/faceSwapEDX-0_6.png");

    edgeDetectionY(fix2, 0.0).write(DATA_DIR "/output/faceSwapEDY-0_0.png");
    edgeDetectionY(fix2, 0.1).write(DATA_DIR "/output/faceSwapEDY-0_1.png");
    edgeDetectionY(fix2, 0.3).write(DATA_DIR "/output/faceSwapEDY-0_3.png");
    edgeDetectionY(fix2, 0.6).write(DATA_DIR "/output/faceSwapEDY-0_6.png");
    edgeDetectionY(fix2, 0.9).write(DATA_DIR "/output/faceSwapEDY-0_9.png");
    edgeDetectionY(fix2, 0.95).write(DATA_DIR "/output/faceSwapEDY-0_95.png");
    edgeDetectionY(fix2, 0.5).write(DATA_DIR "/output/faceSwapEDY-0_5.png");
    edgeDetectionY(fix2, 0.6).write(DATA_DIR "/output/faceSwapEDY-0_6.png");
}

void toonShaderDemo() {
    Image imageB(DATA_DIR "/output/faceSwapHealedRight.png");
    Image imageA(DATA_DIR "/output/faceSwapHealedRight.png");
    // Image imageB(DATA_DIR "/input/girlTwo.jpg");
    // Image imageA(DATA_DIR "/input/girlTwo.jpg");
    Image omega(DATA_DIR "/input/fullFaceWeights.png");

    // Image edge = edgeDetection(imageA, 0.09);

    Eigen::Vector2i offset;
    offset[0] = 0;
    offset[1] = 0;

    // edgeDetection(imageA, 0.03).write(DATA_DIR "/output/faceSwapED-0_03.png");
    // edgeDetection(imageA, 0.05).write(DATA_DIR "/output/faceSwapED-0_05.png");
    // edgeDetection(imageA, 0.1).write(DATA_DIR "/output/faceSwapED-0_1.png");
    // edgeDetection(imageA, 0.2).write(DATA_DIR "/output/faceSwapED-0_2.png");
    // edgeDetection(imageA, 0.3).write(DATA_DIR "/output/faceSwapED-0_3.png");
    // edgeDetection(imageA, 0.4).write(DATA_DIR "/output/faceSwapED-0_4.png");
    // edgeDetection(imageA, 0.5).write(DATA_DIR "/output/faceSwapED-0_5.png");
    // edgeDetection(imageA, 0.6).write(DATA_DIR "/output/faceSwapED-0_6.png");
    // edgeDetection(imageA, 0.7).write(DATA_DIR "/output/faceSwapED-0_7.png");
    // edgeDetection(imageA, 0.8).write(DATA_DIR "/output/faceSwapED-0_8.png");
    // edgeDetection(imageA, 0.9).write(DATA_DIR "/output/faceSwapED-0_9.png");

    // Image results = seamlessPoissonCloningToon(imageA, imageB, omega, edge, offset);
    // results.write(DATA_DIR "/output/faceWHOA.png");
    // results.setZero();

    Image results = seamlessPoissonCloningToon(imageA, imageB, omega, 0.1, offset);
    results.write(DATA_DIR "/output/faceSwapToonResultsValFaceOnly_1.png");
    results.setZero();

    results = seamlessPoissonCloningToon(imageA, imageB, omega, 0.05, offset);
    results.write(DATA_DIR "/output/faceSwapToonResultsValFaceOnly_05.png");
    results.setZero();

    results = seamlessPoissonCloningToon(imageA, imageB, omega, 0.15, offset);
    results.write(DATA_DIR "/output/faceSwapToonResultsValFaceOnly_15.png");
    results.setZero();

    // results = seamlessPoissonCloningToon(imageA, imageB, omega, 0.4, offset);
    // results.write(DATA_DIR "/output/faceSwapToonResultsValFaceOnly_4.png");
    // results.setZero();
    //
    // results = seamlessPoissonCloningToon(imageA, imageB, omega, 0.5, offset);
    // results.write(DATA_DIR "/output/faceSwapToonResultsValFaceOnly_5.png");
    // results.setZero();
    //
    // results = seamlessPoissonCloningToon(imageA, imageB, omega, 0.6, offset);
    // results.write(DATA_DIR "/output/faceSwapToonResultsValFaceOnly_6.png");
    // results.setZero();
    //
    // results = seamlessPoissonCloningToon(imageA, imageB, omega, 0.7, offset);
    // results.write(DATA_DIR "/output/faceSwapToonResultsValFaceOnly_7.png");
    // results.setZero();
    //
    // results = seamlessPoissonCloningToon(imageA, imageB, omega, 0.8, offset);
    // results.write(DATA_DIR "/output/faceSwapToonResultsValFaceOnly_8.png");
    // results.setZero();
    //
    // results = seamlessPoissonCloningToon(imageA, imageB, omega, 0.9, offset);
    // results.write(DATA_DIR "/output/faceSwapToonResultsValFaceOnly_9.png");
    // results.setZero();
}

void polarBearExample() {
    cout << "Polar Bear Example" << endl;
    // TODO
}

void testImageWrite() {
    Image img(DATA_DIR "/input/lucario.jpg");

    img.debugWrite();
}

void quickSparseMatrixVerify() {
    Eigen::SparseMatrix<double> sparse;

    sparse.resize(6, 6);

    cout << "SPARSE:" << endl;
    cout << sparse << endl;

    sparse.insert(1, 0) = 4;

    cout << "ROW SPARSE" << endl;
    cout << sparse << endl;
    cout << endl;

    Eigen::MatrixXd dense;

    dense.resize(6, 6);
    // dense.setZero();

    cout << "DENSE:" << endl;
    cout << dense << endl;

    dense(1, 0) = 4;

    cout << "ROW DENSE" << endl;
    cout << dense << endl;
    cout << endl;
}

Image naiveCombine(const Image& imageA, const Image& imageB, const Image& omega) {
    Image results(imageB.cols(), imageB.rows(), 3);

    for (int i = 0; i < imageB.rows(); i++) {
        for (int j = 0; j < imageB.cols(); j++) {
            for (int c = 0; c < 3; c++) {
                if (omega.smartAccess(j, i, c) > 0.1) results.set(j, i, c, imageA.smartAccess(j, i, c));
                else results.set(j, i, c, imageB.smartAccess(j, i, c));
            }
        }
    }

    return results;
}

void cellShadeTest() {
    Image imageB(DATA_DIR "/output/faceSwapHealedRight.png");
    Image imageB2(DATA_DIR "/output/faceSwapHealedRight.png");
    Image imageB3(DATA_DIR "/output/faceSwapHealedRight.png");
    Image imageB4(DATA_DIR "/output/faceSwapHealedRight.png");
    Image imageB5(DATA_DIR "/output/faceSwapHealedRight.png");
    imageB.cellShade(1.0 / 2.0);
    imageB2.cellShade(1.0 / 5.0);
    imageB3.cellShade(1.0 / 10.0);
    imageB4.cellShade(1.0 / 15.0);
    imageB5.cellShade(1.0 / 20.0);

    imageB.write(DATA_DIR "/output/cellShadeExample_2.png");
    imageB2.write(DATA_DIR "/output/cellShadeExample_5.png");
    imageB3.write(DATA_DIR "/output/cellShadeExample_10.png");
    imageB4.write(DATA_DIR "/output/cellShadeExample_15.png");
    imageB5.write(DATA_DIR "/output/cellShadeExample_20.png");
}

void cellShadeDemo() {
    Image imageB(DATA_DIR "/input/girlOne.jpg");
    Image imageA(DATA_DIR "/input/girlTwo.jpg");
    Image omega(DATA_DIR "/input/girlTwoWeights2.png");
    Eigen::Vector2i offset;
    offset[0] = 4;
    offset[1] = 25;

    imageA.cellShade(1.0 / 10.0);

    printGradient(imageA, DATA_DIR "/output/cellShadeGradient.png", 0.5);

    // imageB.cellShade(1.0 / 12.0);

    Image results = seamlessPoissonCloningToon(imageA, imageB, omega, 0.1, offset);
    results.write(DATA_DIR "/output/cellGradDemoImage.png");

    // Image swap = seamlessPoissonCloning(imageA, imageB, omega, offset);
    // swap.write(DATA_DIR "/output/cellShadeDemo.png");
}

void cellShadeDemo2() {
    Image imageB(DATA_DIR "/output/faceSwapHealedRight.png");
    Image imageA(DATA_DIR "/output/faceSwapHealedRight.png");
    // Image imageB(DATA_DIR "/input/girlTwo.jpg");
    // Image imageA(DATA_DIR "/input/girlTwo.jpg");
    Image omega(DATA_DIR "/input/fullFaceWeights.png");

    // Image edge = edgeDetection(imageA, 0.09);

    Eigen::Vector2i offset;
    offset[0] = 0;
    offset[1] = 0;

    imageA.cellShade(1.0 / 12.0);

    edgeDetection(imageA, 0.1).write(DATA_DIR "/output/cellED-0_5.png");
    // edgeDetection(imageA, 0.05).write(DATA_DIR "/output/faceSwapED-0_05.png");
    // edgeDetection(imageA, 0.1).write(DATA_DIR "/output/faceSwapED-0_1.png");
    // edgeDetection(imageA, 0.2).write(DATA_DIR "/output/faceSwapED-0_2.png");
    // edgeDetection(imageA, 0.3).write(DATA_DIR "/output/faceSwapED-0_3.png");
    // edgeDetection(imageA, 0.4).write(DATA_DIR "/output/faceSwapED-0_4.png");
    // edgeDetection(imageA, 0.5).write(DATA_DIR "/output/faceSwapED-0_5.png");
    // edgeDetection(imageA, 0.6).write(DATA_DIR "/output/faceSwapED-0_6.png");
    // edgeDetection(imageA, 0.7).write(DATA_DIR "/output/faceSwapED-0_7.png");
    // edgeDetection(imageA, 0.8).write(DATA_DIR "/output/faceSwapED-0_8.png");
    // edgeDetection(imageA, 0.9).write(DATA_DIR "/output/faceSwapED-0_9.png");

    // Image results = seamlessPoissonCloningToon(imageA, imageB, omega, edge, offset);
    // results.write(DATA_DIR "/output/faceWHOA.png");
    // results.setZero();

    Image results = seamlessPoissonCloningToon(imageA, imageB, omega, 0.1, offset);
    results.write(DATA_DIR "/output/cellGradTest_1_5.png");
    results.setZero();

    results = seamlessPoissonCloningToon(imageA, imageB, omega, 0.05, offset);
    results.write(DATA_DIR "/output/cellGradTest_05_5.png");
    results.setZero();

    results = seamlessPoissonCloningToon(imageA, imageB, omega, 0.15, offset);
    results.write(DATA_DIR "/output/cellGradTest_15_5.png");
    results.setZero();
}

void dogNaive() {
    Image imageB(DATA_DIR "/input/water.jpg");
    Image imageA(DATA_DIR "/input/doggo.jpg");
    Image omega(DATA_DIR "/input/dogweights.png");

    Image dog = naiveCombine(imageA, imageB, omega);
    dog.write(DATA_DIR "/output/naiveDog.png");

    printGradient(imageA, DATA_DIR "/output/naiveDogGradient.png", 0.5);
}

void inClassDemo() {
    Image imageB(DATA_DIR "/input/water.jpg");
    Image imageA(DATA_DIR "/input/doggo.jpg");
    Image omega(DATA_DIR "/input/dogweights2.png");
    Eigen::Vector2i offset;
    offset[0] = 100;
    offset[1] = 0;

    Image results = seamlessPoissonCloning(imageA, imageB, omega, offset);
    results.write(DATA_DIR "/output/doggoInClassDemo.png");
}

int main(int argc, char* argv[]) {
    // quickSparseMatrixVerify();
    //
    // std::cout << "PHOTO OF MY BINS" << std::endl;
    //
    // testImageWrite();
    //
    // Eigen::MatrixXd m(2,2);
    // m(0,0) = 3;
    // m(1,0) = 2.5;
    // m(0,1) = -1;
    // m(1,1) = m(1,0) + m(0,1);
    // std::cout << m << std::endl;
    //
    // inClassExample();
    // //
    // twodExample();
    //
    // polarBearExample();
    //
    // imageExample();

    // dogDemo();

    // sparseDogDemo();

    // sparseDogLogDemo();

    // sparseFoxDemo();
    // sparseFoxLogDemo();

    // sparseFaceSwapDemo();
    // toonShaderDemo();

    // cellShadeTest();
    // cellShadeDemo();
    // cellShadeDemo2();

    // dogNaive();

    inClassDemo();

    return 0;
}

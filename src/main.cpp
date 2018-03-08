#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>

using namespace std;

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
        sum -= 2.0 * (newImg(i, j) - oldImg(i - 1, j));
    } else if (smartAccess(i - 1, j, omega) == 1.0) {
        // cout << "N: " << newImg(i - 1, j) << endl;
        sum -= 2.0 * (newImg(i, j) - newImg(i - 1, j));
    } else {
        cout << "MAJOR ERROR: OUT OF EXPECTED BOUNDS NORTH" << endl;
    }

    // cout << "SUM POST N: " << sum << endl;

    // South Box
    if (smartAccess(i + 1, j, boundary) == 1.0) {
        // cout << "S: " << oldImg(i + 1, j) << endl;
        sum += 2.0 * (oldImg(i + 1, j) - newImg(i, j));
    } else if (smartAccess(i + 1, j, omega) == 1.0) {
        // cout << "S: " << newImg(i + 1, j) << endl;
        sum += 2.0 * (newImg(i + 1, j) - newImg(i, j));
    } else {
        cout << "MAJOR ERROR: OUT OF EXPECTED BOUNDS SOUTH" << endl;
    }

    // cout << "SUM POST S: " << sum << endl;

    // East Box
    if (smartAccess(i, j + 1, boundary) == 1.0) {
        // cout << "E: " << oldImg(i, j + 1) << endl;
        sum += 2.0 * (oldImg(i, j + 1) - newImg(i, j));
    } else if (smartAccess(i, j + 1, omega) == 1.0) {
        // cout << "E: " << newImg(i, j + 1) << endl;
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
        sum -= 2.0 * (newImg(i, j) - oldImg(i, j - 1));
    } else if (smartAccess(i, j - 1, omega) == 1.0) {
        // cout << "W: " << newImg(i, j - 1) << endl;
        sum -= 2.0 * (newImg(i, j) - newImg(i, j + 1));
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

    imageA << 5, 4, 1, 3, 1, 5, 3,
              8, 9, 6, 7, 3, 2, 4,
              3, 5, 1, 3, 1, 5, 6,
              2, 3, 1, 5, 1, 7, 0,
              3, 4, 2, 4, 5, 6, 9,
              1, 2, 4, 2, 6, 3, 7,
              9, 8, 5, 7, 4, 2, 3;

    // imageB:
    // 6 3 2 2 4 4 1
    // 4 6 3 1 3 5 1
    // 3 4 7 2 5 8 2
    // 7 7 5 3 3 3 2
    // 8 2 8 0 6 2 9
    // 1 4 3 8 3 3 4
    // 2 5 9 3 2 1 1

    imageB << 6, 3, 2, 2, 4, 4, 1,
              4, 6, 3, 1, 3, 5, 1,
              3, 4, 7, 2, 5, 8, 2,
              7, 7, 5, 3, 3, 3, 2,
              8, 2, 8, 0, 6, 2, 9,
              1, 4, 3, 8, 3, 3, 4,
              2, 5, 9, 3, 2, 1, 1;

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

    // cout << "Boundary: " << endl << boundaryMap << endl;
    //
    // cout << "REPLACE MAP: " << endl << replaceMap << endl;

    cout << "Gradiant sum at (2,2): " << calculateSumOfGradiant(2, 2, imageB, imageA, boundaryMap, replaceMap) << endl;

    // cout << "OMEGA (2, 3): " << replaceMap(0, 1) << endl;

    // cout << endl << replaceMap << endl;

    Eigen::VectorXd points = solveMissingPoints(imageB, imageA, boundaryMap, replaceMap);
}

int main(int argc, char* argv[]) {
    std::cout << "PHOTO OF MY BINS" << std::endl;

    Eigen::MatrixXd m(2,2);
    m(0,0) = 3;
    m(1,0) = 2.5;
    m(0,1) = -1;
    m(1,1) = m(1,0) + m(0,1);
    std::cout << m << std::endl;

    inClassExample();

    twodExample();

    return 0;
}

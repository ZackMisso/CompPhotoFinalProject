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

Eigen::MatrixXd calculateBoundaryMap(Eigen::MatrixXd map) {
    Eigen::MatrixXd boundary;
    // TODO
    return boundary;
}

Eigen::SparseMatrix<double> calculateBoundaryMap(Eigen::SparseMatrix<double> map) {
    Eigen::SparseMatrix<double> boundary;
    // TODO
    return boundary;
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

    

    // imageB:
    // 6 3 2 2 4 4 1
    // 4 6 3 1 3 5 1
    // 3 4 7 2 5 8 2
    // 7 7 5 3 3 3 2
    // 8 2 8 0 6 2 9
    // 1 4 3 8 3 3 4
    // 2 5 9 3 2 1 1

    // replace map:
    // 0 0 0 0 0 0 0
    // 0 0 0 0 0 0 0
    // 0 0 1 1 1 0 0
    // 0 0 1 1 1 0 0
    // 0 0 1 1 1 0 0
    // 0 0 0 0 0 0 0
    // 0 0 0 0 0 0 0

    // boundary map:
    // 0 0 0 0 0 0 0
    // 0 1 1 1 1 1 0
    // 0 1 0 0 0 1 0
    // 0 1 0 0 0 1 0
    // 0 1 0 0 0 1 0
    // 0 1 1 1 1 1 0
    // 0 0 0 0 0 0 0
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

    return 0;
}

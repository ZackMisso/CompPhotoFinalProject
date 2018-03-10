#pragma once

#include <poly/common.h>

class Image {
public:
    Image();
    Image(int w, int h, int d);
    Image(const std::string& filename);
    Image(const Image& other);

    double smartAccess(int x, int y, int z) const;
    double smartAccess(int x, int y, Eigen::Vector2i offset, int z) const;
    void set(int x, int y, int z, double value);

    bool read(const std::string& filename);
    bool write(const std::string& filename) const;
    bool debugWrite() const;

    void resize(int cols, int rows, int ch);
    void setZero();

    void toLog();
    void toLin();
    double minNonZero();

    void printRed();
    void printGreen();
    void printBlue();

    int getWidth() const;
    int getHeight() const;
    int getDepth() const;
    int cols() const;
    int rows() const;
protected:
    Eigen::MatrixXd red;
    Eigen::MatrixXd green;
    Eigen::MatrixXd blue;

    int width;
    int height;
    int depth;
};

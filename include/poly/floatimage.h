#pragma once

#include <poly/common.h>

class Image {
public:
    Image();
    Image(int w, int h, int d);
    Image(const std::string& filename);
    Image(const Image& other);

    double smartAccess(int x, int y, int z) const;

    bool read(const std::string& filename);
    bool write(const std::string& filename) const;
    bool debugWrite() const;

    void resize(int cols, int rows, int ch);

    int getWidth() const;
    int getHeight() const;
    int getDepth() const;
protected:
    Eigen::MatrixXd red;
    Eigen::MatrixXd green;
    Eigen::MatrixXd blue;

    int width;
    int height;
    int depth;
};

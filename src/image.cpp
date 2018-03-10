//
//     This class's IO was taken from code written by Wojciech Jarosz
//

#include <poly/image.h>


#define STB_IMAGE_IMPLEMENTATION

#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "stb_image_write.h"    // for stbi_write_bmp, stbi_write_hdr, stbi...

unsigned char floatToByte(float in)
{
	return int(255.0f * clamp(in, 0.0f, 1.0f));
}

float byteToFloat(const unsigned char in)
{
	return in / 255.0f;
}

std::string getExtension(const std::string &filename)
{
	if (filename.find_last_of(".") != std::string::npos)
		return filename.substr(filename.find_last_of(".") + 1);
	return "";
}

Image::Image() {
    width = 0;
    height = 0;
    depth = 0;
}

Image::Image(int w, int h, int c) {
    width = w;
    height = h;
    depth = 3;

    red = Eigen::MatrixXd(height, width);
    green = Eigen::MatrixXd(height, width);
    blue = Eigen::MatrixXd(height, width);

    red.setZero();
    green.setZero();
    blue.setZero();
}

Image::Image(const std::string& filename) {
    // cout << "READING" << endl;
    read(filename);
    // cout << "FINISHED READING" << endl;
    depth = 3;
}

Image::Image(const Image& other) {
    // TODO
    cout << "NOT YET IMPLEMENTED" << endl;
}

double Image::smartAccess(int x, int y, int z) const {
    if (z >= 3 || z < 0) return 0.0;
    if (x >= width || x < 0) return 0.0;
    if (y >= height || y < 0) return 0.0;

    if (z == 0) return red(y, x);
    if (z == 1) return green(y, x);
    if (z == 2) return blue(y, x);

    return 0.0;
}

double Image::smartAccess(int x, int y, Eigen::Vector2i offset, int z) const {
    int xx = x + offset[0];
    int yy = y + offset[1];

    if (z >= 3 || z < 0) cout << "WHAT" << endl;
    if (xx >= width || xx < 0) cout << "WHAT" << endl;
    if (yy >= height || yy < 0) cout << "WHAT" << endl;

    if (z >= 3 || z < 0) return 0.0;
    if (xx >= width || xx < 0) return 0.0;
    if (yy >= height || yy < 0) return 0.0;

    if (z == 0) return red(yy, xx);
    if (z == 1) return green(yy, xx);
    if (z == 2) return blue(yy, xx);

    return 0.0;
}

void Image::set(int x, int y, int z, double value) {
    if (z >= 3 || z < 0) return;
    if (x >= width || x < 0) return;
    if (y >= height || y < 0) return;

    if (z == 0) red(y, x) = value;
    if (z == 1) green(y, x) = value;
    if (z == 2) blue(y, x) = value;
}

void Image::setZero() {
    red.setZero();
    green.setZero();
    blue.setZero();
}

void Image::printRed() {
    cout << red << endl;
}

void Image::printBlue() {
    cout << blue << endl;
}

void Image::printGreen() {
    cout << green << endl;
}

void Image::toLog() {
    double minf = minNonZero();

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            red(i, j) = log10(max(red(i, j), minf));
            green(i, j) = log10(max(green(i, j), minf));
            blue(i, j) = log10(max(blue(i, j), minf));
        }
    }
}

void Image::toLin() {
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            red(i, j) = pow(10.0, red(i, j));
            green(i, j) = pow(10.0, green(i, j));
            blue(i, j) = pow(10.0, blue(i, j));
        }
    }
}

double Image::minNonZero() {
    double minNon = 2.0;

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            if (red(i, j) < minNon && red(i, j) > 0.0) minNon = red(i, j);
            if (green(i, j) < minNon && green(i, j) > 0.0) minNon = green(i, j);
            if (blue(i, j) < minNon && blue(i, j) > 0.0) minNon = blue(i, j);
        }
    }

    return minNon;
}

bool Image::read(const std::string& filename) {
    int n, w, h;

	try
	{
		if (stbi_is_hdr(filename.c_str()))
		{
			float *floatPixels = stbi_loadf(filename.c_str(), &w, &h, &n, 3);
			if (floatPixels)
			{
				resize(w, h, 3);

				for (int x = 0; x < w; x++) {
					for (int y = 0; y < h; y++) {
                        red(y, x) = floatPixels[3 * (x + y * w) + 0];
                        green(y, x) = floatPixels[3 * (x + y * w) + 1];
                        blue(y, x) = floatPixels[3 * (x + y * w) + 2];
                    }
                }

				stbi_image_free(floatPixels);
				return true;
			}
			else
				throw runtime_error("Could not load HDR image.");
		}
		else
		{
			unsigned char *bytePixels = stbi_load(filename.c_str(), &w, &h, &n, 3);
			if (bytePixels)
			{
				resize(w, h, 3);

				for (int x = 0; x < w; x++) {
					for (int y = 0; y < h; y++) {
                        red(y, x) = byteToFloat(bytePixels[3 * (x + y * w) + 0]);
                        green(y, x) = byteToFloat(bytePixels[3 * (x + y * w) + 1]);
                        blue(y, x) = byteToFloat(bytePixels[3 * (x + y * w) + 2]);
                    }
                }

				stbi_image_free(bytePixels);
				return true;
			}
			else
				throw runtime_error("Could not load LDR image.");
		}
	}
	catch (const exception &e)
	{
		cerr << "Image decoder error in FloatImage::read(...) for file: \"" << filename << "\":\n\t"
		     << stbi_failure_reason() << endl;
		return false;
	}
}

bool Image::write(const std::string& filename) const {
    if (getDepth() != 1 && getDepth() != 3 && getDepth() != 4) {
        cout << "CHANNEL ERROR" << endl;
        return false;
    }

    // cout << endl << "RED:" << endl << red << endl;

	string extension = getExtension(filename);
    // cout << "EXTENSION: " << extension << endl;
    // cout << "WIDTH: " << width << endl;
    // cout << "HEIGHT: " << height << endl;
	transform(extension.begin(),
	          extension.end(),
	          extension.begin(),
	          ::tolower);
	try
	{
		if (extension == "hdr")
		{
			// stbi_write_hdr expects color channels for a single pixel to be adjacent in memory
			vector<float> floatPixels(height * width * 3, 1.0f);
			for (int x = 0; x < width; x++)
				for (int y = 0; y < height; y++)
					for (int c = 0; c < getDepth(); c++)
						floatPixels[c + x * 3 + y * 3 * width] = smartAccess(x, y, c);

			if (!stbi_write_hdr(filename.c_str(), width, height, getDepth(), &floatPixels[0]))
				throw runtime_error("Could not write HDR image.");
		}
		else if (extension == "png" ||
			extension == "bmp" ||
			extension == "tga" ||
			extension == "jpg" || extension == "jpeg")
		{
			int outputChannels = 4;
			vector<unsigned char> bytePixels(height * width * outputChannels, 255);
			int c;
			for (int x = 0; x < width; x++)
				for (int y = 0; y < height; y++)
				{
					for (c = 0; c < getDepth(); c++) {
						bytePixels[c + x * outputChannels + y * outputChannels * width] =
							floatToByte(smartAccess(x, y, c));
                            // cout << smartAccess(x, y, c) << endl;
                        }

					for (; c < 3; c++)
						// Only executes when there is one channel
						bytePixels[c + x * outputChannels + y * outputChannels * width] =
							floatToByte(smartAccess(x, y, 0));
				}

			if (extension == "png")
			{
				if (!stbi_write_png(filename.c_str(), width, height,
				                    outputChannels, &bytePixels[0],
				                    sizeof(unsigned char) * width * outputChannels))
					throw runtime_error("Could not write PNG image.");
			}
			else if (extension == "bmp")
			{
				if (!stbi_write_bmp(filename.c_str(), width, height, outputChannels, &bytePixels[0]))
					throw runtime_error("Could not write BMP image.");
			}
			else if (extension == "tga")
			{
				if (!stbi_write_tga(filename.c_str(), width, height, outputChannels, &bytePixels[0]))
					throw runtime_error("Could not write TGA image.");
			}
			else if (extension == "jpg" || extension == "jpeg")
			{
				if (!stbi_write_jpg(filename.c_str(), width, height, outputChannels, &bytePixels[0], 100))
					throw runtime_error("Could not write JPG image.");
			}
		}
		else
			throw invalid_argument("Could not determine desired file type from extension.");
	}
	catch (const exception &e)
	{
		// if there's an error, display it
		cerr << "Error in FloatImage::write(...) for file:  \"" << filename << "\":\n\t" << e.what() << endl;
	}

	return true;
}

void Image::resize(int cols, int rows, int ch) {
    red.resize(rows, cols);
    green.resize(rows, cols);
    blue.resize(rows, cols);

    red.setZero();
    green.setZero();
    blue.setZero();

    width = cols;
    height = rows;
    depth = 3;
}

bool Image::debugWrite() const {
    ostringstream ss;
	ss << DATA_DIR "/output/" << "TEST" << ".png";
	string filename = ss.str();
    cout << "WRITING TO: " << filename << endl;
	return write(filename);
}

int Image::getWidth() const {
    return width;
}

int Image::getHeight() const {
    return height;
}

int Image::getDepth() const {
    return depth;
}

int Image::cols() const {
    return width;
}

int Image::rows() const {
    return height;
}

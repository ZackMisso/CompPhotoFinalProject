// /*!
//     \file floatimage.cpp
//     \brief Contains the implementation of a floating-point image class with an arbitrary number of channels
//     \author Wojciech Jarosz
//
//     CS 89/189 Computational Aspects of Digital Photography C++ basecode.
// */
//
// #include "floatimage.h"
// #include "utils.h"
// #include <math.h>
// #include <iostream>
// #include <sstream>
// #include <algorithm>
//
// #define STB_IMAGE_IMPLEMENTATION
// #include "stb_image.h"
//
// #define STB_IMAGE_WRITE_IMPLEMENTATION
// #include "stb_image_write.h"     // for stbi_write_bmp, stbi_write_hdr, stbi...
//
// using namespace std;
//
// // anonymous (unnamed) namespaces are used in C++ to make functions and variable
// // local to the file and not pollute the global namespace. This is the same as
// // the C way of having a static global variable or static function
// namespace
// {
//
// string getExtension(const string &filename)
// {
// 	if (filename.find_last_of(".") != string::npos)
// 		return filename.substr(filename.find_last_of(".") + 1);
// 	return "";
// }
//
// float byteToFloat(const unsigned char in)
// {
// #ifndef ANY_ASSIGNMENT
// 	// alternate conversion based on Andrew Kensler's blog
// //    return (in+0.5f)/(256.0f);
// #endif
// 	return in / 255.0f;
// }
//
// unsigned char floatToByte(float in)
// {
// #ifndef ANY_ASSIGNMENT
// 	// alternate conversion based on Andrew Kensler's blog
// //    return clamp(int(256.0f*in), 0, 255);
// #endif
// 	return int(255.0f * clamp(in, 0.0f, 1.0f));
// }
//
// #ifndef ASSIGNMENT_1
// void compareXYDimensions(const FloatImage &im1, const FloatImage &im2)
// {
// 	if (im1.sizeX() != im2.sizeX() ||
// 		im1.sizeY() != im2.sizeY())
// 		throw MismatchedDimensionsException();
// }
// #endif
//
// } // end namespace
//
//
// FloatImage::FloatImage() : Array3D<float>()
// {
// 	// empty
// }
//
// FloatImage::FloatImage(int width, int height, int channels) :
// 	Array3D<float>(width, height, channels)
// {
// 	// empty
// }
//
// FloatImage::FloatImage(const string &filename) : Array3D<float>()
// {
// 	read(filename);
// }
//
// FloatImage::FloatImage(const FloatImage &in) : Array3D<float>()
// {
// 	resize(in.width(), in.height(), in.depth());
// 	m_data = in.m_data;
// }
//
// FloatImage &FloatImage::operator=(const FloatImage &in)
// {
// 	m_data = in.m_data;
// 	m_sizeX = in.m_sizeX;
// 	m_sizeY = in.m_sizeY;
// 	m_sizeZ = in.m_sizeZ;
// 	m_strideZ = in.m_strideZ;
// 	return *this;
// }
//
// #ifdef ASSIGNMENT_4
// #ifndef ASSIGNMENT_3
// // A4: Safe Accessor that will return a black pixel (clampToEdge = false) or the
// // nearest pixel value (clampToEdge = true) when indexing out of the bounds of the image
// float FloatImage::smartAccessor(int x, int y, int z, bool clampToEdge) const
// {
// 	return 0.0f; // CHANGEME
// }
// #endif
// #else
// // A4: Safe Accessor that will return a black pixel (clampToEdge = false) or the
// // nearest pixel value (clampToEdge = true) when indexing out of the bounds of the image
// float FloatImage::smartAccessor(int x, int y, int z, bool clampToEdge) const
// {
// 	if (clampToEdge)
// 	{
// 		x = clamp(x, 0, width() - 1);
// 		y = clamp(y, 0, height() - 1);
// 	}
// 	else
// 	{
// 		if (x < 0 || x >= width() || y < 0 || y >= height())
// 			return 0.0f;
// 	}
//
// 	return smartAccess(x, y, z);
// }
// #endif
//
// void FloatImage::clear(const vector<float> &channelValues)
// {
// 	// check z bounds
//
// 	for (int z = 0; z < sizeZ(); ++z)
// 		for (int y = 0; y < sizeY(); ++y)
// 			for (int x = 0; x < sizeX(); ++x)
// 				smartAccess(x, y, z) = channelValues[z];
// }
//
// #ifndef ASSIGNMENT_1
// FloatImage &FloatImage::operator+=(float s)
// {
// 	for (int i = 0; i < size(); ++i)
// 		smartAccess(i) += s;
// 	return *this;
// }
//
// FloatImage &FloatImage::operator*=(float s)
// {
// 	for (int i = 0; i < size(); ++i)
// 		smartAccess(i) *= s;
// 	return *this;
// }
//
// FloatImage &FloatImage::operator/=(float s)
// {
// 	if (s == 0.0f)
// 		throw DivideByZeroException();
// 	return *this *= 1.0f / s;
// }
//
// FloatImage &FloatImage::operator+=(const FloatImage &other)
// {
// 	compareXYDimensions(*this, other);
//
// 	if (this->channels() == other.channels())
// 	{
// 		for (int i = 0; i < size(); ++i)
// 			smartAccess(i) += other(i);
// 	}
// 	else
// 	{
// 		if (other.channels() != 1)
// 			throw MismatchedDimensionsException();
//
// 		// if other is a grayscale image, add it to each channel
// 		for (int c = 0; c < channels(); c++)
// 			for (int y = 0; y < height(); y++)
// 				for (int x = 0; x < width(); x++)
// 					smartAccess(x, y, c) += other(x, y, 0);
// 	}
// 	return *this;
// }
//
// FloatImage &FloatImage::operator-=(const FloatImage &other)
// {
// 	compareXYDimensions(*this, other);
//
// 	if (this->channels() == other.channels())
// 	{
// 		for (int i = 0; i < size(); ++i)
// 			smartAccess(i) -= other(i);
// 	}
// 	else
// 	{
// 		if (other.channels() != 1)
// 			throw MismatchedDimensionsException();
//
// 		// if other is a grayscale image, subtract it from each channel
// 		for (int c = 0; c < channels(); c++)
// 			for (int y = 0; y < height(); y++)
// 				for (int x = 0; x < width(); x++)
// 					smartAccess(x, y, c) -= other(x, y, 0);
// 	}
// 	return *this;
// }
//
// FloatImage &FloatImage::operator*=(const FloatImage &other)
// {
// 	compareXYDimensions(*this, other);
//
// 	if (this->channels() == other.channels())
// 	{
// 		for (int i = 0; i < size(); ++i)
// 			smartAccess(i) *= other(i);
// 	}
// 	else
// 	{
// 		if (other.channels() != 1)
// 			throw MismatchedDimensionsException();
//
// 		// if other is a grayscale image, multiply each channel by it
// 		for (int c = 0; c < channels(); c++)
// 			for (int y = 0; y < height(); y++)
// 				for (int x = 0; x < width(); x++)
// 					smartAccess(x, y, c) *= other(x, y, 0);
// 	}
// 	return *this;
// }
//
// FloatImage &FloatImage::operator/=(const FloatImage &other)
// {
// 	compareXYDimensions(*this, other);
//
// 	if (this->channels() == other.channels())
// 	{
// 		for (int i = 0; i < size(); ++i)
// 		{
// 			if (other(i) == 0.0f)
// 				throw DivideByZeroException();
// 			smartAccess(i) /= other(i);
// 		}
// 	}
// 	else
// 	{
// 		if (other.channels() != 1)
// 			throw MismatchedDimensionsException();
//
// 		// if other is a grayscale image, divide each channel by it
// 		for (int c = 0; c < channels(); c++)
// 			for (int y = 0; y < height(); y++)
// 				for (int x = 0; x < width(); x++)
// 				{
// 					if (other(x, y, 0) == 0.0f)
// 						throw DivideByZeroException();
// 					smartAccess(x, y, c) /= other(x, y, 0);
// 				}
// 	}
// 	return *this;
// }
//
// FloatImage operator*(float s, const FloatImage &other)
// {
// 	FloatImage ret(other);
// 	return ret *= s;
// }
//
// FloatImage operator/(float s, const FloatImage &other)
// {
// 	FloatImage ret(other.width(), other.height(), other.channels());
// 	for (int i = 0; i < ret.size(); ++i)
// 	{
// 		if (other(i) == 0.0f)
// 			throw DivideByZeroException();
// 		ret(i) = s / other(i);
// 	}
// 	return ret;
// }
//
// FloatImage operator+(float s, const FloatImage &other)
// {
// 	FloatImage ret(other);
// 	return ret += s;
// }
//
// FloatImage operator-(float s, const FloatImage &other)
// {
// 	FloatImage ret(other);
// 	for (int i = 0; i < ret.size(); ++i)
// 		ret(i) = s - ret(i);
// 	return ret;
// }
//
// FloatImage FloatImage::operator+(float s) const
// {
// 	FloatImage ret(*this);
// 	return ret += s;
// }
//
// FloatImage FloatImage::operator*(float s) const
// {
// 	FloatImage ret(*this);
// 	return ret *= s;
// }
//
// FloatImage FloatImage::operator+(const FloatImage &other) const
// {
// 	FloatImage ret(*this);
// 	return ret += other;
// }
//
// FloatImage FloatImage::operator-(const FloatImage &other) const
// {
// 	FloatImage ret(*this);
// 	return ret -= other;
// }
//
// FloatImage FloatImage::operator*(const FloatImage &other) const
// {
// 	FloatImage ret(*this);
// 	return ret *= other;
// }
//
// FloatImage FloatImage::operator/(const FloatImage &other) const
// {
// 	FloatImage ret(*this);
// 	return ret /= other;
// }
//
// float FloatImage::min(int c) const
// {
// 	float mn = smartAccess(0, 0, c);
//
// 	for (int y = 0; y < height(); y++)
// 		for (int x = 0; x < width(); x++)
// 			mn = std::min(mn, smartAccess(x, y, c));
//
// 	return mn;
// }
//
// float FloatImage::min() const
// {
// 	float mn = 0;
// 	for (int c = 0; c < channels(); ++c)
// 		mn = std::min(mn, min(c));
// 	return mn;
// }
//
// float FloatImage::max(int c) const
// {
// 	float mx = smartAccess(0, 0, c);
//
// 	for (int y = 0; y < height(); y++)
// 		for (int x = 0; x < width(); x++)
// 			mx = std::max(mx, smartAccess(x, y, c));
//
// 	return mx;
// }
//
// float FloatImage::max() const
// {
// 	float mx = 0;
// 	for (int c = 0; c < channels(); ++c)
// 		mx = std::max(mx, max(c));
// 	return mx;
// }
// #endif
//
// bool FloatImage::read(const string &filename)
// {
// 	int n, w, h;
//
// 	try
// 	{
// 		if (stbi_is_hdr(filename.c_str()))
// 		{
// 			float *floatPixels = stbi_loadf(filename.c_str(), &w, &h, &n, 3);
// 			if (floatPixels)
// 			{
// 				resize(w, h, 3);
//
// 				for (int x = 0; x < w; x++)
// 					for (int y = 0; y < h; y++)
// 						for (int c = 0; c < 3; c++)
// 							smartAccess(x, y, c) = floatPixels[3 * (x + y * w) + c];
//
// 				stbi_image_free(floatPixels);
// 				return true;
// 			}
// 			else
// 				throw runtime_error("Could not load HDR image.");
// 		}
// 		else
// 		{
// 			unsigned char *bytePixels = stbi_load(filename.c_str(), &w, &h, &n, 3);
// 			if (bytePixels)
// 			{
// 				resize(w, h, 3);
//
// 				for (int x = 0; x < w; x++)
// 					for (int y = 0; y < h; y++)
// 						for (int c = 0; c < 3; c++)
// 							smartAccess(x, y, c) = byteToFloat(bytePixels[3 * (x + y * w) + c]);
//
// 				stbi_image_free(bytePixels);
// 				return true;
// 			}
// 			else
// 				throw runtime_error("Could not load LDR image.");
// 		}
// 	}
// 	catch (const exception &e)
// 	{
// 		cerr << "Image decoder error in FloatImage::read(...) for file: \"" << filename << "\":\n\t"
// 		     << stbi_failure_reason() << endl;
// 		return false;
// 	}
// }
//
// bool FloatImage::write(const string &filename) const
// {
// 	if (channels() != 1 && channels() != 3 && channels() != 4)
// 		throw ChannelException();
//
// 	string extension = getExtension(filename);
// 	transform(extension.begin(),
// 	          extension.end(),
// 	          extension.begin(),
// 	          ::tolower);
// 	try
// 	{
// 		if (extension == "hdr")
// 		{
// 			// stbi_write_hdr expects color channels for a single pixel to be adjacent in memory
// 			vector<float> floatPixels(height() * width() * 3, 1.0f);
// 			for (int x = 0; x < width(); x++)
// 				for (int y = 0; y < height(); y++)
// 					for (int c = 0; c < channels(); c++)
// 						floatPixels[c + x * 3 + y * 3 * width()] = smartAccess(x, y, c);
//
// 			if (!stbi_write_hdr(filename.c_str(), width(), height(), channels(), &floatPixels[0]))
// 				throw runtime_error("Could not write HDR image.");
// 		}
// 		else if (extension == "png" ||
// 			extension == "bmp" ||
// 			extension == "tga" ||
// 			extension == "jpg" || extension == "jpeg")
// 		{
// 			int outputChannels = 4;
// 			vector<unsigned char> bytePixels(height() * width() * outputChannels, 255);
// 			int c;
// 			for (int x = 0; x < width(); x++)
// 				for (int y = 0; y < height(); y++)
// 				{
// 					for (c = 0; c < channels(); c++)
// 						bytePixels[c + x * outputChannels + y * outputChannels * width()] =
// 							floatToByte(smartAccess(x, y, c));
//
// 					for (; c < 3; c++)
// 						// Only executes when there is one channel
// 						bytePixels[c + x * outputChannels + y * outputChannels * width()] =
// 							floatToByte(smartAccess(x, y, 0));
// 				}
//
// 			if (extension == "png")
// 			{
// 				if (!stbi_write_png(filename.c_str(), width(), height(),
// 				                    outputChannels, &bytePixels[0],
// 				                    sizeof(unsigned char) * width() * outputChannels))
// 					throw runtime_error("Could not write PNG image.");
// 			}
// 			else if (extension == "bmp")
// 			{
// 				if (!stbi_write_bmp(filename.c_str(), width(), height(), outputChannels, &bytePixels[0]))
// 					throw runtime_error("Could not write BMP image.");
// 			}
// 			else if (extension == "tga")
// 			{
// 				if (!stbi_write_tga(filename.c_str(), width(), height(), outputChannels, &bytePixels[0]))
// 					throw runtime_error("Could not write TGA image.");
// 			}
// 			else if (extension == "jpg" || extension == "jpeg")
// 			{
// 				if (!stbi_write_jpg(filename.c_str(), width(), height(), outputChannels, &bytePixels[0], 100))
// 					throw runtime_error("Could not write JPG image.");
// 			}
// 		}
// 		else
// 			throw invalid_argument("Could not determine desired file type from extension.");
// 	}
// 	catch (const exception &e)
// 	{
// 		// if there's an error, display it
// 		cerr << "Error in FloatImage::write(...) for file:  \"" << filename << "\":\n\t" << e.what() << endl;
// 	}
//
// 	return true;
// }
//
// int FloatImage::s_debugWriteNumber = 0;
//
// bool FloatImage::debugWrite() const
// {
// 	ostringstream ss;
// 	ss << DATA_DIR "/output/" << s_debugWriteNumber << ".png";
// 	string filename = ss.str();
// 	s_debugWriteNumber++;
// 	return write(filename);
// }

#include <poly/floatimage.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"     // for stbi_write_bmp, stbi_write_hdr, stbi...

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
}

Image::Image(const std::string& filename) {
    // TODO
    cout << "NOT YET IMPLEMENTED" << endl;
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

	string extension = getExtension(filename);
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
					for (c = 0; c < getDepth(); c++)
						bytePixels[c + x * outputChannels + y * outputChannels * width] =
							floatToByte(smartAccess(x, y, c));

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
}

bool Image::debugWrite() const {
    ostringstream ss;
	ss << DATA_DIR "/output/" << "TEST" << ".png";
	string filename = ss.str();
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

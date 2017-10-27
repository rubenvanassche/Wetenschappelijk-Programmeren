//
// Created by Ruben Van Assche on 24/10/17.
//
#include <string>
#include <fstream>
#include <list>

#ifndef DATA_FITTING_POINTSWRITER_H
#define DATA_FITTING_POINTSWRITER_H

/**
 * Writes a collection of points to a file
 */
class PointsWriter {
public:
    PointsWriter(std::string filename);

    void header(std::string one, std::string two);
    void line(double one, double two);

    virtual ~PointsWriter();

private:
    std::ofstream stream;
    void writeHeader(std::list<std::string>& items);
    void writeLine(std::list<double>& items);
};


#endif //DATA_FITTING_POINTSWRITER_H

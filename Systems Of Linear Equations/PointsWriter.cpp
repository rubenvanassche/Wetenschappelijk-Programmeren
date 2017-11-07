//
// Created by Ruben Van Assche on 24/10/17.
//

#include "PointsWriter.h"

template <typename Iter>
Iter next(Iter iter)
{
    return ++iter;
}

PointsWriter::PointsWriter(std::string filename) : stream(filename, std::ofstream::trunc) {
    this->stream.precision(20);
}

void PointsWriter::header(std::string one, std::string two) {
    std::list<std::string> items;
    items.push_back(one);
    items.push_back(two);

    this->writeHeader(items);
}

void PointsWriter::line(double one, double two) {
    std::list<double> items;
    items.push_back(one);
    items.push_back(two);

    this->writeLine(items);
}

void PointsWriter::writeHeader(std::list<std::string>& items) {
    for(auto it = items.begin();it != items.end();it++){
        this->stream << *it;
        if(std::next(it) != items.end()){
            this->stream << " ";
        }else{
            this->stream << "\n";
        }
    }
}

void PointsWriter::writeLine(std::list<double>& items) {
    for(auto it = items.begin();it != items.end();it++){
        this->stream << *it;
        if(std::next(it) != items.end()){
            this->stream << " ";
        }else{
            this->stream << "\n";
        }
    }
}

PointsWriter::~PointsWriter() {
    this->stream.close();
}

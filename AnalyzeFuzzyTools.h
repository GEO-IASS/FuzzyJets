#ifndef ANALYZEFUZZYTOOLS_H // header guard
#define ANALYZEFUZZYTOOLS_H

#include <vector>
#include <string>
#include <stdint.h>

#include <TColor.h>

namespace StyleTypes {
    enum HistOptions {
        LINE = 1,
        FILL = 2,
        DASHED = 4,
        STRIPED = 8
    };
}

class CanvasHelper {
public:
    std::string xLabel;
    std::string yLabel;
    std::string title;

    size_t width;
    size_t height;

    CanvasHelper(std::string xLabel, std::string yLabel, std::string title,
                 size_t width, size_t height) {
        this->xLabel = xLabel;
        this->yLabel = yLabel;
        this->title = title;
        this->width = width;
        this->height = height;
    }
};

class HistHelper {
public:
    std::string fileName;
    std::string branchName;
    std::string legendLabel;

    size_t      ticks;

    double      minEdge, maxEdge;
    size_t      nBins;

    Color_t color;

    StyleTypes::HistOptions options;

    HistHelper(std::string fileName, std::string branchName,
               std::string legendLabel, size_t ticks, double minEdge,
               double maxEdge, size_t nBins, StyleTypes::HistOptions options,
               Color_t color) {
        this->fileName = fileName;
        this->branchName = branchName;
        this->legendLabel = legendLabel;
        this->ticks = ticks;
        this->options = options;
        this->minEdge = minEdge;
        this->maxEdge = maxEdge;
        this->nBins = nBins;
        this->color = color;
    }
};


template<typename T>
std::vector<T>
loadSingleBranch(std::string const& file,
                 std::string const& branch);

template<typename T>
std::vector<std::vector<T> >
loadBranches(std::vector<std::string> const& files,
             std::vector<std::string> const& branches);

template<typename T>
void prettyHist(std::vector<HistHelper> const& histdecs,
                CanvasHelper const& cdec);

#endif

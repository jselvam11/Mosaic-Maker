/**
 * @file maptiles.cpp
 * Code for the maptiles function.
 */

#include <iostream>
#include <map>
#include "maptiles.h"
//#include "cs225/RGB_HSL.h"

using namespace std;


Point<3> convertToXYZ(LUVAPixel pixel) {
    return Point<3>( pixel.l, pixel.u, pixel.v );
}

MosaicCanvas* mapTiles(SourceImage const& theSource,
                       vector<TileImage>& theTiles)
{
    /**
     * @todo Implement this function!
     */
    MosaicCanvas *finalImage = new MosaicCanvas(theSource.getRows(), theSource.getColumns());

    vector<Point<3>> tileColors;


    for (TileImage i : theTiles) {
        double tileC [3] = {i.getAverageColor().l, i.getAverageColor().u, i.getAverageColor().v};
        Point<3> p = Point<3>(tileC);
        tileColors.push_back(p);

    }

    KDTree<3> *tileTree = new KDTree<3>(tileColors);

    map<Point<3>, unsigned> tilesMap;

    for (unsigned i = 0; i < tileColors.size(); i++) {
        tilesMap.insert(pair<Point<3>, unsigned>(tileColors[i], i));
    }

    for (int i = 0; i < finalImage->getRows(); i++) {
        for (int j = 0; j < finalImage->getColumns(); j++) {
            Point<3> regionColor = Point<3>(theSource.getRegionColor(i,j).l, theSource.getRegionColor(i,j).u, theSource.getRegionColor(i,j).v);
            Point<3> neighborColor = tileTree->findNearestNeighbor(regionColor);
            finalImage->setTile(i, j, &theTiles[tilesMap[neighborColor]]);
        }
    }

    return finalImage;
}


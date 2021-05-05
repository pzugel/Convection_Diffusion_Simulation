#include <iostream>
#include <fstream>	
#include <sstream>
#include <set>
#include <vector>

class VTKReader { 
	public:
		double minXValue = 0;
		double minYValue = 0;
		double maxXValue = 0;
		double maxYValue = 0;
		int numXCoords = 0;
		int numYCoords = 0;
		int numXCoordsTotal = 0;
		int numYCoordsTotal = 0;
		int numValues = 0;
		int numPoints = 0;
		int numCells = 0;

		double xOffset = 0;
		double yOffset = 0;

		std::set<double> xCoordsSet = {};
		std::set<double> yCoordsSet = {};

		std::vector<double> xCoords = {};
		std::vector<double> yCoords = {};
		std::vector<double> values = {};
		
		bool init(const char* filename);
		void get_VTKData(double (&returnDataArray)[100][100]);

	protected:
		std::string VTUOrderedPoints = "";
};

bool VTKReader::init(const char* filename){

	//Reinitialize to default value
	minXValue = 0;
	minYValue = 0;
	maxXValue = 0;
	maxYValue = 0;
	numXCoords = 0;
	numYCoords = 0;
	numXCoordsTotal = 0;
	numYCoordsTotal = 0;
	numValues = 0;
	numPoints = 0;
	numCells = 0;
	xCoordsSet = {};
	yCoordsSet = {};
	xCoords = {};
	yCoords = {};
	values = {};

	std::string VTUCoords;
	std::string VTUValues;

	std::ifstream VTUFile(filename);
	if(!VTUFile.good())
		return false;

	for(std::string line; getline(VTUFile,line ); )
	{
		if (line.compare("        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">") == 0)
		{
			std::string pointsLine;
			getline(VTUFile, pointsLine);
			VTUCoords += pointsLine;
		}
		if (line.compare("        <DataArray type=\"Float32\" Name=\"c\" NumberOfComponents=\"1\" format=\"ascii\">") == 0)
		{
			std::string valuesLine;
			getline(VTUFile, valuesLine);
			VTUValues += valuesLine;
		}
	}
	VTUFile.close();

	//---------------------------------
	//-- Get Coordinates
	//---------------------------------

	int pos = 0;
	std::istringstream pointStream(VTUCoords);
	while (pointStream)
	{
		double coord;
		pointStream >> coord;

		if(pos==0)
		{
			xCoordsSet.insert(coord);
			xCoords.push_back(coord);
			if(coord>this->maxXValue){this->maxXValue = coord;}
			if(coord<this->minXValue){this->minXValue = coord;}
		}
		else if(pos==1)
		{
			yCoordsSet.insert(coord);
			yCoords.push_back(coord);
			if(coord>this->maxYValue){this->maxYValue = coord;}
			if(coord<this->minYValue){this->minYValue = coord;}		
		}
		pos += 1;
		if(pos>2)
		{
			pos=0;
		}
	}

	this->numXCoords = xCoordsSet.size();
	this->numYCoords = yCoordsSet.size();
	this->numXCoordsTotal = xCoords.size();
	this->numYCoordsTotal = yCoords.size();

	this->numPoints = this->numXCoords*this->numYCoords;
	this->numCells = (this->numXCoords-1)*(this->numYCoords-1);

	this->xOffset = (this->maxXValue - this->minXValue) / this->numXCoords;
	this->yOffset = (this->maxYValue - this->minYValue) / this->numYCoords;

	//---------------------------------
	//-- Get Values
	//---------------------------------

	std::istringstream valueStream(VTUValues);
	while (valueStream)
	{
		double singleValue;
		valueStream >> singleValue;
		this->values.push_back(singleValue);
	};
	this->numValues = values.size();

	return true;
}

void VTKReader::get_VTKData(double (&returnDataArray)[100][100]){
	double pointArray [3][this->numPoints];
	
	for(int i=0; i<this->numPoints; i++){
		pointArray[0][i] = this->xCoords[i]*(this->numXCoords-1);
		pointArray[1][i] = this->yCoords[i]*(this->numYCoords-1);		
		pointArray[2][i] = this->values[i];
		/*
		std::cout.precision(6);
		std::cout << "pointArray " << i << ": " << pointArray[0][i] << " " << pointArray[1][i] << " " << pointArray[2][i] << std::endl;
		*/
	}
	
	for(int i=0; i<this->numPoints; i++){
		returnDataArray[(int) pointArray[0][i]][(int) pointArray[1][i]] = pointArray[2][i];
	}

	/*
	std::cout << "maxXValue: " << this->maxXValue << "\n";
	std::cout << "minXValue: " << this->minXValue << "\n";
	std::cout << "maxYValue: " << this->maxYValue << "\n";
	std::cout << "minYValue: " << this->minYValue << "\n";
	std::cout << "numXCoords: " << this->numXCoords << "\n";
	std::cout << "numYCoords: " << this->numYCoords << "\n";
	std::cout << "numPoints: " << this->numPoints << "\n";
	std::cout << "numCells: " << this->numCells << "\n";

	std::cout << "xOffset: " << this->xOffset << "\n";
	std::cout << "yOffset: " << this->yOffset << "\n";
	*/
	
}


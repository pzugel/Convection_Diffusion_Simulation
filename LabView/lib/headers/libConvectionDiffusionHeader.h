int wrapperInit(
		int dim, 
		const char* fileName,
		const char* discType, 
		const char* upwindType, 
		const char* subsets,
		bool setDiffusionValue,
		double diffusionValue,
		double diffusionMatrix[3][3],	
		double reactionRateValue,
		double sourceValue,
		double massScaleValue,
		double velocityVector[3][1],
		const char* precondName,
		const char* solverName,
		int numRefs,
		int maxNumIter,
		double absTol);

void saveVTKFile(const char* path);
char* displayLog();
void clearLog();
double getVTK_minXValue();
double getVTK_minYValue();
double getVTK_maxXValue();
double getVTK_maxYValue();
int getVTK_numXCoords();
int getVTK_numYCoords();
int getVTK_numPoints();
int getVTK_numCells();
double getVTK_xOffset();
double getVTK_yOffset();
void getVTK_Data(double (&vtkdata)[100][100]);

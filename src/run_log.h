#include <cstring>
#include <string>

class RunLog{
	private :
		std::string runLogString = "*********** Run Logger ***********\n\n";
	public :
		void writeLog(std::string add);
		void clearLog();
		char* displayLog();
};

void RunLog::writeLog(std::string add){
	this->runLogString += add;
}

void RunLog::clearLog(){
	this->runLogString = "*********** Run Logger ***********\n\n";
}

char* RunLog::displayLog(){
	char * cstr = new char [this->runLogString.length()+1];
  	std::strcpy (cstr, this->runLogString.c_str());
	return cstr;
}


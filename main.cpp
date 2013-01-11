#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <boost/lexical_cast.hpp>
#include <queue>
#include <vector>
#include <map>
#include "Tissue.h"
using namespace std;

#define MAX_X_POS  50
#define MAX_Y_POS  50
#define MAX_Z_POS  50
#define MIN_X_POS -50
#define MIN_Y_POS -50
#define MIN_Z_POS -50
#define NUM_SIDES   6
#define MAX_ANTIBODY_STRENGTH  100
#define	MIN_ANTIBODY_STRENGTH    0 
#define MAX_INFECTION_STRENGTH 100
#define	MIN_INFECTION_STRENGTH   1 
const string logFileName = "infectionLog.txt";

class TissueReactor; 

struct simStats_t {
    TissueReactor *reactor;     
    int totalInfectedCells;  
    int totalInfectionAttempts; 
    int totalStrengthDelta;	/* Sum of difference between disease strength and antibody strength. */ 
    int totalCytotoxicCells; 
    int totalHelperCells; 
    int infectionSpread;	/* Volume of the smallest rectangular box containing all infected cells. */
    int longestInfectionPath; 
};

struct qElement_t{
    CellMembrane::Side side;  /* The side of the cell from which the infection is trying to penetrate */
    Cell::Coordinates coord;  /* The coordinates of the cell. */
    int BFSdepth;	      /* To keep track of the degree of separation from the infection root. */ 
};

map<string, Tissue::Ptr> tissueMap; /* Data structure to keep track of all the tissues created */ 
map<string, simStats_t *> statsMap; /* Data structure to keep track of all the stats of each tissue */ 

class TissueReactor : public Tissue::Notifiee {
    public:
	void onCellNew(Cell::Ptr cell){
	    AntibodyStrength abStrength; 
	    if(cell->cellType() == Cell::cytotoxicCell()){
		abStrength.valueIs(MAX_ANTIBODY_STRENGTH); 
		numCytotoxicCells_++;
	    } 
	    if(cell->cellType() == Cell::helperCell()){
		abStrength.valueIs(MIN_ANTIBODY_STRENGTH); 
		numHelperCells_++; 
	    }
		
	    CellMembrane::Side sides[] = {CellMembrane::north(), CellMembrane::south(), CellMembrane::east(), CellMembrane::west(), CellMembrane::up(), CellMembrane::down()};
	    for(int i = 0; i < NUM_SIDES; i++){
		CellMembrane::Ptr membrane = cell->membraneNew(cell->name(), sides[i]);
		membrane->antibodyStrengthIs(abStrength);
	    }
	}
	
	void onCellDel(Cell::Ptr cell){
	    CellMembrane::Side sides[] = {CellMembrane::north(), CellMembrane::south(), CellMembrane::east(), CellMembrane::west(), CellMembrane::up(), CellMembrane::down()};
	    for(int i = 0; i < NUM_SIDES; i++){
	        CellMembrane::Ptr membrane = cell->membraneDel(sides[i]);
	    }
	    if(cell->cellType() == Cell::cytotoxicCell()){
		numCytotoxicCells_--;
	    } 
	    if(cell->cellType() == Cell::helperCell()){
		numHelperCells_--; 
	    }
	}
	
	static TissueReactor * TissueReactorIs(Tissue *t){
	    TissueReactor *m = new TissueReactor(t);
	    return m; 
	}
	
	int numCytotoxicCells(){
	    return numCytotoxicCells_; 
	}

	int numHelperCells(){
	    return numHelperCells_; 
	}

    protected:
	TissueReactor(Tissue *t) : Tissue::Notifiee(){
	    numCytotoxicCells_	= 0; 
	    numHelperCells_	= 0; 
	    notifierIs(t);
	}

    private: 
	int numCytotoxicCells_ ; 
	int numHelperCells_; 

}; 

void explodeLine(vector<string> & words, const string line); /* Explode the line: split it by spaces */
void parseLine(const string textline);
void readInputFile(ifstream &infile);
int  parseXYZ(vector<string> &fields, int xpos, int ypos, int zpos, int &x, int &y, int &z);
bool areValidCoordinates(int x, int y, int z);
void decodeCellCommand(vector<string> &fields, const string cmd);
void decodeTissueCommand(vector<string> &fields, const string cmd);
void newCellCommand(const string cellType, const string tissueName, int x, int y, int z );
void cloneCellsCommand(const string tissueName, const string direction);
void shiftedCoordinates(Cell::Coordinates &shiftedCoord, Cell::Coordinates coord, string direction); 
void cloneOneCellCommand(Cell::Coordinates refCoord, Tissue::Ptr tissue, const string direction);
void infectionSimulation(Tissue::Ptr tissue, simStats_t *stats, Cell::Coordinates infectionOriginCoord, string direction, int infectionStrength);
void enqueueNeighbors(queue<qElement_t> &cellQueue, Cell::Coordinates coord, int neighborDepth); 
void computeTissueStateStats(Tissue::Ptr tissue, simStats_t *stats); 
void statsInit(simStats_t *stats);
void printStats(simStats_t *stats);

/*
  The main takes in one input, the file name with the rules.
  The rules are then executed and the appropriate statistics are printed
  to the console.
*/
int main(int argc, const char* argv[]) {
    ifstream infile(argv[1]);
    if(infile.fail()){
      cout << "Error reading file" << endl;
      return 1;
    }
    ofstream logFile(logFileName.c_str());
    cerr << "Simulation: Running configuration file " << argv[1] << endl;

    readInputFile(infile);    
    return 0;
}

void readInputFile(ifstream &infile){
    string textline;
    while(!infile.eof()){
	getline(infile, textline);
	if(textline[0] == '#' || textline.empty()) 
	    continue;
	cerr << "COMMAND = " << textline << endl;
	parseLine(textline);
	textline.clear();
    }
}

void parseLine(const string textline){
    vector<string> words; 
    explodeLine(words, textline);
    
    if ( words[0] == "Tissue" ) {
	decodeTissueCommand(words, textline);
    } 
    else if (words[0] == "Cell") {
	decodeCellCommand(words, textline);
    } 
    else {
	 /* Should not come here, since the 0th word 
	    must be either Tissue or Cell */
	return; 
    }
    words.clear();
}

/*Return 0 on success and -1 on failure*/
int parseXYZ(vector<string> &fields, int xpos, int ypos, int zpos, int &x, int &y, int &z){
    try {
        x = boost::lexical_cast<int>(fields[xpos]);  
        y = boost::lexical_cast<int>(fields[ypos]);  
        z = boost::lexical_cast<int>(fields[zpos]);  
    } 
    catch( boost::bad_lexical_cast const & ) {
	return -1; 
    }
    return 0; 
}

bool areValidCoordinates(int x, int y, int z){
    if(x < MIN_X_POS || x > MAX_X_POS)
	return false;
    if(y < MIN_Y_POS || y > MAX_Y_POS)
	return false;
    if(z < MIN_Z_POS || z > MAX_Z_POS)
	return false;
    
    return true; 
}

void newCellCommand(const string cellType, const string tissueName, int x, int y, int z ){
    if(!areValidCoordinates(x,y,z)){
        cerr << "\tEXCEPTION: The coordinates to create the new cell ("<< x << "," << y << "," << z << ") are not valid. Skipping command." << endl; 
        return; 
    }
    Cell::Coordinates coord = Cell::Coordinates();
    coord.x = x;
    coord.y = y;
    coord.z = z;
    
    Tissue::Ptr tissue = tissueMap[tissueName];
    
    Cell::Ptr cell = tissue->cell(coord);
    if(cell.ptr() != NULL){
    	cerr << "\tEXCEPTION: A Cell is already present at the coordinates ("<< x << "," << y << "," << z << "). Skipping command." << endl; 
    	return;
    }
    
    if(cellType == "cytotoxicCell") {
	tissue->cellIs(Cell::CellNew(coord, tissue.ptr(), Cell::cytotoxicCell() ));
	Cell::Ptr cell = tissue->cell(coord);
    } else if (cellType == "helperCell"){ 
	tissue->cellIs(Cell::CellNew(coord, tissue.ptr(), Cell::helperCell() ));
	Cell::Ptr cell = tissue->cell(coord);
    } else {
	cerr << "\tEXCEPTION: Unrecognized cell type " << cellType << ".Skipping command." << endl;
    }
}

void shiftedCoordinates(Cell::Coordinates &shiftedCoord, Cell::Coordinates coord, string direction){
    shiftedCoord.x = coord.x;
    shiftedCoord.y = coord.y;
    shiftedCoord.z = coord.z;
    if(direction == "north"){
	shiftedCoord.y++;
    } else if (direction == "south"){
	shiftedCoord.y--;
    } else if (direction == "east"){
	shiftedCoord.x++;
    } else if (direction == "west"){
	shiftedCoord.x--;
    } else if (direction == "up") {
	shiftedCoord.z++;
    } else {
	shiftedCoord.z--;
    }
} 

void cloneCellsCommand(const string tissueName, const string direction){
    Tissue::Ptr tissue = tissueMap[tissueName];
    vector<Cell::Coordinates> parents; /* "parents" is a list of the currently existing cells in the tissue. */
    for(Tissue::CellIteratorConst cell = tissue->cellIterConst(); cell; ++cell){
	parents.push_back(cell->location()); 
    }
    for(unsigned int i = 0; i < parents.size(); i++){
       	cloneOneCellCommand(parents[i], tissue, direction);
    }
}

void cloneOneCellCommand(Cell::Coordinates refCoord, Tissue::Ptr tissue, const string direction){
    if(!areValidCoordinates(refCoord.x, refCoord.y, refCoord.z)){
        cerr << "\tEXCEPTION: The reference coordinates ("<< refCoord.x << "," << refCoord.y << "," << refCoord.z << ") are not valid. Skipping command." << endl;
        return; 
    }
    Cell::Ptr ref = tissue->cell(refCoord);
    if(ref.ptr() == NULL) {
        /* No cell exists here to clone from. Skip. */
	cerr << "\tEXCEPTION: No cell here ("<< refCoord.x << "," << refCoord.y << "," << refCoord.z << ") to clone from. Skipping." << endl;
        return;
    } else {
        Cell::Coordinates cloneCoord = Cell::Coordinates();
        shiftedCoordinates(cloneCoord, refCoord, direction); 
	if(!areValidCoordinates(cloneCoord.x, cloneCoord.y, cloneCoord.z)){
    	    cerr << "\tEXCEPTION: The clone location coordinates ("<< cloneCoord.x << "," << cloneCoord.y << "," << cloneCoord.z << ") not valid. Skipping." << endl; 
    	    return; 
    	}
        /* coord updated to have the coordinates of the 
           position where the cell is to be cloned. */
        Cell::Ptr temp = tissue->cell(cloneCoord);
        if(temp.ptr() != NULL){
            /* A cell already exists in this location, can't clone into this location */
	    cerr << "\tEXCEPTION: Cell already exists in the destination coordinates ("<< cloneCoord.x << "," << cloneCoord.y << "," << cloneCoord.z << "). Skipping command." << endl;
	    return; 	
        } else {
    	/* clone the cell : type, health, membrane strength */
	   tissue->cellIs(Cell::CellNew(cloneCoord, tissue.ptr(), ref->cellType() ));
    	   Cell::Ptr clone = tissue->cell(cloneCoord);
    	   clone->healthIs(ref->health());
	   CellMembrane::Side sides[] = {CellMembrane::north(), CellMembrane::south(), CellMembrane::east(), CellMembrane::west(), CellMembrane::up(), CellMembrane::down()};
    	   for(int i = 0; i < NUM_SIDES; i++){
    	       CellMembrane::Ptr refMembrane = ref->membrane(sides[i]);
    	       CellMembrane::Ptr cloneMembrane = clone->membrane(sides[i]);
    	       cloneMembrane->antibodyStrengthIs(refMembrane->antibodyStrength());
    	   }
        }
    }
}

void enqueueNeighbors(queue<qElement_t> &cellQueue, Cell::Coordinates coord, int neighborDepth){ 
    /* The order in which the cells are enqueued is 
       north, east, south, west, up, down */
    struct qElement_t next; 
    next.BFSdepth = neighborDepth; 

    shiftedCoordinates(next.coord, coord, "north"); /* Coordinates of the cell to be infected */
    next.side = CellMembrane::south(); /* Direction of infection attempt w.r.t the cell being infected. */
    cellQueue.push(next);

    shiftedCoordinates(next.coord, coord, "east");
    next.side = CellMembrane::west();
    cellQueue.push(next);

    shiftedCoordinates(next.coord, coord, "south");
    next.side = CellMembrane::north();
    cellQueue.push(next);

    shiftedCoordinates(next.coord, coord, "west");
    next.side = CellMembrane::east();
    cellQueue.push(next);

    shiftedCoordinates(next.coord, coord, "up");
    next.side = CellMembrane::down();
    cellQueue.push(next);

    shiftedCoordinates(next.coord, coord, "down");
    next.side = CellMembrane::up();
    cellQueue.push(next);
}

void computeTissueStateStats(Tissue::Ptr tissue, simStats_t *stats){
    stats->totalCytotoxicCells = stats->reactor->numCytotoxicCells(); 
    stats->totalHelperCells = stats->reactor->numHelperCells(); 
    stats->totalInfectedCells = 0; 

    int maxX = MIN_X_POS, maxY = MIN_Y_POS, maxZ = MIN_Z_POS, minX = MAX_X_POS, minY = MAX_Y_POS, minZ = MAX_Z_POS; 
    bool pristine = true; 
    for(Tissue::CellIterator cell = tissue->cellIter(); cell; ++cell){
        if(cell->health() == Cell::infected()){
	    stats->totalInfectedCells++; 
            pristine = false; 
            Cell::Coordinates coord = cell->location(); 
            if(coord.x > maxX) maxX = coord.x;
            if(coord.y > maxY) maxY = coord.y;
            if(coord.z > maxZ) maxZ = coord.z;
            if(coord.x < minX) minX = coord.x;
            if(coord.y < minY) minY = coord.y;
            if(coord.z < minZ) minZ = coord.z;
        }			
    }
    if(!pristine){
        stats->infectionSpread = ( maxX - minX + 1 )*( maxY - minY + 1 )*( maxZ - minZ + 1 ); 
    }
}

void infectionSimulation(Tissue::Ptr tissue, simStats_t *stats, Cell::Coordinates infectionOriginCoord, string direction, int infectionStrength){
    queue<qElement_t> cellQueue;  
    
    qElement_t root;  
    root.side = CellMembrane::SideInstance(direction); 
    root.coord = infectionOriginCoord; 
    root.BFSdepth = 0; 
    cellQueue.push(root);

    while(!cellQueue.empty()){
	qElement_t curr = cellQueue.front();
	cellQueue.pop();
	Cell::Ptr cell = tissue->cell(curr.coord);
	CellMembrane::Side side = curr.side;
	int depth = curr.BFSdepth; 
	stats->longestInfectionPath = depth; 
	if(cell.ptr() == NULL){
	    continue;
	} else {
	    if(cell->health() == Cell::infected()){
	        continue; 
	    } else { /* the cell is healthy, make an infection attempt */
		stats->totalInfectionAttempts++; 
	        CellMembrane::Ptr membrane = cell->membrane(side);  
		int antibodyStrength = (membrane->antibodyStrength()).value(); 
		stats->totalStrengthDelta += (infectionStrength - antibodyStrength); 
	        if( antibodyStrength < infectionStrength ){
		   /* Infect this cell and prepare to propagate */
		    cell->healthIs(Cell::infected());
		    enqueueNeighbors(cellQueue, cell->location(),depth + 1); 
	        }
	    }
	}
    }
    computeTissueStateStats(tissue, stats); 
}

void statsInit(simStats_t *stats){
    stats->totalInfectedCells	   = 0;   
    stats->totalInfectionAttempts  = 0; 
    stats->totalStrengthDelta	   = 0; 
    stats->totalCytotoxicCells	   = 0; 
    stats->totalHelperCells	   = 0; 
    stats->infectionSpread	   = 0; 
    stats->longestInfectionPath	   = 0; 
}

void printStats(simStats_t *stats){
    int a,b,c,d,e,f,g;
    a = stats->totalInfectedCells;  
    b = stats->totalInfectionAttempts;
    c = stats->totalStrengthDelta;   
    d = stats->totalCytotoxicCells; 
    e = stats->totalHelperCells; 
    f = stats->infectionSpread;      
    g = stats->longestInfectionPath; 
    cout << a << " " << b << " " << c << " " << d << " " << e << " " << f << " " << g << endl; 
} 

void decodeTissueCommand(vector<string> &fields, const string cmd) {
    string tissueName; 
    string operation; 
    string direction; 
    int x, y, z, infectionStrength, status;

    if(fields[1] == "tissueNew") {
        tissueName.assign(fields[2]);
	
	Tissue::Ptr tissue = Tissue::TissueNew(tissueName); 
	tissueMap[tissueName] = tissue; 
		
	/* Declare an instance of the tissue reactor and 
	   link it to the newly created tissue.*/
	TissueReactor *m = TissueReactor::TissueReactorIs(tissue.ptr());

	simStats_t *simStats = new simStats_t;
	statsInit(simStats);
	simStats->reactor = m; 
	statsMap[tissueName] = simStats; 

	tissueName.clear();
    } else {
        tissueName.assign(fields[1]);	
        operation.assign(fields[2]);
	bool validTissue = false; 
	for(map<string, Tissue::Ptr>::iterator it = tissueMap.begin(); it != tissueMap.end(); it++) {
	    if(it->first == tissueName)
		validTissue = true;
    	} 

	if( operation == "cytotoxicCellNew") {
	    if(!validTissue){
	        cerr << "\tEXCEPTION: Non-existant Tissue name = " << tissueName << ". Skipping Command." << endl;  
	        return; 
	    }
	    status = parseXYZ(fields, 3,4,5,x,y,z);
	    if(status == -1)
		return; 
	    newCellCommand("cytotoxicCell", tissueName, x,y,z);
    	} 
	else if (operation == "helperCellNew") {
	    if(!validTissue){
	        cerr << "\tEXCEPTION: Non-existant Tissue name = " << tissueName << ". Skipping Command." << endl;  
	        return; 
	    }
    	    status = parseXYZ(fields, 3,4,5,x,y,z);
	    if(status == -1)
		return; 
	    newCellCommand("helperCell", tissueName, x,y,z);
	}
	else if ( operation == "infectionStartLocationIs") {
	    simStats_t *stats = new simStats_t; 
	    statsInit(stats);
	    if(!validTissue){
	        cerr << "\tEXCEPTION: Non-existant Tissue name = " << tissueName << ". Skipping Command." << endl;  
		printStats(stats); 
	        return; 
	    }
	    Tissue::Ptr tissue = tissueMap[tissueName];
	   
	    /* Start working with a fresh statistics data structure for this run. */
	    stats->reactor = statsMap[tissueName]->reactor; 
	    delete statsMap[tissueName]; 

	    computeTissueStateStats(tissue, stats); 

	    status = parseXYZ(fields, 3,4,5,x,y,z);	
	    if(status == -1 || !areValidCoordinates(x, y, z)){
    	        cerr << "\tEXCEPTION: The infection start location coordinates ("<< x << "," << y << "," << z << ") are not valid. Skipping command." << endl; 
		printStats(stats);
    	        return; 
    	    }
	    Cell::Coordinates infectionOriginCoord = Cell::Coordinates();
    	    infectionOriginCoord.x = x; 
    	    infectionOriginCoord.y = y;
    	    infectionOriginCoord.z = z;
    	    direction.assign(fields[6]); 
    	    try {
    	        infectionStrength = boost::lexical_cast<int>(fields[7]);		 
    	    }
    	    catch (boost::bad_lexical_cast const &) {
    	        cerr << "\tEXCEPTION! Unable to cast the infection strength string to integer. Skipping command." << endl;
		printStats(stats);
		return; 
    	    }
	    if(infectionStrength > MAX_INFECTION_STRENGTH || infectionStrength < MIN_INFECTION_STRENGTH){
	        cerr << "\tEXCEPTION: Infection strength = " << infectionStrength << "is not between allowed bounds. Skipping command." << endl;
		printStats(stats);
		return; 
	    } 
	    infectionSimulation(tissue, stats, infectionOriginCoord, direction, infectionStrength);
	    statsMap[tissueName] = stats; 
	    printStats(stats);
    	}
    	else if ( operation == "infectedCellsDel" ) {
	    if(!validTissue){
	        cerr << "\tEXCEPTION: Non-existant Tissue name = " << tissueName << ". Skipping Command." << endl;  
	        return; 
	    }
    	    /* call the infected cell delete function*/
	    Tissue::Ptr tissue = tissueMap[tissueName];
	    vector<string> infectedCellList; 
	    for(Tissue::CellIterator it = tissue->cellIter(); it ; ++it) {
		if (it->health() == Cell::infected()){
		    infectedCellList.push_back(it->name());
		}
	    }
	    for(unsigned int i = 0; i < infectedCellList.size(); i++){
		tissue->cellDel(infectedCellList[i]);
	    }
    	}
    	else if ( operation ==  "cloneCellsNew" ) {
	    if(!validTissue){
	        cerr << "\tEXCEPTION: Non-existant Tissue name = " << tissueName << ". Skipping Command." << endl;  
	        return; 
	    }
    	    direction.assign(fields[3]);	
	    cloneCellsCommand(tissueName, direction);
    	}
	else {
    	    /* Should not come here, since the 6th field in  
    	       the command must be either cloneNew or membrane */ 
    	    return; 
	}	
    }
}

void decodeCellCommand(vector<string> & fields,const string cmd){
    string component;	/* Tissue or Cell */
    string tissueName;
    int x, y, z, abStrength, status; 
    string direction;  /* north, south, east, west, up, down */
    string operation;  

    component.assign(fields[0]);
    tissueName.assign(fields[1]);
    operation.assign(fields[5]);
    direction.assign(fields[6]);
   
    bool validTissue = false; 
    for(map<string, Tissue::Ptr>::iterator it = tissueMap.begin(); it != tissueMap.end(); it++) {
        if(it->first == tissueName)
    	validTissue = true;
    } 
    if(!validTissue){
	cerr << "\tEXCEPTION: Non-existant Tissue name = " << tissueName << ". Skipping Command." << endl;  
	return; 
    }

    Tissue::Ptr tissue = tissueMap[tissueName]; 

    status = parseXYZ(fields, 2,3,4,x,y,z);
    if(status == -1) 
	return;
    if(!areValidCoordinates(x, y, z)){
        cerr << "\tEXCEPTION: The supplied coordinates ("<< x << "," << y << "," << z << ") are not valid. Skipping command." << endl;
        return; 
    }
    Cell::Coordinates coord = Cell::Coordinates();
    coord.x = x; 
    coord.y = y;
    coord.z = z;

    if(operation == "cloneNew"){
	
	cloneOneCellCommand(coord, tissue, direction);
    }
    else if (operation == "membrane") {
	try {
	    abStrength = boost::lexical_cast<int>(fields[8]);		 
	}
	catch (boost::bad_lexical_cast const &) {
	    cerr << "\tEXCEPTION! Unable to cast antibody strength string to integer. Skipping command." << endl;
	    return; 
	}
	
	if(abStrength > MAX_ANTIBODY_STRENGTH || abStrength < MIN_ANTIBODY_STRENGTH){
	    cerr << "\tEXCEPTION: Antibody strength = " << abStrength << " is not between allowed bounds. Skipping command." << endl;  
	    return; 
	}
	
	Cell::Ptr cell = tissue->cell(coord);
	CellMembrane::Side side = CellMembrane::SideInstance(direction); 	
	CellMembrane::Ptr membrane = cell->membrane(side);
	AntibodyStrength strength; 
	strength.valueIs(abStrength);
	membrane->antibodyStrengthIs(strength);
    }
    else{
         /* Should not come here, since the 6th field in  
            the command must be either cloneNew or membrane */ 
	return; 
    }    
}

void explodeLine(vector<string> & words, const string line){
    string word;
    int wordStartPos = 0;
    int wordStopPos  = 0; 
    bool explosionDone = false;

    while(!explosionDone){
        size_t spacePos = line.find_first_of(' ',wordStartPos);
        if(spacePos != string::npos){
            wordStopPos = spacePos - 1;
        } else {
            /* no more spaces in the line, all the words are done. */
            wordStopPos = line.size();
            explosionDone = true;
        }
        int numChars = wordStopPos - wordStartPos + 1;
        word.assign(line, wordStartPos, numChars);
        words.push_back(word);
        word.clear();
        wordStartPos = spacePos + 1; 
    }
}

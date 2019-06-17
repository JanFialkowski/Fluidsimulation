#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <climits>
#include <random>

using namespace std;

typedef vector< vector<double> >   vec2Ddouble;
typedef vector< vector<unsigned> > vec2Dunsigned;

typedef mt19937 MyRNG;            //Set the Mersenne-Twister pseudorandom number generator
                                  //as the type of my random number generator

#define SQ(a) ((a)*(a))

/**** Class Simulation ***/

class Simulation {

  public:
    Simulation();           //Constructor
    ~Simulation(){};        //Deconstructor
    /*** Public Functions ***/

    //Initialization of the simulation
    void initialize(string & path);
    //Computation
    void update();
    //Output info 
    void output(unsigned long int step, string path);

    unsigned long int get_nStep(){ return _nStep; } 
    unsigned long int get_outputPeriod(){ return _outputPeriod; }
    
  private:
    unsigned        _nPar;            //Number of particles
    vec2Ddouble        _r;            //Particle positions
    vec2Ddouble        _f;            //Forces due to particle interactions
    vec2Ddouble       _vr;            //Translational noise
    vector<double>     _u;            //The interaction potential energy
                             
    double           _eps;            //The interaction strength (in unit of thermal energy)
    double         _sigma;            //The size of a particle   (length unit)
    double            _rc;            //The cut-off radius
                             
    double            _gt;            //The translational drag constant gamma_t (unit of drag constant)
    double            _Dt;            //The translational diffusion constant
    double            _gr;            //The rotational drag constant
    double            _Dr;            //The rotational diffusion constant
 
    double            _kBT;           //The thermal energy (energy unit)
                             
    double            _v0;            //Propulsion speed
    vec2Ddouble        _e;            //Particle orientation
    vector<double>   _phi;            //Orientational angle
    vector<double>    _wr;            //Rotational noise
                             
    bool     _ifUseNbrListAlgorithm;  //whether to employ neighbor list algorithm
    vec2Dunsigned            _nList;  //The neighbor list for implementing
                                      //Verlet list algorithm
    vec2Ddouble               _rOld;  //Particle positions when _nList was updated.
    double                      _rL;  //Distance criterion 

    unsigned long int        _nStep;  //Number of simulation steps
    unsigned long int _outputPeriod;  //Output info every <outputPeriod> steps
    double                      _dt;  //the time interval in one simulation step
                                      //(in unit of tau = _sigma^2/_Dt)

    vector<double>         _boxSize;  //The size of the simulation box

    MyRNG                      _rng;  //My random number generator
    MyRNG::result_type    _seed_val;  //The random seed
    normal_distribution<double> _normal_dist;  // N(mean, stddeviation)

    //private functions
    void initVariables(string & path);
    void initializeConfig();

    void computeInteractions();
    //Different types of interactions
    void            Lennard_Jones(unsigned i, unsigned j);
    vector<double>  calcForceLennard_Jones(vector<double> imgD,
                                           double r2, 
                                           double sigma_dev_r_6,
                                           double sigma_dev_r_12);
    void            Weeks_Chandler_Andersen(unsigned i, unsigned j);
    vector<double>  calcForceWCA(vector<double> imgD,
                                 double r2, 
                                 double sigma_dev_r_6,
                                 double sigma_dev_r_12);
    void computeThermalNoise();

    //void computePtlPos();
    void computePtlPosAndOrnt();   //compute particles' positions and orientations

    void outputConfig(unsigned long int step, string path);
    double calcEnergyPerPtl();

    //Neighbor list
    void updateNbr();
    bool checkIfUpdateNbr();

    //Handling periodic boundary conditions
    vector<double> calcImgDisplacement(unsigned i, unsigned j);  //rij = rj - ri
};

//Useful functions
double dotProduct2D(vector<double> p, vector<double> q);

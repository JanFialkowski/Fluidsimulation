#include "LJ.h"

int main(int argc, char* argv[])
{
  string path;

  if(argc != 2){
    cout << "Usage: " << argv[0] << " <directory>" << endl;
    exit(-1);
  }

  path = argv[1];
  initialize(path);

  for(unsigned long int i = 0 ; i < _nStep ; ++i){ // i = 0 - (_nStep - 1)
    if(i%_outputPeriod == 0){
      cout << "Current time step = " << i << endl;
      output(i, path);
    }

    update();
  }

  if(_nStep % _outputPeriod == 0){  //Output info for the last step
    cout << "Current time step = " << _nStep << endl;
    output(_nStep, path);
  }

  return 0;
}

void initialize(string & path)
{ 

  cout << "==== Initializing the Simulation ====" << endl;

  //Parse input file and assign variables
  initVariables(path);

  //Generate an initial configuration
  initializeConfig();
}

void initVariables(string & path)
{
  ifstream iFile;
  string key, key2;
 
  if(path.back() != '/')
    path = path + '/';

  string fileName = path + "input.dat";

  iFile.open(fileName.c_str());
  if(!iFile.is_open()){
    cout << "Failed to open input file <<" << fileName << ">>" << endl;
  }
  else{
    while(iFile.good()){
      iFile >> key;
      if(key == "nPar"){
        iFile >> _nPar;
        cout << "nPar = " << _nPar << endl;
      }
      else if(key == "kBT"){
        iFile >> _kBT;
        cout << "kBT = " << _kBT << endl;
      }
      else if(key == "ifUseNbrList"){
        iFile >> key2;
        if(key2 == "yes"){
          _ifUseNbrListAlgorithm = true; 
          cout << "Neighbor List Algorithm is employed." << endl;
        }
        else if(key2 == "no"){
          _ifUseNbrListAlgorithm = false; 
          cout << "Simulation is performed without introducing Nbr List." << endl;
        }
        else{
          cout << "Error found in initVariables(string & path)!" << endl;
          cout << "When reading ifUseNbrList, options other than <yes> or <no> are not allowed,";
          cout << "and it is case sensitive." << endl;
          cout << "The wrong input option is " << key2 << endl;
          cout << "Program exiting ... " << endl;
          exit(-1);
        }
      }
      else if(key == "nStep"){
        iFile >> _nStep;
        cout << "nStep = " << _nStep << endl;
      }
      else if(key == "outputPeriod"){
        iFile >> _outputPeriod;
        cout << "outputPeriod = " << _outputPeriod << endl;
      }
      else if(key == "dt"){
        iFile >> _dt;
        cout << "dt = " << _dt << endl;
      }
      else if(key == "boxSize"){
        iFile >> _boxSize[0];
        iFile >> _boxSize[1];
        cout << "boxSize = {" << _boxSize[0] << ", " << _boxSize[1] << "}" << endl;
      }
      else if(key == ""){
        continue;
      }
      else if(key == "END_OF_FILE"){
        break;
      }
      else{
        cout << "key value " << key << " unknown." << endl;
        cout << "Error found in initVariables(string & path)!" << endl;
        cout << "Program exiting ..." << endl;
        exit(-1);
      }
    }
  }
    
  iFile.close();

  //Allocate sizes and assign values to variables
  _r.assign(_nPar, vector<double>(2, 0.));
  _f.assign(_nPar, vector<double>(2, 0.));
  _vr.assign(_nPar, vector<double>(2, 0.));

  _u.assign(_nPar, 0.);

  _eps = 1.; //unit
  _sigma = 1.; //unit
  _rc = 2.5*_sigma;

  _gt = 1.; //unit
  _Dt = _kBT/_gt;

  _nList.assign( _nPar - 1, vector<unsigned>(1, 0));
  for(unsigned i = 0 ; i < _nPar - 1 ; ++i) //upper triangle
    _nList[i].assign(_nPar - i, 0);
  
  _rOld.assign(_nPar, vector<double>(2, 0.));

  _rL = 3.5*_sigma;

  //Initialize random number generator
  _seed_val = time(0);
  _rng.seed(_seed_val);

}

void initializeConfig()
{
  unsigned nn = (unsigned) ceil(sqrt(_nPar));
  vector<double> dL(2, 0.);

  for(unsigned i = 0 ; i < 2 ; ++i)
    dL[i] = _boxSize[i]/nn;

  for(unsigned i = 0 ; i < _nPar ; ++i){
    _r[i][0] = (0.5 + (i/nn))*dL[0];
    _r[i][1] = (0.5 + (i%nn))*dL[1];
    _rOld[i][0] = _r[i][0];
    _rOld[i][1] = _r[i][1];
  }
}


//Computation
void update()
{
  computeInteractions();
  computeThermalNoise();
  computePtlPos();

  if(_ifUseNbrListAlgorithm){
    if(checkIfUpdateNbr()){
      updateNbr();
    }
  }

}

void computeInteractions()
{
  for(unsigned i = 0 ; i < _nPar ; ++i){
    _u[i] = 0.;
    for(unsigned j = 0 ; j < 2 ; ++j)
      _f[i][j] = 0.;
  }

  for(unsigned i = 0 ; i < _nPar ; ++i){
    if(_ifUseNbrListAlgorithm){
      if( i == _nPar - 1)
        continue;
      for(unsigned j, k = 1 ; k <= _nList[i][0] ; ++k){
        j = _nList[i][k];
        Lennard_Jones(i,j);
      }
    }
    else{
      for(unsigned j = i + 1; j < _nPar ; ++j){
        Lennard_Jones(i,j);
      }
    }
  }
}

void computeThermalNoise()
{
  static double preFactor = sqrt(2*_Dt/_dt); //set to zero to turn off the noise when debugging

  for(unsigned i = 0 ; i < _nPar ; ++i){
    for(unsigned j = 0 ; j < 2 ; ++j){
      _vr[i][j] = preFactor * normal_dist(_rng); 
    }
  }
}

void computePtlPos()
{
  for(unsigned i = 0 ; i < _nPar ; ++i){
    for(unsigned j = 0 ; j < 2 ; ++j)
      _r[i][j] += _dt*(_f[i][j]/_gt+_vr[i][j]);
  }
}

//Different types of interactions
void Lennard_Jones(unsigned i, unsigned j)
{
  vector<double> imgD = calcImgDisplacement(i,j);

  double r2  = SQ(imgD[0]) + SQ(imgD[1]);
  
  if(r2 > _rc*_rc){
    return;
  }
  else{
    double sigma_3 = _sigma*_sigma*_sigma;
    double sigma_dev_r_6  = sigma_3*sigma_3/(r2*r2*r2);
    double sigma_dev_r_12 = sigma_dev_r_6*sigma_dev_r_6;
    
    double Uij =  4*_eps*(sigma_dev_r_12 - sigma_dev_r_6);
    vector<double> fij = calcForceLennard_Jones(imgD, r2, sigma_dev_r_6, sigma_dev_r_12);

    if(SQ(fij[0])+SQ(fij[1])>1e8){
      cout << "Force between " << i << " and " << j << " is too large!" << endl;
      cout << "fij = {" << fij[0] << ", " << fij[1] << "}" << endl;
    }
    
    _u[i] += 0.5*Uij;
    _u[j] += 0.5*Uij;
    
    for(unsigned k = 0 ; k < 2 ; ++k){
      _f[i][k] -= fij[k]; //Newton's third law: fji = -fij
      _f[j][k] += fij[k];
    }
  }

}

vector<double>  calcForceLennard_Jones(vector<double> imgD, 
                                       double r2, 
                                       double sigma_dev_r_6,
                                       double sigma_dev_r_12)
{
  vector<double> force(2,0.);

  //fij = -grad_rij(uij);
  for(unsigned i = 0 ; i < 2 ; ++i)
    force[i] = 24.0*_eps*( 2.0*sigma_dev_r_12 - sigma_dev_r_6 )*imgD[i]/r2;

  return force;
}

void Weeks_Chandler_Andersen(unsigned i, unsigned j)
{
  //TODO for purely repulsive systems 
}

//Output info
void output(unsigned long int step, string path)
{
  double P, E;
  ofstream oFile;
 
  P = calcPressure();
  E = calcEnergyPerPtl();

  string fileName = path + "properties.dat";

  if(step == 0){
    oFile.open(fileName.c_str());
    oFile << "#iStep Energy_per_ptl Pressure" << endl;
  }
  else
    oFile.open(fileName.c_str(), ios::app);
    
  oFile << step << " " << E << " " << P << endl;

  oFile.close();

  outputConfig(step, path);
}

double calcPressure()
{
  double area = _boxSize[0]*_boxSize[1];
  double rho = _nPar / area;
  double p = 0.;
  double r2, sigma_3 = _sigma*_sigma*_sigma;
  double sigma_dev_r_6;
  double sigma_dev_r_12;
  vector<double> imgD, force;

  for(unsigned i = 0 ; i < _nPar ; ++i){
    if(_ifUseNbrListAlgorithm){
      if(i == _nPar -1)
        continue;
      for(unsigned j, k = 1 ; k <= _nList[i][0] ; ++k){
        j = _nList[i][k];
      
        imgD = calcImgDisplacement(i,j);
        r2  = SQ(imgD[0]) + SQ(imgD[1]);
        if(r2 < _rc*_rc){
          sigma_dev_r_6  = sigma_3*sigma_3/(r2*r2*r2);
          sigma_dev_r_12 = sigma_dev_r_6*sigma_dev_r_6;
       
          force = calcForceLennard_Jones(imgD, r2, sigma_dev_r_6, sigma_dev_r_12);
          p += dotProduct2D(imgD, force);
        }
      }
    }
    else{
      for(unsigned j = i + 1; j < _nPar ; ++j){
        imgD = calcImgDisplacement(i,j);
        r2  = SQ(imgD[0]) + SQ(imgD[1]);
        if(r2 < _rc*_rc){
          sigma_dev_r_6  = sigma_3*sigma_3/(r2*r2*r2);
          sigma_dev_r_12 = sigma_dev_r_6*sigma_dev_r_6;
      
          force = calcForceLennard_Jones(imgD, r2, sigma_dev_r_6, sigma_dev_r_12);
          p += dotProduct2D(imgD, force);
        }
      }
    }
  }
  p /= (2*area);
  p += rho*_kBT;

  return p;
}

double calcEnergyPerPtl()
{
  double E = 0.;

  for(unsigned i = 0 ; i < _nPar ; ++i){
    E += _u[i];
  }
  E /= _nPar;

  return E;
}

void outputConfig(unsigned long int step, string path)
{
  ofstream oFile;
 
  string fileName = path + "config.xyz";

  if(step == 0){
    oFile.open(fileName.c_str());
  }
  else
    oFile.open(fileName.c_str(), ios::app);
    
  oFile << _nPar << endl;
  oFile << "Lattice=";
  oFile << "\"";
  oFile << _boxSize[0] << " 0.0 "     << "0.0 ";
  oFile << "0.0 "      << _boxSize[1] << " 0.0 ";
  oFile << "0.0 "      << "0.0 "      << "1.0 ";
  oFile << "\" ";
  oFile << "Properties=species:S:1:pos:R:2 ";
  oFile << "Time=" << step*_dt << endl;

  for(unsigned i = 0 ; i < _nPar ; ++i){
    oFile << "M " << _r[i][0] << " " << _r[i][1] << endl;
  }

  oFile.close();

}

//Neighbor list
void updateNbr()
{
  vector<double> imgD;
  double r2;
  unsigned nNbr;

  for(unsigned i = 0 ; i < _nPar - 1 ; ++i){
    nNbr = _nList[i][0];
    for(unsigned j = 1 ; j<= nNbr ; ++j){  //Set to default
      _nList[i][j] = UINT_MAX;             //2^(16) - 1 = 65535
    }
    _nList[i][0] = 0;

    for(unsigned j = i + 1 ; j< _nPar ; ++j){     //Updating
      imgD = calcImgDisplacement(i, j);
      r2   = SQ(imgD[0])+SQ(imgD[1]);
      if(r2 < _rL*_rL){
        ++_nList[i][0];
        nNbr = _nList[i][0];
        _nList[i][nNbr] = j;
      }
    }
  }
}

bool checkIfUpdateNbr()
{
  double ld2;        //the longest squared displcement
  double secondld2;  //the 2nd longest squared displacement 
  double temp;
  double distance;   //how much it is possible 
                     //for two particle to get closer to each other 
  double dShell;     //shell thinkness
  secondld2 = 0.;
  dShell = _rL - _rc;


  ld2 = SQ(_r[0][0]-_rOld[0][0])+SQ(_r[0][1]-_rOld[0][1]);
  for( unsigned i = 1 ; i < _nPar ; ++i){
    temp = SQ(_r[i][0]-_rOld[i][0])+SQ(_r[i][1]-_rOld[i][1]);
    if(temp > ld2){
      secondld2 = ld2;
      ld2 = temp;
    }
    else if(temp > secondld2){
      secondld2 = temp;
    }

  }

  distance = sqrt(ld2) + sqrt(secondld2);
  if(distance >= 0.9*dShell){
    if(distance >= 0.999*dShell){
      cout << "Warning: Nbr List Algorithm is very close to break down." << endl;
      cout << "Once it breaks down, the simulation is not numerically correct or the system may crash." << endl;
      cout << "Please make _dt smaller or _rL larger to sovle this issue." << endl;
    }
    for(unsigned i = 0 ; i < _nPar ; ++i){
      for(unsigned j = 0 ; j < 2 ; ++j)
        _rOld[i][j] = _r[i][j];
    }
    return true;
  }
  else{
    return false;
  }
}

//Handling periodic boundary conditions
vector<double> calcImgDisplacement(unsigned i, unsigned j)
{
  vector<double> imgD(2, 0.);
  double d, b;

  for(unsigned k = 0 ; k < 2 ; ++k){
    d = _r[j][k] - _r[i][k];
    b = _boxSize[k];
    imgD[k] = d - floor(d/b + 0.5)*b;
  }

  return imgD;
}

//Useful functions
double dotProduct2D(vector<double> p, vector<double> q)
{
  return p[0]*q[0] + p[1]*q[1];
}

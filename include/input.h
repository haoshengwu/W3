#ifndef INPUT_H__
#define INPUT_H__


//  Define a class which contain all the input parameters:
//    
//  equlibrium_file: the equlibrium file, which is assumed to be .geqdsk file.
//  
//  others will be added step by step.

typedef struct {
  char equilibrium_file[50];
  char topology_type[20];
  //The estimation of X-point position. Current it only consider ONE X point, will be updated.
  double xpt_estimation[2];  // unit is meter.

}  InputPara;

void init_inputpara(InputPara* input);
void print_inputpara(InputPara* input);

#endif


#ifndef __INPUT_H__
#define __INPUT_H__


//  Define a class which contain all the input parameters:
//    
//  equlibrium_file: the equlibrium file, which is assumed to be .geqdsk file.
//  
//  others will be added step by step.

typedef struct {

  char equlibrium_file[50];

}  InputPara;

void init_inputpara(InputPara* input);
void print_inputpara(InputPara* input);

#endif

#include <iostream>


// argc è il counter degli argomenti che passi all'eseguibile
// argv sono gli argomenti
int main (int argc, char ** argv)
  {
    std::cout << "Input parameters from command line:" << std::endl ;
    std::cout << "You decided to pass a string of "<< argc <<" parameters: "<< std::endl ;

    for (int i=0; i<argc; i++) {
      std::cout << "Parameter"<< i << ": " << argv[i]<< std::endl;
    }
    
    return 0 ;
  }    

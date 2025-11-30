#include <iostream>

int main() {
  int num = 3;
  int *ptr = & num;

  std::cout<<"num:   "<<num<<std::endl;
  std::cout<<"ptr:   "<<ptr<<std::endl;

  int pippo = *ptr;

  std::cout<<"pippo: "<<pippo<<std::endl;
  num = 5;
  std::cout<<"num:   "<<num<<std::endl;
  std::cout<<"ptr:   "<<ptr<<std::endl;
  std::cout<<"*ptr:  "<<*ptr<<std::endl;
  std::cout<<"pippo: "<<pippo<<std::endl;

  return 0;
}

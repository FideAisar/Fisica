#include <iostream>

int main() {
  int num = 5;
  int *ptr = & num ;

  std::cout << "value   : " << num << std::endl ;
  std::cout << "address : " << &num << std::endl ;
  std::cout << "pointer : " << ptr << std::endl ;
  std::cout << "value   : " << *ptr << std::endl ;

  return 0;
}

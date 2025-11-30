#include <iostream>
using namespace std;


int main() {

  int x;
  cout << "Number: ";
  cin >> x;

  switch (x%2) {
  case 0:
    cout << "even" << endl;
    break;

  case 1:
    cout << "odd" << endl;
    break;
    
  default:
    cout << "How can I be here?" << endl;
  }

  return 0;
}



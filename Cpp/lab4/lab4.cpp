#include <iostream>
#include <cmath>

using namespace std;

class poly{

 private:
  double *v; // vettore di coefficienti
  int g; // grado del polinomio

 public:
  // costruttore senza parametri:
  poly(){
      v=NULL;
      g=-1;}
  // costruttore con grado e coefficienti:
  poly(double *vv, int gg){
      g=gg;
      v=new double[g+1];
      for(int i1=0;i1<=g;i1++) v[i1]=vv[i1];}
  // distruttore
  ~poly(){
      if(g!=-1)
        delete[] v;}

  // altre funzioni:
  void print();
  void read();
  double evaluate(double);
  void sum(const poly&, poly&);
  void product(const poly&, poly&);
  void derive(poly&);
  double root(double);
  
  void setChangeValues(double*, int);
  void setChangeValue(double, int);
  void prod(const poly&, poly&);
  void derive(const poly&);
};

// stampa polinomio
void poly::print(){
  for(int i1=0;i1<=g;i1++){
    cout << v[i1] << "*x^" << i1;
    if(i1!=g)
      cout << " + ";
  }
  cout << endl;
}

// legge polinomio da tastiera
void poly::read(){
  cout << "Grado del polinomio: ";
  cin >> g;
  cout << "Coefficienti:\n";
  v=new double[g+1];
  for(int i1=0;i1<=g;i1++)
    cin >> v[i1];
}

// calcola il valore del polinomio ad x
double poly::evaluate(double x){
  double s=0.0,pow=1.0;
  for(int i1=0;i1<=g;i1++){
      s+=v[i1]*pow;
      pow*=x;
    }
  return s;
}

// somma di polinomi
// i due polinomi sono descritti da this e p
// risultato e` creato in r
void poly::sum(const poly &p, poly &r)
{
  r.g=g;
  if(p.g>g) r.g=p.g;
  r.v=new double[r.g+1];
  for(int i1=0;i1<=r.g;i1++){
      r.v[i1]=0.0;
      if(i1<=g) r.v[i1]+=v[i1];
      if(i1<=p.g) r.v[i1]+=p.v[i1];
    }
}

void poly::prod(const poly &p, poly &r)
{
  r.v = new double[p.g + g + 1];
  for(int i=0;i<=g;i++){
    for(int j=0;j<=p.g;j++){
      r.v[i + j] += v[i]*p.v[j];
    }
  }
}

void poly::derive(const poly &r)
{
  r.v = new double[g];
  for(int i=0;i<=g;i++)
  {
    r.v[i] = i*v[i];
  }
}

void poly::setChangeValues(double* vv, int gg){
    if(g!=-1) delete[] v;
    this->g = gg;
    this->v = new double[g+1];
    for(int i=0;i<=g;i++){
        v[i] = vv[i];
    }
}
void poly::setChangeValue(double coeff, int intero){
    v[intero] = coeff;
}

int main()
{
  double coeff[]={-1.2,2.0,3.0};
  poly a(coeff,2),b,c,d,e;

  cout << "Il primo polinomio e': ";
  a.print();

  b.read();
  cout << "Il secondo polinomio e': ";
  b.print();

  a.sum(b,c);
  cout << "La somma dei due polinomi e':";
  c.print();
  
  a.prod(b,d);    
  cout << "Il prodotto dei due polinomi e':";
  d.print();
  
  a.derive(e);
  cout << "La derivata di a è':";
  e.print();
  return 0;
}
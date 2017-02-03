#include <iostream>
#include <vector>
#include <queue>
#include <math.h>
#include <assert.h> 

using namespace std;

#ifdef _WIN32
  typedef unsigned long int ulong;
  typedef unsigned short int ushort;
  typedef unsigned int uint;
#elif __APPLE__

#elif __linux__

#endif

// MATRIZ

inline int modulo(int a, int b) // OK. Sempre retorna numero nao-negativo entre [0,b-1]
{ if (b<0) { a=-a; b=-b; }
  int d=a%b;
  if (d<0) d=d+b;
  return d;
}

inline void erro(string s){
  cerr << "ERROR - " << s << endl;
}

inline float gLog2 (float val)
{
#if __ANDROID_API__ < 19
  int * const    exp_ptr = reinterpret_cast <int *> (&val);
  int            x = *exp_ptr;
  const int      log_2 = ((x >> 23) & 255) - 128;
  x &= ~(255 << 23);
  x += 127 << 23;
  *exp_ptr = x;

  val = ((-1.0f/3) * val + 2) * val - 2.0f/3;   // (1)

  return (val + log_2);
#else
  return log2(val);
#endif
}

template<typename T>
class MATRIZ {
 public:
  unsigned rows,cols,total;
  T backgv;
  vector<T> vet; 

  void fill(T e1)
  { for (unsigned i=0; i<total; i++) vet[i]=e1; 
    backgv=e1;
  }

  void resize(unsigned par_nl, unsigned par_nc)
  { rows=par_nl; cols=par_nc; total=par_nl*par_nc;  
    vet.resize(total); 
  }

  void resize(unsigned par_nl, unsigned par_nc, T e1)
  { resize(par_nl,par_nc); fill(e1); 
  }

  //<<< ctors
  MATRIZ<T>(unsigned par_nl, unsigned par_nc)
    : rows(par_nl), cols(par_nc), total(par_nl*par_nc)  
  { vet.resize(total); 
  }

  MATRIZ<T>(): MATRIZ<T>(0,0)  
  { //cout << "default ctor" << endl;
  }

  MATRIZ<T>(unsigned par_nl, unsigned par_nc, T e1): MATRIZ<T>(par_nl,par_nc)
  { fill(e1);
  }

  MATRIZ<T>(unsigned par_nl, unsigned par_nc, initializer_list<T> args)
    : MATRIZ<T>(par_nl,par_nc) 
  { auto it=begin(args);
    for (unsigned i=0; i<total; i++) {
      assert(it!=end(args));
      vet[i]=*it;
      it++;
    }
    assert(it==end(args));
    //cout << "initializer ctor" << endl << *this << endl;
  }

//   ~MATRIZ<T>() // vet vai chamar dtor automaticamente
//   { //cout << "dtor" << endl << *this << endl;
//   }
  
  //<<< copy ctor
  MATRIZ<T>(const MATRIZ<T>& a) : MATRIZ<T>(a.rows,a.cols) 
  { vet=a.vet; // copia vet.
    backgv=a.backgv;
    //cout << "copy ctor" << endl << *this << endl;
  }

  //<<< copy assign
  MATRIZ<T>& operator=(const MATRIZ<T>& a)
  { if (this == &a) {
      //cout << "self copy assign" << endl << *this << endl;
      return *this; // beware of self-assignment: x=x
    }
    if (a.total!=(*this).total) {
      vet.resize(a.total);
    }
    vet=a.vet; // copia vet
    backgv=a.backgv;
    //cout << "copy assign" << endl << *this << endl;
    return *this;
  }

  //<<< move ctor
  MATRIZ<T>(MATRIZ<T>&& a)
  { rows=a.rows; cols=a.cols; total=a.total; backgv=a.backgv;
    vet=move(a.vet);
    a.rows=a.cols=a.total=0;
    //cout << "move ctor" << endl << *this << endl;
  }

  //<<<move assign
  MATRIZ<T>& operator=(MATRIZ<T>&& a)
  { if (this == &a) {
      //cout << "self move assign" << endl << *this << endl;
      return *this; // beware of self-assignment: x=x
    }
    vet=move(a.vet);
    rows=a.rows; cols=a.cols; total=a.total; backgv=a.backgv;
    a.rows=a.cols=a.total=0;
    //cout << "move assign" << endl << *this << endl;
    return *this;
  }

  //<<<funcoes reserve, push, pop
  unsigned capacity() const 
  // Retorna numero de linhas reservadas
  { return vet.capacity() / cols;
  }

  void reserve(unsigned sz)
  // sz e' numero de linhas
  { vet.reserve(sz*cols);
  }

  void push_back(const MATRIZ<T>& a)
  { if (total>0) {
      assert(a.nc()==nc()); 
      vet.resize(total+a.total);
      for (unsigned i=0; i<a.total; i++)
        vet[i+total] = a.vet[i];
      rows += a.rows;
      total += a.total;
    } else {
      rows=a.rows; cols=a.cols; total=a.total; backgv=a.backgv;
      vet=a.vet;
    }
  }

  void pop_back(unsigned nrows=1)
  { assert(nrows>=1 && nrows<=rows);
    rows=rows-nrows;
    total=rows*cols;
    vet.resize(total);
  }

  //<<< Funcoes de acesso
  T& operator() (unsigned l, unsigned c)
  { assert(l<rows && c<cols);
    return vet[l*cols+c]; 
  }
  T& operator() (unsigned i)
  { assert(i<total);
    return vet[i]; 
  }

  T& atf(unsigned i) { // Free - Sem verificacao de indice invalido
    return vet[i];
  }
  T& ate(unsigned i) { // Modo geracao de erro se indexar fora do dominio
    if (i<total) return vet[i];
    else erro("Erro MATRIZ Indice invalido");
  }
  T& atf(unsigned l, unsigned c) { // Free - Sem verificacao de indice invalido
    return vet[l*cols+c];
  }
  T& atn(unsigned l, unsigned c) { // Modo normal com backg
    if (rows<=l || cols<=c) return backgv;
    else return vet[l*cols+c];
  }
  T& ate(unsigned l, unsigned c) { // Modo geracao de erro se indexar fora do dominio
    if (rows<=l || cols<=c)
      erro("Erro MATRIZ Indice invalido");
    return vet[l*cols+c];
  }
  T& atr(int l, int c) { // Modo replicado
    return vet[modulo(l,rows)*cols+modulo(c,cols)];
  }

//   const T& operator() (unsigned l, unsigned c) const
//   { assert(l<rows && c<cols);
//     cout << "const operator(l,c)" << endl;
//     return vet[l*cols+c]; 
//   }

//   const T& operator() (unsigned i) const
//   { assert(i<n);
//     return vet[i]; 
//   }

  unsigned nl() const { return rows; }
  unsigned nc() const { return cols; }
  unsigned n() const { return total; }
  T& backg() { return backgv; }
  const T* data() const { return vet.data(); }
};

template<typename T>
ostream& operator<< (ostream& saida, MATRIZ<T>& x)
{ for (unsigned l=0; l<x.rows; l++) {
    for (unsigned c=0; c<x.cols-1; c++) 
      saida << x(l,c) << ", ";
    saida << x(l,x.cols-1) << ";" << endl;
  }
  return saida;
}

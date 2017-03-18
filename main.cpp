#include "lib/granul.h"

void strSplit(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
}

string getSuffix(string st) 
// Devolve sufixo em minusculas "c:\ab\cd.ext => ext"
{ size_t comp=st.length();
  if (comp<=0 || comp>255){
    printf("Erro sufixo: nome arquivo invalido %s",st.c_str());
    return "";
  }

  size_t ponto=st.find_last_of('.');

  string suf;
  if (ponto==size_t(-1)) suf=""; // Arquivo sem extensao
  else {
    suf = st.substr(ponto+1);
    transform(suf.begin(), suf.end(), suf.begin(),(int (*)(int))tolower);
  }

  return suf;
}

int leChar(FILE* arq)
{ int c=fgetc(arq);
  if (c==26) c=EOF;
  if (c!=EOF) c=tolower(c);
  return c;
}

bool separador(int c)
{ if (c==',' || c==';' || c==' ' || c=='\t' || c=='\n' || c=='\r' || c=='(' || c==')' || c=='=')
    return true;
  else
    return false;
}

bool ident(int c)
{ if ( ('a'<=c && c<='z') || ('0'<=c && c<='9')
       || c=='+' || c=='-' || c=='.' || c=='_') return true;
  else return false;
}

// Comentarios podem vir:
// (1) Com //
// (2) Com /* */
// (3) Com #
string leStr(FILE* arq)
{ string st; int c;

  c=leChar(arq);
  LAB1:
    if (c==EOF) { st="eof"; goto FIM; }
    else if (separador(c)) { c=leChar(arq); goto LAB1; }
    else if (ident(c)) { st=c; c=leChar(arq); goto LAB2; }
    else if (c=='/') { c=leChar(arq); goto LAB3; }
    else if (c=='#') { c=leChar(arq); goto LAB4; }
    else { st="erro"; goto FIM; }
  LAB2:
    if (c==EOF || separador(c)) { goto FIM; }
    else if (ident(c)) { st+=c; c=leChar(arq); goto LAB2; }
    else { st="erro"; goto FIM; }
  LAB3:
    if (c=='/') { c=leChar(arq); goto LAB4; }
    else if (c=='*') { c=leChar(arq); goto LAB5; }
    else { st="erro"; goto FIM; }
  LAB4:
    if (c=='\n') { c=leChar(arq); goto LAB1; }
    else if (c==EOF) { st="eof"; goto FIM; }
    else { c=leChar(arq); goto LAB4; }
  LAB5:
    if (c=='*') { c=leChar(arq); goto LAB6; }
    else if (c!=EOF) { c=leChar(arq); goto LAB5; }
    else { st="erro"; goto FIM; }
  LAB6:
    if (c=='/') { c=leChar(arq); goto LAB1; }
    else if (c!=EOF) { c=leChar(arq); goto LAB5; }
    else { st="erro"; goto FIM; }
  FIM:
    return st;
}

void leNum(FILE* arq, BYTE& b)
{ string st=leStr(arq);
  int i;
  if (sscanf(st.c_str(),"%d",&i)!=1) erro("leNum(BYTE)");
  if (i<0) b=0;
  else if (i>255) b=255;
  else b=i;
}

void leNum(FILE* arq, int& i)
{ string st=leStr(arq);
  if (sscanf(st.c_str(),"%d",&i)!=1) erro("leNum(int)");
}

void leNum(FILE* arq, FLT& f)
{ string st=leStr(arq);
  if (st=="max" || st=="m" || st=="+max" || st=="infinito" || st=="+m") {
    f=std::numeric_limits<float>::infinity(); return;
  }
  if (st=="-max" || st=="-m" || st=="min" || st=="-infinito" || st==".") {
    f=-std::numeric_limits<float>::infinity(); return;
  }
  double d;
  if (sscanf(st.c_str(),"%lf",&d)!=1) erro("leNum(FLT)");
  f=d;
}

void leNum(FILE* arq, DBL& d)
{ string st=leStr(arq);
  if (st=="max" || st=="m" || st=="+max" || st=="infinito" || st=="+m") {
    d=std::numeric_limits<double>::infinity(); return;
  }
  if (st=="-max" || st=="-m" || st=="min" || st=="-infinito" || st==".") {
    d=-std::numeric_limits<double>::infinity(); return;
  }
  if (sscanf(st.c_str(),"%lf",&d)!=1) erro("leNum(double)");
}

void leNum(FILE* arq, bool& b)
{ string st=leStr(arq);
  if (st=="true" || st=="1" || st=="t" || st=="v" || st=="x") b=true;
  else if (st=="false" || st=="0" || st=="f" || st==".") b=false;
  else erro("leNum(bool)");
}


void le(vector<MAXLOCAL>& v, char* fileName){
  FILE* arq=fopen(fileName,"r");
  if (arq==NULL) printf("Can not open file %s\n", fileName);
  int n; leNum(arq,n); v.reserve(n);
  MAXLOCAL m; string st2;
  for (int i=0; i<n; i++) {
    leNum(arq,m.r);
    st2=leStr(arq); if (st2.size()==0) erro("Erro leitura forma"); m.forma=st2[0];
    leNum(arq,m.li);
    leNum(arq,m.ci);
    leNum(arq,m.d1);
    leNum(arq,m.d2);
    leNum(arq,m.ang);
    v.push_back(m);
  }
  fclose(arq);
}

void constroiKernel(string nome, vector<KERNEL>& ker)
{ FILE* arq=fopen(nome.c_str(),"rt");

  if (arq==NULL){
    ostringstream strStream; 
    strStream << "Erro abertura arquivo " <<  nome;
    erro(strStream.str());
    return;
  }
  string st;
  double passoAng=10; int escOitava=5; double escMin=0.5, escMax=1.0;

  st=leStr(arq); 
  while (st!="eof") {
    if (st=="passoang") { 
      leNum(arq,passoAng);
    } else if (st=="escoitava") {
      leNum(arq,escOitava);
      if (!(1<=escOitava)) erro("Erro: deve ter 1<=escOitava");
    } else if (st=="escala" || st=="escmin") {
      leNum(arq,escMin);
      if (!(0<=escMin)) erro("Erro: deve ter 0<=escMin");
    } else if (st=="a" || st=="escmax") {
      leNum(arq,escMax);
      if (!(escMin<=escMax)) erro("Erro: deve ter escMin<=escMax");
    } else if (st=="circulo") {
      double d; leNum(arq,d);

      createKernels(ker, 'c', d, d, escOitava, escMax, escMin, passoAng);

    } else if (st=="quadrado") {
      double d; leNum(arq,d);

      createKernels(ker, 's', d, d, escOitava, escMax, escMin, passoAng);

    } else if (st=="retangulo") {
      double d1,d2; leNum(arq,d1); leNum(arq,d2);

      createKernels(ker, 'r', d1, d2, escOitava, escMax, escMin, passoAng);

    } else if (st=="elipse") {
      double d1,d2; leNum(arq,d1); leNum(arq,d2);

      createKernels(ker, 'e', d1, d2, escOitava, escMax, escMin, passoAng);

    } else {
      cout << "Erro: Palavra desconhecida " << st << endl;
    }
    st=leStr(arq);
  }
  fclose(arq);
}

void kernel(int argc, char** argv){ 

  if (argc!=3) {
    printf("Kernel kernel.cfg kernel.pgm\n");
    perror("Erro: Numero de argumentos invalido");
    return;
  }

  vector<KERNEL> ker;

  constroiKernel(argv[1],ker);
  printf("Numero de kernels=%lu\n",ker.size());

  int mxH = 0, mxW = 0;
  int cntH = 0, cntW = 0;

  for (unsigned i=0; i<ker.size(); i++) {
    printf("forma=%c d1=%g d2=%g ang=%g",
            ker[i].shape,ker[i].d1,ker[i].d2,ker[i].ang);
    printf("\n");

    if(i%5 == 0){
      cntH = 0;
      mxW += cntW + 2;
      cntW = 0;
    }

    cntH += ker[i].imgG.rows + 1;

    if(mxH < cntH) {
      mxH = cntH;
    }

    if(cntW < ker[i].imgG.cols) {
      cntW = ker[i].imgG.cols;
    }
  }

  if(cntW > 0){
    mxW += cntW + 2;
  }

  Mat_<GRY> hor(mxH, mxW);

  hor = 0;

  mxH = 0, mxW = 0;
  cntH = 0, cntW = 0;

  for (unsigned i=0; i<ker.size(); i++) {
    if(i%5 == 0){
      cntH = 0;
      cntW += mxW + 2;
      mxW = 0;
    }

    for(int l = 0; l < ker[i].imgG.rows; l++){
      for(int c = 0; c < ker[i].imgG.cols; c++){
        hor.at<GRY>(cntH+l, cntW+c) = ker[i].imgG.at<GRY>(l, c);
      }
    }

    cntH += ker[i].imgG.rows + 1;

    if(mxW < ker[i].imgG.cols) {
      mxW = ker[i].imgG.cols;
    }

  }

  imwrite(argv[2], hor);   
}

void correla(int argc, char** argv){ 

  if (argc<4 || 6<argc) {
    printf("Correla imagem.pgm kernel.cfg imagem.ho1 [minCorr] [maxDist]\n");
    printf("  maxLocal com correlacao<minCorr e' eliminado\n");
    printf("  um maxLocal com distancia<=maxDist de outro maxLocal e' eliminado\n");
    printf("  Default: minCorr=0.001 maxDist=2\n");
    perror("Erro: Numero de argumentos invalido");
    return;
  }

  double minCorr=0.001; 
  if (argc>=5) minCorr = std::stof(string(argv[4]));

  int maxDist=2; 
  if (argc>=6) maxDist = std::stoi(string(argv[5]));

  printf("Lendo %s...\n",argv[1]);
  Mat_<GRY> ent; ent = imread( argv[1], 0);

  printf("Lendo %s e construindo kernels...\n", argv[2]);
  vector<KERNEL> ker;
  constroiKernel(argv[2], ker);
  printf("Numero de kernels=%lu\n",ker.size());

  printf("Achando maximos locais...\n");
  vector<MAXLOCAL> u,v;

  correlationInBatch(ent, v, ker, minCorr );

  printf("\n");

  printf("Ordenando os %lu maximos locais em ordem decrescente de correlacao...\n",v.size());
  sortMAXLOCAL(v,0,v.size()-1);

  printf("Eliminando os maximos locais muito proximos...\n");
  removeCloser(v,u,maxDist);

  printf("Imprimindo os %lu maximos locais\n",u.size());
  printTextMAXLOCAL(u, argv[3]);
}

void mostra(int argc, char** argv)
{ if (argc!=4 && argc!=5) {
    printf("Mostra imagem.pgm imagem.ho1 saida.ppm [n/d]\n");
    printf("  Mostra os maximos locais sem filtrar\n");
    printf("  n=normal d=detalhado\n");
    perror("Erro: Numero de argumentos invalido");
  }

  printf("Lendo %s...\n",argv[1]);
  Mat ent; ent = imread( argv[1]);

  printf("Lendo %s...\n",argv[2]);
  vector<MAXLOCAL> v; le(v,argv[2]);

  printf("Numero de maximos locais originais=%lu\n",v.size());

  char modo='n';
  if (argc>=5) modo=tolower(argv[4][0]);
  if (modo!='n' && modo!='d') erro("Erro modo invalido");

  printf("Tracando contornos...\n");
  Mat sai = ent.clone(); 
  if (modo=='n') printMAXLOCAL(v, sai);
  else printMAXLOCAL(v,sai, true, true);

  imwrite(argv[3], sai);
}

void filtra(int argc, char** argv)
{ if (argc<4) {
    printf("Filtra imagem.pgm imagem.ho1 saida.ppm|saida.ho2 parametros\n");
    printf("  parametros: default   minCorrDef maxInterDef (0.1 e 0.3)\n");
    printf("              circulo   minCorrCir maxInterCir\n");
    printf("              retangulo minCorrRet maxInterRet\n");
    printf("              quadrado  minCorrQua maxInterQua\n");
    printf("              elipse    minCorrEli maxInterEli\n");
    perror("Erro: Numero de argumentos invalido");
    return;
  }

  double minCorrDef=0.1;
  double minCorrCir=-1;
  double minCorrRet=-1;
  double minCorrQua=-1;
  double minCorrEli=-1;
  double maxInterDef=0.3; 
  double maxInterCir=-1; 
  double maxInterRet=-1; 
  double maxInterQua=-1; 
  double maxInterEli=-1; 
  if (argc>=5) {
    //ISTR arq; arq.arg2str(argc-4,&argv[4]);
    //string st=leStr(arq);
    std::vector<std::string> vStr;
    strSplit(string(argv[4]), ' ', vStr);
    uint count = 0;

    string st;

    while (vStr.size() > count) {
      st = vStr[count];

      if (st=="default") { 
        count++;
        minCorrDef = std::stod(vStr[count]); minCorrCir=minCorrRet=minCorrQua=minCorrEli=minCorrDef;
        count++;
        maxInterDef = std::stod(vStr[count]); maxInterCir=maxInterRet=maxInterQua=maxInterEli=maxInterDef;

      } else if (st=="circulo") { 
        count++;
        minCorrCir = std::stod(vStr[count]);
        count++;
        maxInterCir = std::stod(vStr[count]);

      } else if (st=="retangulo") { 
        count++;
        minCorrRet = std::stod(vStr[count]);
        count++;
        maxInterRet = std::stod(vStr[count]);

      } else if (st=="quadrado") { 
        count++;
        minCorrQua = std::stod(vStr[count]);
        count++;
        maxInterQua = std::stod(vStr[count]);

      } else if (st=="elipse") { 
        count++;
        minCorrEli = std::stod(vStr[count]);
        count++;
        maxInterEli = std::stod(vStr[count]);
      } else {
         printf("Erro: Palavra desconhecida %s",st.c_str());
      }
      count++;
    }
  }
  if (minCorrCir<0) minCorrCir=minCorrDef;
  if (minCorrRet<0) minCorrRet=minCorrDef;
  if (minCorrQua<0) minCorrQua=minCorrDef;
  if (minCorrEli<0) minCorrEli=minCorrDef;
  if (maxInterCir<0) maxInterCir=maxInterDef;
  if (maxInterRet<0) maxInterRet=maxInterDef;
  if (maxInterQua<0) maxInterQua=maxInterDef;
  if (maxInterEli<0) maxInterEli=maxInterDef;
  
  printf("Lendo %s...\n",argv[1]);
  Mat_<GRY> ent; ent = imread( argv[1], 0);

  printf("Lendo %s...\n",argv[2]);
  vector<MAXLOCAL> v;
  le(v,argv[2]);
  printf("Numero de maximos locais originais=%lu\n",v.size());

  printf("Filtrando os maximos locais com interseccao>maxInter...\n");
  vector<MAXLOCAL> w;
  sift(v, w, minCorrCir,maxInterCir,minCorrRet,maxInterRet,minCorrQua,maxInterQua,minCorrEli,maxInterEli, ent.size());
  printf("Numero de maximos locais apos filtro interseccao=%lu\n",w.size());

  if (getSuffix(argv[3])=="ho2") {
    printTextMAXLOCAL(w, argv[3]);
  } else {
    printf("Tracando contornos...\n");
    Mat sai; sai = imread( argv[1]);
    printMAXLOCAL(w, sai);
    imwrite(argv[3], sai);
  }
}

void relat(int argc, char** argv){ 

  if (argc!=3) {
    printf("Relat ent.ho2 sai.rel\n");
    perror("Erro: Numero de argumentos invalido");
    return;
  }

  printf("Lendo %s...\n",argv[1]);
  vector<MAXLOCAL> v;
  le(v,argv[1]);
  printf("Numero de maximos locais originais=%lu\n",v.size());

  map<double,int> circulos;
  map<double,int> retangulos; // retangulos ou quadrados
  for (int i=0; i<int(v.size()); i++) {
    MAXLOCAL m=v[i];
    if (m.forma=='c') {
      double area=cvRound(M_PI*pow ((m.d1/2), 2));
      circulos[area]++; // se nao existir, assume zero?
    } else if (m.forma=='q' || (m.forma=='r')) {
      double area=cvRound(m.d1*m.d2);
      retangulos[area]++; // se nao existir, assume zero?
    }
  }

  FILE* arq=fopen(argv[2],"wt");
  if (arq==0) printf("Erro abertura %s",argv[2]);
  fprintf(arq,"Circulos:\n");
  fprintf(arq,"%6s %6s\n","area","quant");
  for (map<double,int>::iterator it=circulos.begin(); it!=circulos.end(); it++)
    fprintf(arq,"%6.0f %6d\n",(*it).first,(*it).second);
  
  fprintf(arq,"Retangulos e quadrados:\n");
  fprintf(arq,"%6s %6s\n","area","quant");
  for (map<double,int>::iterator it=retangulos.begin(); it!=retangulos.end(); it++)
    fprintf(arq,"%6.0f %6d\n",(*it).first,(*it).second);
  fclose(arq);
}

void msgranul(int argc, char** argv) //It will find by granul correlations from MSER structs
{
  if (argc<4) {
    printf("Msgranul imagem.pgm imagem.ho1 saida.ppm|saida.ho2 parametros\n");
    printf("  parametros: default   minCorrDef maxInterDef minCorrMSER (0.1 0.3 0.9)\n");
    printf("              circulo   minCorrCir maxInterCir minCorrMSER \n");
    printf("              retangulo minCorrRet maxInterRet minCorrMSER \n");
    printf("              quadrado  minCorrQua maxInterQua minCorrMSER \n");
    printf("              elipse    minCorrEli maxInterEli minCorrMSER \n");
    perror("Erro: Numero de argumentos invalido");
    return;
  }

  double minCorrDef=0.1;
  double minCorrCir=-1;
  double minCorrRet=-1;
  double minCorrQua=-1;
  double minCorrEli=-1;
  double minCorrMSER = 0.9;

  double maxInterDef=0.3; 
  double maxInterCir=-1; 
  double maxInterRet=-1; 
  double maxInterQua=-1; 
  double maxInterEli=-1; 

  if (argc>=5) {
    std::vector<std::string> vStr;
    strSplit(string(argv[4]), ' ', vStr);
    uint count = 0;

    string st;

    while (vStr.size() > count) {
      st = vStr[count];

      if (st=="default") { 
        count++;
        minCorrDef = std::stod(vStr[count]); minCorrCir=minCorrRet=minCorrQua=minCorrEli=minCorrDef;
        count++;
        maxInterDef = std::stod(vStr[count]); maxInterCir=maxInterRet=maxInterQua=maxInterEli=maxInterDef;
        count++;
        minCorrMSER = std::stod(vStr[count]);

      } else if (st=="circulo") { 
        count++;
        minCorrCir = std::stod(vStr[count]);
        count++;
        maxInterCir = std::stod(vStr[count]);
        count++;
        minCorrMSER = std::stod(vStr[count]);

      } else if (st=="retangulo") { 
        count++;
        minCorrRet = std::stod(vStr[count]);
        count++;
        maxInterRet = std::stod(vStr[count]);
        count++;
        minCorrMSER = std::stod(vStr[count]);

      } else if (st=="quadrado") { 
        count++;
        minCorrQua = std::stod(vStr[count]);
        count++;
        maxInterQua = std::stod(vStr[count]);
        count++;
        minCorrMSER = std::stod(vStr[count]);

      } else if (st=="elipse") { 
        count++;
        minCorrEli = std::stod(vStr[count]);
        count++;
        maxInterEli = std::stod(vStr[count]);
        count++;
        minCorrMSER = std::stod(vStr[count]);
      } else {
         printf("Erro: Palavra desconhecida %s",st.c_str());
      }
      count++;
    }
  }
  if (minCorrCir<0) minCorrCir=minCorrDef;
  if (minCorrRet<0) minCorrRet=minCorrDef;
  if (minCorrQua<0) minCorrQua=minCorrDef;
  if (minCorrEli<0) minCorrEli=minCorrDef;
  if (maxInterCir<0) maxInterCir=maxInterDef;
  if (maxInterRet<0) maxInterRet=maxInterDef;
  if (maxInterQua<0) maxInterQua=maxInterDef;
  if (maxInterEli<0) maxInterEli=maxInterDef;

  printf("Lendo %s...\n",argv[1]);
  Mat_<GRY> ent; ent = imread( argv[1], 0);
  Mat_<COR> entC; entC = imread(argv[1]);

  printf("Lendo %s...\n",argv[2]);
  vector<MAXLOCAL> v;
  le(v,argv[2]);
  printf("Numero de maximos locais originais=%lu\n",v.size());

  printf("Filtrando os maximos locais com interseccao>maxInter...\n");
  vector<MAXLOCAL> w;
  sift(v, w, minCorrCir,maxInterCir,minCorrRet,maxInterRet,minCorrQua,maxInterQua,minCorrEli,maxInterEli, ent.size());
  printf("Numero de maximos locais originais=%lu\n",w.size());
  printf("Filtrando os maximos locais com MSER...\n");
  v = w;
  w.clear();
  siftMSER(entC, v, w, minCorrMSER);
  printf("Numero de maximos locais apos filtro interseccao=%lu\n",w.size());

  if (getSuffix(argv[3])=="ho2") {
    printTextMAXLOCAL(w, argv[3]);
  } else {
    printf("Tracando contornos...\n");
    Mat sai; sai = imread( argv[1]);
    printMAXLOCAL(w, sai);
    imwrite(argv[3], sai);
  }
}

void msgranul_kmeans(int argc, char** argv) //It will find by granul correlations from MSER structs
{
  if (argc<4) {
    printf("Msgranul_kmeans imagem.pgm imagem.ho1 saida.ppm|saida.ho2 parametros\n");
    printf("  parametros: default   minCorrDef maxInterDef minCorrMSER num_k (0.1 0.3 0.9 5)\n");
    printf("              circulo   minCorrCir maxInterCir minCorrMSER num_k\n");
    printf("              retangulo minCorrRet maxInterRet minCorrMSER num_k\n");
    printf("              quadrado  minCorrQua maxInterQua minCorrMSER num_k\n");
    printf("              elipse    minCorrEli maxInterEli minCorrMSER num_k\n");
    perror("Erro: Numero de argumentos invalido");
    return;
  }

  double minCorrDef=0.1;
  double minCorrCir=-1;
  double minCorrRet=-1;
  double minCorrQua=-1;
  double minCorrEli=-1;
  double minCorrMSER = 0.9;
  double num_k = 5;

  double maxInterDef=0.3; 
  double maxInterCir=-1; 
  double maxInterRet=-1; 
  double maxInterQua=-1; 
  double maxInterEli=-1; 

  if (argc>=5) {
    std::vector<std::string> vStr;
    strSplit(string(argv[4]), ' ', vStr);
    uint count = 0;

    string st;

    while (vStr.size() > count) {
      st = vStr[count];

      if (st=="default") { 
        count++;
        minCorrDef = std::stod(vStr[count]); minCorrCir=minCorrRet=minCorrQua=minCorrEli=minCorrDef;
        count++;
        maxInterDef = std::stod(vStr[count]); maxInterCir=maxInterRet=maxInterQua=maxInterEli=maxInterDef;
        count++;
        minCorrMSER = std::stod(vStr[count]);
        count++;
        num_k = std::stod(vStr[count]);

      } else if (st=="circulo") { 
        count++;
        minCorrCir = std::stod(vStr[count]);
        count++;
        maxInterCir = std::stod(vStr[count]);
        count++;
        minCorrMSER = std::stod(vStr[count]);
        count++;
        num_k = std::stod(vStr[count]);

      } else if (st=="retangulo") { 
        count++;
        minCorrRet = std::stod(vStr[count]);
        count++;
        maxInterRet = std::stod(vStr[count]);
        count++;
        minCorrMSER = std::stod(vStr[count]);
        count++;
        num_k = std::stod(vStr[count]);

      } else if (st=="quadrado") { 
        count++;
        minCorrQua = std::stod(vStr[count]);
        count++;
        maxInterQua = std::stod(vStr[count]);
        count++;
        minCorrMSER = std::stod(vStr[count]);
        count++;
        num_k = std::stod(vStr[count]);

      } else if (st=="elipse") { 
        count++;
        minCorrEli = std::stod(vStr[count]);
        count++;
        maxInterEli = std::stod(vStr[count]);
        count++;
        minCorrMSER = std::stod(vStr[count]);
        count++;
        num_k = std::stod(vStr[count]);
      } else {
         printf("Erro: Palavra desconhecida %s",st.c_str());
      }
      count++;
    }
  }
  if (minCorrCir<0) minCorrCir=minCorrDef;
  if (minCorrRet<0) minCorrRet=minCorrDef;
  if (minCorrQua<0) minCorrQua=minCorrDef;
  if (minCorrEli<0) minCorrEli=minCorrDef;
  if (maxInterCir<0) maxInterCir=maxInterDef;
  if (maxInterRet<0) maxInterRet=maxInterDef;
  if (maxInterQua<0) maxInterQua=maxInterDef;
  if (maxInterEli<0) maxInterEli=maxInterDef;

  printf("Lendo %s...\n",argv[1]);
  Mat_<GRY> ent; ent = imread( argv[1], 0);
  Mat_<COR> entC; entC = imread(argv[1]);

  printf("Lendo %s...\n",argv[2]);
  vector<MAXLOCAL> v;
  le(v,argv[2]);
  printf("Numero de maximos locais originais=%lu\n",v.size());

  printf("Filtrando os maximos locais com interseccao>maxInter...\n");
  vector<MAXLOCAL> w;
  sift(v, w, minCorrCir,maxInterCir,minCorrRet,maxInterRet,minCorrQua,maxInterQua,minCorrEli,maxInterEli, ent.size());
  printf("Numero de maximos locais originais=%lu\n",w.size());
  printf("Filtrando os maximos locais com MSER...\n");
  v = w;
  w.clear();
  siftMSER(entC, v, w, minCorrMSER);
  printf("Numero de maximos locais apos filtro interseccao=%lu\n",w.size());

  printf("Filtrando os maximos locais com Kmeans...\n");
  v = w;
  w.clear();
  chooseKmeans(entC, v, w, num_k);
  printf("Numero de maximos locais apos filtro interseccao=%lu\n",w.size());

  if (getSuffix(argv[3])=="ho2") {
    printTextMAXLOCAL(w, argv[3]);
  } else {
    printf("Tracando contornos...\n");
    Mat sai; sai = imread( argv[1]);
    printMAXLOCAL(w, sai);
    imwrite(argv[3], sai);
  }
}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< main <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
const char about[]=
"< MGranul.exe: Programas para granulometria multi-formas v1.2>\n"
"_______________________________________________________________________________\n";

const char help[]=
"Programas:\n"
"  Kernel   - Gera imagem das mascaras a partir de kernel.cfg\n"
"  Correla  - Maximos locais (.ho1) das correlacoes com mascaras multiformas\n"
"  Mostra   - Mostra maximos locais (.ho1) sem filtrar\n"
"  Filtra   - Filtra maximos locais (.ho1) e mostra\n"
"  Relat    - Le .ho2 e gera .rel\n"
"  Msgranul - Filtra maximos locais (.ho1) e mostra\n"
"  Msgranul_kmeans - Filtra maximos locais (.ho1) e mostra\n"
"...............................................................................\n";

int main(int argc, char** argv)
{
  printf(about);
  if (argc<2) {
    printf("%s",help);
    perror("Erro: Numero de argumentos invalido");
  } else {
    
    clock_t t1= (100*clock()+50)/CLOCKS_PER_SEC;
    string comando=string(argv[1]);
    transform(comando.begin(), comando.end(), comando.begin(),(int (*)(int))tolower);

    if      (comando=="kernel")   kernel(argc-1,&argv[1]);
    else if (comando=="correla")  correla(argc-1,&argv[1]);
    else if (comando=="mostra")   mostra(argc-1,&argv[1]);
    else if (comando=="filtra")   filtra(argc-1,&argv[1]);
    else if (comando=="relat")    relat(argc-1,&argv[1]);
    else if (comando == "msgranul") msgranul(argc - 1, &argv[1]);
    else if (comando == "msgranul_kmeans") msgranul_kmeans(argc - 1, &argv[1]);
    else printf("Erro: Programa inexistente %s",comando.c_str());
    t1 = ((100*clock()+50)/CLOCKS_PER_SEC)-t1;
    printf("Tempo gasto: %ld.%02ld segundos.\n",t1/100,t1%100);
  }
}

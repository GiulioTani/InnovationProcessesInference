/* bcl1: legge Piccolo solo dentro Grande, fa huffmann delle distanze
   e delle lunghezze trovate solo dentro P */
/* fa lettura disassata dei file con lo stesso nome (indipendemente
   dalla directory; cio' permette di calcolare l'autoentropia
   utilizzando una sola directory */
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <math.h>
#include <unistd.h>
#include <fcntl.h>
#include <thread>
#include <vector>
#include <mutex>
#include <queue>

#define numThreads 4
#define LUNGNOME 80
#define LUNGPATHNOME 160
#define MAX_WIN 32768
#define MIN(a, b) (a < b ? a : b)

int PICCOLO;
int GRANDE;
char *RIS;

static int c_cod[256];
static int c_ind[32768];

struct info_file
{
  unsigned short int lung;
  char *buff;
  char path_nome[LUNGPATHNOME];
  char nome[LUNGNOME];
  unsigned short int *indice;
  unsigned short int ultimo[256];
};

struct info_file *G;
int Ng;
int *lun_zip;
char *buff;
int l;
char *nome;
std::queue<int> Queue;
std::mutex lock;

int compara(const void *p1, const void *p2)
{
  if (lun_zip[*(int *)p1] < lun_zip[*(int *)p2])
  {
    return -1;
  };
  if (lun_zip[*(int *)p1] > lun_zip[*(int *)p2])
  {
    return 1;
  };
  return 0;
}

void sorta(int *permu, int n)
{
  /* n e' il numero di punti  */
  int i;
  for (i = 0; i < n; i++)
  {
    permu[i] = i;
  };
  qsort(permu, n, sizeof(int), compara);
}

void codice_cod(void)
{
  int i;

  for (i = 0; i < 8; i++)
  {
    c_cod[i] = 256 + i;
  };
  for (i = 8; i < 16; i++)
  {
    c_cod[i] = 264 + (i - 8) / 2;
  };
  for (i = 16; i < 32; i++)
  {
    c_cod[i] = 268 + (i - 16) / 4;
  };
  for (i = 32; i < 64; i++)
  {
    c_cod[i] = 272 + (i - 32) / 8;
  };
  for (i = 64; i < 128; i++)
  {
    c_cod[i] = 276 + (i - 64) / 16;
  };
  for (i = 128; i < 256; i++)
  {
    c_cod[i] = 280 + (i - 128) / 32;
  };
  return;
}

void codice_ind(void)
{
  int i, j, k;

  for (i = 0; i < 4; i++)
  {
    c_ind[i] = i;
  };
  k = 4;
  for (j = 1; j < 14; j++)
  {
    for (i = k; i < 2 * k; i++)
    {
      c_ind[i] = 2 * j + 2 * i / k;
    };
    k = 2 * k;
  };

  return;
}

int huff_length(int nn, int *freq, int *lungZ)
{
  struct nodo_H
  {
    int fre;
    int figlio;
    int ind;
  } * nod;
  int i, j, k, n, nj;
  int m1, m2, n_nod;
  int *ind;
  int max_lun_cod;
  int *lung;

  n = nn;
  ind = (int *)calloc(n, sizeof(int));
  lung = (int *)calloc(n, sizeof(int));
  nod = (struct nodo_H *)calloc(2 * n - 1, sizeof(struct nodo_H));

  j = 0;
  for (i = 0; i < nn; i++)
  {
    lungZ[i] = 0;
    if (freq[i] != 0)
    {
      ind[j] = j;
      nod[j].fre = freq[i];
      nod[j].figlio = -1;
      nod[j].ind = i;
      lung[j] = -1;
      j++;
    };
  };

  if (j == 0)
  {

    free(lung);
    free(ind);
    free(nod);

    return 0;
  };

  nj = j;
  n = nj;

  n_nod = n;
  while (n > 1)
  {
    if (nod[ind[0]].fre <= nod[ind[1]].fre)
    {
      m1 = nod[ind[0]].fre;
      m2 = nod[ind[1]].fre;
      j = 0;
      k = 1;
    }
    else
    {
      m2 = nod[ind[0]].fre;
      m1 = nod[ind[1]].fre;
      k = 0;
      j = 1;
    };
    for (i = 2; i < n; i++)
    {
      if (nod[ind[i]].fre < m1)
      {
        k = j;
        m2 = m1;
        m1 = nod[ind[i]].fre;
        j = i;
      }
      else
      {
        if (nod[ind[i]].fre < m2)
        {
          m2 = nod[ind[i]].fre;
          k = i;
        };
      };
    };
    nod[n_nod].fre = m1 + m2;
    nod[n_nod].figlio = -1;
    nod[ind[j]].figlio = n_nod;
    nod[ind[k]].figlio = n_nod;
    ind[j] = n_nod;
    ind[k] = ind[n - 1];
    n--;
    n_nod++;
  };
  nod[n_nod].figlio = -1;
  for (i = 0; i < nj; i++)
  {
    k = i;
    while (nod[k].figlio != -1)
    {
      lung[i]++;
      k = nod[k].figlio;
    };
  };

  max_lun_cod = 0;
  for (i = 0; i < nj; i++)
  {
    lungZ[nod[i].ind] = lung[i] + 1;
    if (lung[i] > max_lun_cod)
    {
      max_lun_cod = lung[i];
    };
  };
  free(lung);
  free(ind);
  free(nod);
  return max_lun_cod;
}

void worker()
{

  int i, h;
  int p, pmax, m_lung, q, qmax, q_lung, lung;
  int *freq_cod;
  int *freq_ind;
  int *len_cod;
  int *len_ind;

  freq_ind = (int *)calloc(30, sizeof(int));
  len_ind = (int *)calloc(30, sizeof(int));
  freq_cod = (int *)calloc(284, sizeof(int));
  len_cod = (int *)calloc(284, sizeof(int));
  while (true)
  {
    lock.lock();
    if (Queue.empty())
    {
      lock.unlock();
      break;
    }
    h = Queue.front();
    Queue.pop();
    lock.unlock();
    printf("%5.2f%%\r", 100. * h / Ng);
    /* azzera le frequenze 16384 margh*/
    for (i = 0; i < 284; i++)
    {
      freq_cod[i] = 0;
    };
    for (i = 0; i < 30; i++)
    {
      freq_ind[i] = 0;
    };
    i = 1;
    while (i < l)
    {
      m_lung = 1;
      q_lung = 1;
      p = G[h].ultimo[(unsigned char)buff[i]];
      q = G[h].ultimo[(unsigned char)buff[i + 1]];
      pmax = 0;
      qmax = 0;
      while (p != 0)
      {
        lung = 1;
        /* pezza a colori per disassare la ricerca nello stesso file*/

        if (strcmp((char *)G[h].nome, nome) != 0 || i != p)
        {
          while ((buff[i + lung] == G[h].buff[p + lung]) && (p + lung <= G[h].lung) && (i + lung <= l))
          {
            lung++;
          };
        };
        if (lung > m_lung)
        {
          m_lung = lung;
          pmax = p;
        };
        if (lung > 2)
        {
          lung = 1;
          while (q > pmax + 1)
          {
            lung = 1;
            if (strcmp((char *)G[h].nome, nome) != 0 || i + 1 != q)
            {
              while ((buff[i + 1 + lung] == G[h].buff[q + lung]) && (q + lung <= G[h].lung) && (i + 1 + lung <= l))
              {
                lung++;
              };
            };
            if (lung > q_lung)
            {
              q_lung = lung;
              qmax = q;
            };
            q = G[h].indice[q];
          };
        };
        p = G[h].indice[p];
      };
      /* ricerca residua su q */
      if (m_lung > 1)
      {
        lung = 1;
        while (q != 0)
        {
          lung = 1;
          /* pezza a colori per disassare la ricerca nello stesso file*/
          if (strcmp((char *)G[h].nome, nome) != 0 || i + 1 != q)
          {
            while ((buff[i + 1 + lung] == G[h].buff[q + lung]) && (q + lung <= G[h].lung) && (i + 1 + lung <= l))
            {
              lung++;
            };
          };
          if (lung > q_lung)
          {
            q_lung = lung;
            qmax = q;
          };
          q = G[h].indice[q];
        };
      };
      if (q_lung < m_lung + 1)
      {
        if (m_lung > 2)
        {
          m_lung = MIN(m_lung, 258);
          freq_cod[c_cod[m_lung - 3]]++;
          freq_ind[c_ind[G[h].lung - pmax]]++;
          i += m_lung;
        }
        else
        {
          freq_cod[(unsigned char)buff[i]]++;
          i++;
        };
      }
      else
      {
        q_lung = MIN(q_lung, 258);
        freq_cod[(unsigned char)buff[i]]++;
        freq_cod[c_cod[q_lung - 3]]++;
        freq_ind[c_ind[G[h].lung - qmax]]++;
        i += q_lung + 1;
      };
    };
    while (i <= l)
    {
      freq_cod[(unsigned char)buff[i]]++;
      i++;
    };

    huff_length(30, freq_ind, len_ind);
    huff_length(284, freq_cod, len_cod);

    /* aggiunge extra bit */
    for (i = 1; i < 6; i++)
    {
      len_cod[260 + 4 * i] += i;
      len_cod[260 + 4 * i + 1] += i;
      len_cod[260 + 4 * i + 2] += i;
      len_cod[260 + 4 * i + 3] += i;
    };

    for (i = 1; i < 14; i++)
    {
      len_ind[2 + 2 * i] += i;
      len_ind[2 + 2 * i + 1] += i;
    };

    lock.lock();
    lun_zip[h] = 0;
    for (i = 0; i < 284; i++)
    {
      lun_zip[h] += (freq_cod[i] * len_cod[i]);
    };
    for (i = 0; i < 30; i++)
    {
      lun_zip[h] += (freq_ind[i] * len_ind[i]);
    };
    lock.unlock();
  }
  free(freq_cod);
  free(len_cod);
  free(freq_ind);
  free(len_ind);
}
int calcola(/* char *buff,int l,int *dist, char *nome */)
{
  int h;
  std::vector<std::thread> workers;
  /* std::mutex lock;
  std::queue<int> Queue; */
  printf("INIZIO CALCOLA Ng = %d\n", Ng);
  Queue = std::queue<int>();
  for (h = 0; h < Ng; Queue.push(h++))
    ;
  for (auto i = 0; i < numThreads; i++)
    workers.push_back(std::thread(worker)); //,std::ref(dist),std::ref(buff),l,std::ref(nome)));
  for (auto &t : workers)
    if (t.joinable())
      t.join();
  return 1;
}

int carica_file_p(char *nome_dir_p)
{
  struct dirent **namelist;
  struct stat stbuf;
  char *tmp_nome;
  /* char *buff;
  char *nome; */
  char *nome_ris;
  int n, i, j, k, c, ifp, /* l, */ err;
  int *permu;
  FILE *fp_out;

  /* legge la directory di nome nome_dir_p; */
  tmp_nome = (char *)malloc(sizeof(char[LUNGPATHNOME]));
  nome = (char *)malloc(sizeof(char[LUNGNOME]));
  nome_ris = (char *)malloc(sizeof(char[LUNGPATHNOME]));

  //do not remove previous results
  strcpy(tmp_nome, "rm -rf ");
  strcat(tmp_nome, nome_dir_p);
  strcat(tmp_nome, RIS);
  strcat(tmp_nome, " ; mkdir ");
  strcat(tmp_nome, nome_dir_p);
  strcat(tmp_nome, RIS);
  system(tmp_nome);

  n = scandir(nome_dir_p, &namelist, 0, alphasort);
  i = 0;

  buff = (char *)calloc(PICCOLO + 1, sizeof(char));

  if (n == -1)
  {
    printf("ERRORE in lettura della dir dei file piccoli %s\n", nome_dir_p);
    exit(1);
  };

  codice_cod();
  codice_ind();
  lun_zip = (int *)calloc(Ng, sizeof(int));
  permu = (int *)calloc(Ng, sizeof(int));

  while (i <= n - 1)
  {
    strcpy(tmp_nome, nome_dir_p);
    strcat(tmp_nome, namelist[i]->d_name);
    stat(tmp_nome, &stbuf);
    err = 0;
    if (S_ISREG(stbuf.st_mode) == 1 && memcmp(namelist[i]->d_name, ".", 1) != 0)
    {
      l = (unsigned short int)MIN(PICCOLO, (int)stbuf.st_size);
      ifp = open(tmp_nome, O_RDONLY);
      if (ifp == -1)
      {
        err = 2;
      };
      k = read(ifp, buff + 1, l);
      close(ifp);
      if (k == 0)
      {
        err = 3;
      };
      if (k != l)
      {
        printf("Attenzione: file %s byte letti %d su %d \n",
               namelist[i]->d_name, k, l);
        l = k;
      };
      if (err != 0)
      {
        printf("Errore %d file %s, SALTATO\n", err, namelist[i]->d_name);
      }
      else
      {
        printf("Calcolo per %s\n", namelist[i]->d_name);
        strcpy(nome, namelist[i]->d_name);
        strcpy(nome_ris, nome_dir_p);
        strcat(nome_ris, RIS);
        strcat(nome_ris, "/Contig.");
        strcat(nome_ris, namelist[i]->d_name);
        if ((fp_out = fopen(nome_ris, "r")))
        {
          fclose(fp_out);
          continue;
        }
        calcola(/* buff,l,lun_zip,nome */);
        /* nel vettore lun_zip ci sono le contig */
        sorta(permu, Ng);
        fp_out = fopen(nome_ris, "w");
        /*c = MIN(1000,Ng);*/
        c = Ng;
        for (j = 0; j < c; j++)
        {
          fprintf(fp_out, "%-30s %f\n",
                  G[permu[j]].nome,
                  (float)lun_zip[permu[j]] / ((float)l));
        };
        fclose(fp_out);
      };
    };
    i++;
  };

  free(lun_zip);
  free(permu);
  free(buff);
  free(tmp_nome);
  free(nome);
  free(nome_ris);

  return 1;
}
int carica_file_g(char *nome_dir_g)
{
  struct dirent **namelist;
  struct stat stbuf;
  char *tmp_nome;
  int n, i, j, k, ifp, err;

  /* legge la directory di nome nome_dir_g;
     alloca un numero sufficiente di strutture info_file,
     seleziona i file regolari e ne memorizza nome, path+nome
     lunghezza */

  tmp_nome = (char *)malloc(sizeof(char[LUNGPATHNOME]));
  n = scandir(nome_dir_g, &namelist, 0, alphasort);
  i = 0;
  Ng = 0;

  if (n == -1)
  {
    printf("ERRORE in lettura della directory %s\n", nome_dir_g);
    exit(1);
  };

  /* n e' il numero di file nella directory, 
     Ng e' variabile globale che contiene il numero di testi= file regolari
     che non siano NASCOSTI (non cominciano con .)
  */

  G = (struct info_file *)calloc(n, sizeof(struct info_file));

  while (i <= n - 1)
  {
    strcpy(tmp_nome, nome_dir_g);
    strcat(tmp_nome, namelist[i]->d_name);
    stat(tmp_nome, &stbuf);
    err = 0;
    if (S_ISREG(stbuf.st_mode) == 1 && memcmp(namelist[i]->d_name, ".", 1) != 0)
    {
      G[Ng].lung = (unsigned short int)MIN(GRANDE, (int)stbuf.st_size);
      strcpy(G[Ng].path_nome, tmp_nome);
      strcpy(G[Ng].nome, namelist[i]->d_name);
      /* alloco la parte della struttura necessaria. 
        pieno di WARNIG non capisco qual e' il modo corretto */
      G[Ng].indice = (unsigned short int *)
          calloc((int)G[Ng].lung + 1, sizeof(unsigned short int));
      G[Ng].buff = (char *)
          calloc((int)G[Ng].lung + 1, sizeof(char));
      if (G[Ng].indice == NULL || G[Ng].buff == NULL)
      {
        err = 1;
      };
      /* carica il file e ne costruisce la copia a forma di indici e 
        gli indici delle ultime occorrenze dei caratteri. */
      ifp = open(tmp_nome, O_RDONLY);
      if (ifp == -1)
      {
        err = 2;
      };
      k = read(ifp, G[Ng].buff + 1, (int)G[Ng].lung);
      if (k == 0)
      {
        err = 3;
      }
      else
      {
        if (k != G[Ng].lung)
        {
          printf("Attenzione: file %s byte letti %d su %d \n",
                 namelist[i]->d_name, k, G[Ng].lung);
          G[Ng].lung = k;
        };
        for (j = 0; j < 256; j++)
        {
          G[Ng].ultimo[j] = 0;
        };
        for (j = 1; j <= (int)G[Ng].lung; j++)
        {
          G[Ng].indice[j] = G[Ng].ultimo[(unsigned char)G[Ng].buff[j]];
          G[Ng].ultimo[(unsigned char)G[Ng].buff[j]] = j;
        };
      };
      if (err > 0)
      {
        /* errore 1 (grave) non alloca 
          errore 2 (strano) non apre il file, forse protetto?
          errore 3 (non errore) il file e' vuoto
        */
        printf("Errore %d, file %s; SALTATO\n", err, namelist[i]->d_name);
      }
      else
      {
        Ng++;
      };
      close(ifp);
    };
    i++;
  };

  free(tmp_nome);
  return 1;
}

int main() // int   argc,  char *argv[] )
{
  char *direct;
  int i;
  int *numero;
  char *pic;
  char *gra;

  direct = (char *)calloc(LUNGPATHNOME, sizeof(char));

  pic = (char *)calloc(LUNGPATHNOME, sizeof(char));
  gra = (char *)calloc(LUNGPATHNOME, sizeof(char));
  RIS = (char *)calloc(LUNGPATHNOME, sizeof(char));
  numero = (int *)malloc(sizeof(int));

  printf("\nInserisci dimensione dei GRANDI\n");

  while (scanf("%d", numero) < 0)
    ;
  GRANDE = *numero;

  printf("\nInserisci dimensione dei PICCOLI\n");

  while (scanf("%d", numero) < 0)
    ;
  PICCOLO = *numero;

  if (PICCOLO > MAX_WIN)
  {
    PICCOLO = MAX_WIN;
  };
  if (GRANDE > MAX_WIN)
  {
    GRANDE = MAX_WIN;
  };

  sprintf(gra, "%d", GRANDE);
  sprintf(pic, "%d", PICCOLO);
  strcat(RIS, "ris_bcl+");
  strcat(RIS, gra);
  strcat(RIS, "-");
  strcat(RIS, pic);
  /*  strcpy(direct,"/home/dario/Tree/");*/
  printf("Inserisci directory dei file grandi\n");
  while (scanf("%s", direct) < 0)
    ;
  i = (int)strlen(direct);
  if (memcmp(direct + i, "/", 1) != 0)
  {
    strcat(direct, "/");
  };
  carica_file_g(direct);

  strcpy(direct, "");
  printf("\nCaricati i file grandi;\ninserisci directory dei file piccoli\n");

  while (scanf("%s", direct) < 0)
    ;
  i = (int)strlen(direct);
  if (memcmp(direct + i, "/", 1) != 0)
  {
    strcat(direct, "/");
  };

  carica_file_p(direct);

  free(direct);

  return 1;
}

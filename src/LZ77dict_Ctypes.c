#include "lib/LZ77dict_Ctypes.h"
#include <dirent.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#define _LMIN_ 3 // minimum word length is _LMIN_
#define MIN(a, b) (a < b ? a : b)
#define BUFF_STEP 1000

struct info_file
{
  int lung;
  char *buff;
  unsigned int *indice;
  unsigned int ultimo[256];
};
struct char_buffer
{
  char *buff;
  int lung, maxl;
};
int set_window;

int calcola(struct info_file G, struct char_buffer out);
/*Window < 0 1: search on copy, = 0 2: same fragment infinite window, > 0 3:
 * same fragment finite window*/

int LZ77dict(char *buffer, int window, char *out_p) //, int len
{
  struct info_file G;
  struct char_buffer OUT;
  int j;
  set_window = window;
  G.lung = strlen(buffer);
  OUT.lung = 0;
  OUT.maxl = G.lung * 2;
  /* alloco la parte della struttura necessaria.pieno di WARNIG non capisco qual
   * e' il modo corretto */
  G.indice = (unsigned int *)calloc((int)G.lung + 1, sizeof(unsigned int));
  G.buff = (char *)calloc((int)G.lung + 2, sizeof(char));
  OUT.buff = out_p;

  strcpy(G.buff + 1, buffer);
  for (j = 0; j < 256; j++)
  {
    G.ultimo[j] = 0;
  }
  for (j = 1; j <= (int)G.lung; j++)
  {
    G.indice[j] = G.ultimo[(unsigned char)G.buff[j]];
    G.ultimo[(unsigned char)G.buff[j]] = j;
  }
  return calcola(G, OUT);
}

int calcola(struct info_file G, struct char_buffer out)
{
  int p, pmax, m_lung, q, q_lung, lung;
  int i, window, finestra;
  int p_loop, pm_loop, q_loop, finite_window = set_window > 0;
  char *lettera;
  if (finite_window)
  { // opzione 3: finestra finita
    if (set_window < G.lung - 1)
      window = set_window;
    else
      window = G.lung - 1;
  }
  else
    window = G.lung - 1; // opzioni 1: copia e 2: finestra infinita
  i = 1;
  while (i < G.lung)
  {
    lettera = G.buff;
    m_lung = 1;
    q_lung = 1;
    p_loop = pm_loop = q_loop = 0;
    if (set_window >= 0)
    { // this should activate search on same sequence for
      // 3: finte sequance and 2: infinite sequence
      p = G.indice[i];
      q = G.indice[i + 1];
    }
    else
    { // otherwise 1: looks in a copy
      p = G.ultimo[(unsigned char)G.buff[i]];
      q = G.ultimo[(unsigned char)G.buff[i + 1]];
    }
    if (finite_window && !p)
    {
      p = G.ultimo[(unsigned char)G.buff[i]];
      p_loop = G.lung;
    }
    if (finite_window && !q)
    {
      q = G.ultimo[(unsigned char)G.buff[i + 1]];
      q_loop = G.lung;
    }
    pmax = 0;
    finestra = i - window;
    while (p != 0 &&
           (set_window < 0 ||
            p - p_loop > finestra))
    { // everywhere if I'm looking on a copy,
      // within the window otherwise.
      lung = 1;
      /* pezza a colori per disassare la ricerca nello stesso file always true
       * if using a window*/
      if (i != p)
      {
        while ((G.buff[i + lung] == G.buff[p + lung]) && (p + lung <= G.lung) &&
               (i + lung <= G.lung) && (lung < i - p + p_loop))
          lung++;
      }
      if (lung > m_lung)
      {
        m_lung = lung;
        pmax = p;
        pm_loop = p_loop;
      }
      if (lung > 2)
      {
        lung = 1;
        while (q + pm_loop - q_loop > pmax + 1)
        {
          lung = 1;
          if (i + 1 != q)
          {
            while ((G.buff[i + 1 + lung] == G.buff[q + lung]) &&
                   (q + lung <= G.lung) && (i + 1 + lung <= G.lung) &&
                   (lung < i + 1 - q + q_loop))
              lung++;
          }
          if (lung > q_lung)
            q_lung = lung;
          q = G.indice[q];
          if (finite_window && !q)
          {
            q = G.ultimo[(unsigned char)G.buff[i + 1]];
            q_loop = G.lung;
          }
        }
      }
      p = G.indice[p];
      if (finite_window && !p)
      {
        p = G.ultimo[(unsigned char)G.buff[i]];
        p_loop = G.lung;
      }
    }
    /* ricerca residua su q */
    if (m_lung > 1)
    {
      lung = 1;
      while (q != 0 && (set_window < 0 || q - q_loop > finestra + 1))
      {
        lung = 1;
        if (i + 1 != q)
        {
          while ((G.buff[i + 1 + lung] == G.buff[q + lung]) &&
                 (q + lung <= G.lung) && (i + 1 + lung <= G.lung) &&
                 (lung < i + 1 - q + q_loop))
          {
            lung++;
          }
        }
        if (lung > q_lung)
        {
          q_lung = lung;
        }
        q = G.indice[q];
        if (finite_window && !q)
        {
          q = G.ultimo[(unsigned char)G.buff[i + 1]];
          q_loop = G.lung;
        }
      }
    }
    /* output of matching sequence */
    if (q_lung < m_lung + 1)
    {
      if (m_lung >= _LMIN_)
      {
        if (out.maxl <= out.lung + m_lung + 2)
          return EXIT_FAILURE;
        strncat(out.buff, lettera + i, m_lung);
        strcat(out.buff, "\n");
        out.lung += m_lung + 1;
        i += m_lung;
      }
      else
        i++;
    }
    else
    {
      if (q_lung >= _LMIN_)
      {
        if (out.maxl <= out.lung + q_lung + 2)
          return EXIT_FAILURE;
        strncat(out.buff, lettera + i + 1, q_lung);
        strcat(out.buff, "\n");
        out.lung += q_lung + 1;
        i += q_lung + 1;
      }
      else
        i++; // to get also the non matching parts I have to do something here
    }
  }
  strcat(out.buff, "");
  return EXIT_SUCCESS;
}

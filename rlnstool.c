//rlnstool version 1.0 Oct 26, 2004
//copyright 2004 Mark G. Arnold All rights reserved


/*
NOTICE AND DISCLAIMER

No use of this software is authorized except under
this NOTICE AND DISCLAIMER. Downloading, use, copying, modification,
and/or distribution of this software in its original or modified form implies
your acceptance of this NOTICE and DISCLAIMER.

Permission to use, copy, modify and distribute this software and its
documentation, for any purpose and without fee is hereby granted, provided that
this permission notice appears prominently in supporting documentation and must
be viewed prior to use, copying, modification or distribution. You are
responsible for any modifications to the software which you make and notice that
the software has been modified must be appended to this notice prior to further
use, copying, modification or distribution.

This software is provided AS IS without warranty of any kind,
including without limitation the warranties that the software is non-infringing,
merchantable, or fit for a particular purpose, including high risk activities.
The entire risk as to the quality and performance of the software is born by
you. Should the software prove defective in any respect, you and not the
developers, nor any parties associated with the developers, assume the entire
cost of any service and repair.

This software may be subject to the Export Control Laws of the United States
of America. It is your responsibility to determine the applicable laws and
regulations and comply with them.

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <dir.h>

#include "rlnstool.h"

#ifndef ZSCALE
#define ZSCALE 16.0
#endif

#ifndef BASE
#define BASE 2.0
#endif

#ifndef RLNSPATH
#define RLNSPATH "c:/RLNS"
#endif


/* for checking against CRT in Koren p. 253-254
#define MAXMODULI 3
#define MODULUS0 2
#define MODULUS1 3
#define MODULUS2 7
#define MODULUS3 1
*/

#define M (MODULUS0*MODULUS1*MODULUS2*MODULUS3)
#define MAXPACKBITSLOW3 7
//#define MAXPACKBITS 7 (for MAXMODULI 3)
#define MAXPACKBITS 10


//also, MAXMODULI,MAXPACKBITS and M needs to be consistant with init_residue_constants

#define MIN (-(MODULUS3/2  )*MODULUS2*MODULUS1*MODULUS0)
#define MAX ( (MODULUS3/2+1)*MODULUS2*MODULUS1*MODULUS0-1)
#define MIN_BIPART (-(MODULUS2/2  )*MODULUS1*MODULUS0)
#define MAX_BIPART ( (MODULUS2/2+1)*MODULUS1*MODULUS0-1)
#define MAXPACK (1<<MAXPACKBITS)
#define MAXPACKLOW3 (1<<MAXPACKBITSLOW3)



//void testmin(){printf("%d %d %d\n",MAX,MIN,MAX-MIN+1); getchar();}

typedef struct
{
  int m[MAXMODULI];
} residuetype;

typedef struct
{
  residuetype x;
  int xsign;
} residueLNStype;

typedef residuetype bin2restype[M];
typedef int res2bintype[MAXPACK];
typedef residuetype res2restype[MAXPACK];
typedef residueLNStype resLNS2resLNStype[2*MAXPACK];

bin2restype bin2res;
res2bintype res2bin;

residuetype moduli;
residuetype bitwidth;
residuetype right_bitpos;
residuetype left_bitpos;

residuetype zero;
residuetype one;
residuetype mone;
residuetype maxres;
residuetype minres;



void init_residue_constants()
{
  int i;

  moduli.m[0] = MODULUS0;
  moduli.m[1] = MODULUS1;
  moduli.m[2] = MODULUS2;
  moduli.m[3] = MODULUS3;

  right_bitpos.m[0] = 0;
  for (i=0; i<MAXMODULI; i++)
  {
    bitwidth.m[i] = (int) ceil(log(moduli.m[i])/log(2.0));
    zero.m[i] = 0;
    one.m[i] = 1;
    mone.m[i] = moduli.m[i]-1;
    left_bitpos.m[i] =  right_bitpos.m[i] + bitwidth.m[i] - 1;
    if (i+1 != MAXMODULI)
      right_bitpos.m[i+1] =  left_bitpos.m[i] + 1;
  }
  right_bitpos.m[MAXMODULI-1] = 0;
  for (i = MAXMODULI-1 ; i>=0; i--)
  {
    printf("%d\n",i);
    left_bitpos.m[i] =  right_bitpos.m[i] + bitwidth.m[i] - 1;
    if (i != 0)
      right_bitpos.m[i-1] =  left_bitpos.m[i] + 1;
  }

  minres = zero;
  minres.m[0] = 1;

  maxres = mone;
  maxres.m[0] = 0;
}



int res2pack(residuetype r)
{
  int i,t;

  t = 0;
  for (i=0; i<MAXMODULI; i++)
  {
    t = t << bitwidth.m[i];
    t = t | r.m[i];
  }
  return t;
}



void add_residue(residuetype * s, residuetype a, residuetype b)
{
  int i;

  for (i=0; i<MAXMODULI; i++)
    (*s).m[i] = (a.m[i] + b.m[i]) % moduli.m[i];
}



void neg_residue(residuetype * r, residuetype a)
{
  int i;

  for (i=0; i<MAXMODULI; i++)
    (*r).m[i] = (moduli.m[i] - a.m[i]) % moduli.m[i];
}



void sub_residue(residuetype * s, residuetype a, residuetype b)
{
   residuetype t;

   neg_residue(&t, b);
   add_residue(s, a, t);
}



void print_residue(residuetype a)
{
  int i;

  printf("(");
  for (i=0; i<MAXMODULI; i++)
    printf("%d ",a.m[i]);
  printf(") ");
}


void fprint_residue(FILE * f, residuetype a)
{
  int i;

  fprintf(f, "(");
  for (i=0; i<MAXMODULI; i++)
    fprintf(f, "%d ",a.m[i]);
  fprintf(f, ") ");
}


int mul_inverse(int y, int modu)
{
  int x;

  for (x = 1; x<modu; x++)
  {
    if ( ((x*y)%modu) == 1)
      return x;
  }
  return -1;
}



void init_residue_tables()
{
  int cnt;
  residuetype rcnt;

  rcnt = zero;
  for (cnt = 0; cnt <= MAX; cnt++)
  {
    bin2res[cnt] = rcnt;
    res2bin[res2pack(rcnt)] = cnt;
    add_residue(&rcnt, rcnt, one);
  }

  rcnt = mone;
  for (cnt = -1; cnt >= MIN; cnt--)
  {
    bin2res[cnt+M] = rcnt;
    res2bin[res2pack(rcnt)] = cnt;
    add_residue(&rcnt, rcnt, mone);
  }
}

#define BIN2RES(x) bin2res[x+( x<0 ? M : 0 )]
#define RES2BIN(x) res2bin[res2pack(x)]


res2restype sb_lookup_tbl;
res2restype db_lookup_tbl;
resLNS2resLNStype sbdb;


void init_residue_lookup()
{
  int yi, zi;
  residuetype yres, zres;
  double y, z;

  for (zi = MIN; zi <= MAX; zi++)
  {
    z = zi / ZSCALE;
    y = log(1 + pow(BASE, z))/log(BASE);
    printf("%f sb=%f\n",z,y);
    yi = floor(y*ZSCALE + 0.5);
    zres = BIN2RES(zi);
    yres = BIN2RES(yi);
    sb_lookup_tbl[res2pack(zres)] = yres;

    sbdb[res2pack(zres)].x = yres;
    sbdb[res2pack(zres)].xsign = 0;

    if (zi != 0)
      y = log(fabs(1-pow(BASE, z)))/log(BASE);
    else
      y = MIN/ZSCALE;
    printf("%f db=%f\n",z,y);
    yi = floor(y*ZSCALE+0.5);
    zres = BIN2RES(zi);
    yres = BIN2RES(yi);
    db_lookup_tbl[res2pack(zres)] = yres;
    sbdb[MAXPACK + res2pack(zres)].x = yres;
    sbdb[MAXPACK + res2pack(zres)].xsign = zi>0;
  }
  //getchar();
}


//verilog generating routines

int check_moduli_count = 0;

void print_check_moduli_define(FILE * f, int insidemodule)
{
  int i;

  fprintf(f, "\n");
  fprintf(f,"`ifdef  MODULI");
  for (i=0; i<MAXMODULI; i++)
    fprintf(f, "_%d", moduli.m[i]);
  fprintf(f, "\n");
  fprintf(f, "`else\n");
  if (!insidemodule)
    fprintf(f, "  module moduli_mismatch_error%d;\n", check_moduli_count);
  fprintf(f, "    initial $display(%cModuli mismatch error %d%c);\n", '"', check_moduli_count , '"');
  check_moduli_count++;
  if (!insidemodule)
    fprintf(f, "  endmodule\n");
  fprintf(f, "`endif\n");

  fprintf(f, "\n");
}


void print_check_moduli_define_onehot(FILE * f, int insidemodule)
{
  int i;

  fprintf(f, "\n");
  fprintf(f,"`ifdef  MODULI_ONEHOT");
  for (i=0; i<MAXMODULI; i++)
    fprintf(f, "_%d", moduli.m[i]);
  fprintf(f, "\n");
  fprintf(f, "`else\n");
  if (!insidemodule)
    fprintf(f, "  module moduli_mismatch_error%d;\n", check_moduli_count);
  check_moduli_count++;
  fprintf(f, "    initial $display(%cModuli mismatch error%c);\n", '"', '"');
  if (!insidemodule)
    fprintf(f, "  endmodule\n");
  fprintf(f, "`endif\n");

  fprintf(f, "\n");
}


void print_add_function(FILE * f, int modu)
{
  int modu_bits;
  modu_bits = (int) ceil(log(modu)/log(2.0));
  fprintf(f,"  function [%d:0] addmod%d;\n",modu_bits-1,modu);
  //fprintf(f,"    input a,b;\n");
  fprintf(f,"    input [%d:0] a, b;\n",modu_bits-1);
  fprintf(f,"    reg [%d:0] sum, diff;\n",modu_bits);
  fprintf(f,"    begin\n");
  fprintf(f,"      //addmod%d = (a+b)%c%d;\n",modu,'%',modu);
  fprintf(f,"      sum = a + b;\n");
  fprintf(f,"      diff = sum - %d;\n", modu);
  fprintf(f,"      if (diff[%d])\n",modu_bits);
  fprintf(f,"        addmod%d = sum;\n",modu);
  fprintf(f,"      else\n");
  fprintf(f,"        addmod%d = diff;\n",modu);
  fprintf(f,"    end\n");
  fprintf(f,"  endfunction\n");
  fprintf(f,"\n");
}

void print_sub_function(FILE * f, int modu)
{
  int modu_bits;
  modu_bits = (int) ceil(log(modu)/log(2.0));
  fprintf(f,"  function [%d:0] submod%d;\n",modu_bits-1,modu);
  //fprintf(f,"    input a,b;\n");
  fprintf(f,"    input [%d:0] a, b;\n",modu_bits-1);
  fprintf(f,"    reg [%d:0] sum, diff;\n",modu_bits);
  fprintf(f,"    begin\n");
  fprintf(f,"      //submod%d = (a+%d-b)%c%d;\n",modu,modu,'%',modu);
  fprintf(f,"      diff = a - b;\n");
  fprintf(f,"      sum = diff + %d;\n", modu);
  fprintf(f,"      if (diff[%d])\n",modu_bits);
  fprintf(f,"        submod%d = sum;\n",modu);
  fprintf(f,"      else\n");
  fprintf(f,"        submod%d = diff;\n",modu);
  fprintf(f,"    end\n");
  fprintf(f,"  endfunction\n");
  fprintf(f,"\n");
}


// only good for x<2*modu
void print_mod_function(FILE * f, int modu)
{
  int modu_bits;
  modu_bits = (int) ceil(log(modu)/log(2.0));
  fprintf(f,"  function [%d:0] mod%d;\n",modu_bits-1,modu);
  fprintf(f,"    input [%d:0] x;\n",modu_bits);
  fprintf(f,"    reg [%d:0] diff;\n",modu_bits);
  fprintf(f,"    begin\n");
  fprintf(f,"      diff = x -  %d;\n", modu);
  fprintf(f,"      if (diff[%d])\n",modu_bits);
  fprintf(f,"        mod%d = x;\n",modu);
  fprintf(f,"      else\n");
  fprintf(f,"        mod%d = diff;\n",modu);
  fprintf(f,"    end\n");
  fprintf(f,"  endfunction\n");
  fprintf(f,"\n");
}


/*
    z0res.m[0] = zres.m[0];
    z0res.m[1] = zres.m[0] % MODULUS1;
    z0res.m[2] = zres.m[0] % MODULUS2;
    sub_residue(&z12res, zres, z0res);
    t2addr = (MODULUS2 + zres.m[2] - t01[zres.m[0]][zres.m[1]])%MODULUS2;
    t12res = t12[z12res.m[1]][z12res.m[2]];
    t02res = t02[zres.m[0]][t2addr];
    add_residue(&approxres, t12res, t02res);
*/

/*

void print_t01_function(FILE * f)
{
  int z0, z1;
  int w0, w1;

  w0 = left_bitpos.m[0] - right_bitpos.m[0] + 1;
  w1 = left_bitpos.m[1] - right_bitpos.m[1] + 1;

  fprintf(f,"  function [%d:0] lookup_t01;\n", left_bitpos.m[2]);
  fprintf(f,"    input [%d:0] z01;\n",left_bitpos.m[0]-right_bitpos.m[1]);
  fprintf(f,"    begin\n");
  fprintf(f,"      case(z01)\n");

  for (z0 = 0; z0 < MODULUS0; z0++)
    for (z1 = 0; z1 < MODULUS1; z1++)
      fprintf(f,"           {%d'h%x,%d'h%x} : lookup_t01 = %d; //t01[%d][%d]\n",
		     w0, z0, w1, z1, t01[z0][z1], z0, z1);

  fprintf(f,"           default : lookup_t01 = 'bx;\n");
  fprintf(f,"      endcase\n");
  fprintf(f,"    end\n");
  fprintf(f,"  endfunction\n");
  fprintf(f,"\n");
}


void print_t12_function(FILE * f)
{
  int z1, z2;
  int w1, w2;

  w1 = left_bitpos.m[1] - right_bitpos.m[1] + 1;
  w2 = left_bitpos.m[2] - right_bitpos.m[2] + 1;

  fprintf(f,"  function [%d:0] lookup_t12;\n", MAXPACKBITS-1);
  fprintf(f,"    input [%d:0] z12;\n",left_bitpos.m[1]-right_bitpos.m[2]);
  fprintf(f,"    begin\n");
  fprintf(f,"      case(z12)\n");

  for (z1 = 0; z1 < MODULUS1; z1++)
    for (z2 = 0; z2 < MODULUS2; z2++)
    {
      fprintf(f,"           {%d'h%x,%d'h%x} : lookup_t12 = 'h%x; //t12[%d][%d]=%f=",
		     w1, z1, w2, z2, res2pack(t12[z1][z2]), z1, z2, t12_double[z1][z2]);
      fprint_residue(f, t12[z1][z2]);
      fprintf(f,"\n");
    }
  fprintf(f,"           default : lookup_t12 = 'bx;\n");
  fprintf(f,"      endcase\n");
  fprintf(f,"    end\n");
  fprintf(f,"  endfunction\n");
  fprintf(f,"\n");
}


void print_t02_function(FILE * f)
{
  int z0, z2;
  int w0, w2;

  w0 = left_bitpos.m[0] - right_bitpos.m[0] + 1;
  w2 = left_bitpos.m[2] - right_bitpos.m[2] + 1;

  fprintf(f,"  function [%d:0] lookup_t02;\n", MAXPACKBITS-1);
  fprintf(f,"    input [%d:0] z02;\n", w2+w0-1);
  fprintf(f,"    begin\n");
  fprintf(f,"      case(z02)\n");

  for (z0 = 0; z0 < MODULUS0; z0++)
    for (z2 = 0; z2 < MODULUS2; z2++)
    {
      fprintf(f,"           {%d'h%x,%d'h%x} : lookup_t02 = 'h%x; //t02[%d][%d]=%f=",
		     w0, z0, w2, z2, res2pack(t02[z0][z2]), z0, z2, t02_double[z0][z2]);
      fprint_residue(f, t02[z0][z2]);
      fprintf(f,"\n");
    }
  fprintf(f,"           default : lookup_t02 = 'bx;\n");
  fprintf(f,"      endcase\n");
  fprintf(f,"    end\n");
  fprintf(f,"  endfunction\n");
  fprintf(f,"\n");
}


void print_bipart_residue()
{
  int i;
  FILE * f;
  int w0;


  f = fopen(("bipart3.v"),"w");
  fprintf(f,"// moduli [");
  for (i=0; i<MAXMODULI; i++)
    fprintf(f, "% d", moduli.m[i]);
  fprintf(f, "]\n");
  print_check_moduli_define(f,0);
  fprintf(f, "module bipart_residue(s,z);\n");
  fprintf(f, "  input [%d:0] z;\n", MAXPACKBITS-1);
  fprintf(f, "  output [%d:0] s;\n", MAXPACKBITS-1);
  fprintf(f, "  wire [%d:0] s;\n", MAXPACKBITS-1);
  fprintf(f, "  wire [%d:0] z12;\n", right_bitpos.m[0]-1);
  fprintf(f, "  wire [%d:0] z2;\n", left_bitpos.m[2]);
  fprintf(f, "  wire [%d:0] t01;\n", left_bitpos.m[2]);
  fprintf(f, "  wire [%d:0] t12;\n", MAXPACKBITS-1);
  fprintf(f, "  wire [%d:0] t02;\n", MAXPACKBITS-1);

  fprintf(f, "\n");
  for (i=0; i<MAXMODULI; i++)
    print_add_function(f,moduli.m[i]);
  for (i=1; i<MAXMODULI; i++)
    print_sub_function(f,moduli.m[i]);
  for (i=1; i<MAXMODULI; i++)
    print_mod_function(f, moduli.m[i]);

  print_t01_function(f);
  print_t12_function(f);
  print_t02_function(f);

  w0 = left_bitpos.m[0]-right_bitpos.m[0]-1;
  for (i=1; i<MAXMODULI; i++)
  {
    if (moduli.m[i] > moduli.m[0])
      fprintf(f, "  assign z12[%d:%d] = submod%d(z[%d:%d], z[%d:%d]);\n",
		 left_bitpos.m[i], right_bitpos.m[i],
		 moduli.m[i],
		 left_bitpos.m[i], right_bitpos.m[i],
		 left_bitpos.m[0], right_bitpos.m[0]);
    else if (2*moduli.m[i] > moduli.m[0])
      fprintf(f, "  assign z12[%d:%d] = submod%d(z[%d:%d], mod%d(z[%d:%d]));\n",
		 left_bitpos.m[i], right_bitpos.m[i],
		 moduli.m[i],
		 left_bitpos.m[i], right_bitpos.m[i],
		 moduli.m[i],
		 left_bitpos.m[0], right_bitpos.m[0]);
    else fprintf(f,"***error: ratio of moduli > 4***\n");
  }
  fprintf(f, "  assign t01 = lookup_t01(z[%d:%d]);\n",
		 left_bitpos.m[0], right_bitpos.m[1]);
  fprintf(f, "  assign z2 = submod%d(z[%d:%d], t01);\n",
		 moduli.m[2],
		 left_bitpos.m[2], right_bitpos.m[2]);
  fprintf(f, "  assign t02 = lookup_t02({z[%d:%d],z2});\n",
		 left_bitpos.m[0], right_bitpos.m[0]);
  fprintf(f, "  assign t12 = lookup_t12(z12);\n");

  for (i=0; i<MAXMODULI; i++)
  {
    fprintf(f, "  assign s[%d:%d] = addmod%d(t12[%d:%d], t02[%d:%d]);\n",
		 left_bitpos.m[i], right_bitpos.m[i],
		 moduli.m[i],
		 left_bitpos.m[i], right_bitpos.m[i],
		 left_bitpos.m[i], right_bitpos.m[i]);
  }


  fprintf(f, "endmodule\n");
  fprintf(f, "\n");

  fclose(f);
}
*/


void print_add_1hot_function(FILE * f, int modu)
{
  int onehotbits, s, a, b;

  onehotbits = modu;

  fprintf(f,"  function [%d:0] add1hmod%d;\n",onehotbits-1,modu);
  fprintf(f, "  input [%d:0] a,b;\n", onehotbits-1);
  fprintf(f, "  reg [%d:0] s;\n", onehotbits-1);
  fprintf(f, "\n");
  fprintf(f, "    begin\n");

    for (s=0; s<modu; s++)
    {
      fprintf(f, "      s[%d] =", s);
      for (a=0; a<modu; a++)
      {
	b = (modu + s - a) % modu;
	fprintf(f, " a[%d]&b[%d]", a, b);
	if (a != modu-1)
	  fprintf(f, " |");
      }
      fprintf(f, ";\n");
    }
  fprintf(f, "      add1hmod%d = s;\n",modu);
  fprintf(f, "    end\n");


  fprintf(f,"  endfunction\n");
  fprintf(f,"\n");
}


void print_sub_1hot_function(FILE * f, int modu)
{
  int onehotbits, s, a, b;

  onehotbits = modu;

  fprintf(f,"  function [%d:0] sub1hmod%d;\n",onehotbits-1,modu);
  fprintf(f, "  input [%d:0] a,b;\n", onehotbits-1);
  fprintf(f, "  reg [%d:0] s;\n", onehotbits-1);
  fprintf(f, "\n");
  fprintf(f, "  begin\n");

    for (s=0; s<modu; s++)
    {
      fprintf(f, "      s[%d] =", s);
      for (b=0; b<modu; b++)
      {
	a = (s + b) % modu;
	fprintf(f, " a[%d]&b[%d]", a, b);
	if (b != modu-1)
	  fprintf(f, " |");
      }
      fprintf(f, ";\n");
    }
  fprintf(f, "      sub1hmod%d = s;\n",modu);
  fprintf(f, "    end\n");


  fprintf(f,"  endfunction\n");
  fprintf(f,"\n");
}


void print_naive_add_residueLNS()
{
  int i;
  FILE * f;



  f = fopen(("nad2rlns.v"),"w");
  fprintf(f,"// moduli [");
  for (i=0; i<MAXMODULI; i++)
    fprintf(f, "% d", moduli.m[i]);
  fprintf(f, "]\n");
  print_check_moduli_define(f,0);

  fprintf(f, "//naive algorithm\n");
  fprintf(f, "module add_residueLNS(s,a,b);\n");
  fprintf(f, "  input [%d:0] a,b;\n", MAXPACKBITS);
  fprintf(f, "  output [%d:0] s;\n", MAXPACKBITS);
  fprintf(f, "  wire [%d:0] z, sbdbz;\n", MAXPACKBITS);
  fprintf(f, "\n");
  fprintf(f, "  div_residueLNS d1(z, a, b);\n");
  fprintf(f, "  sbdb_logic     t1(sbdbz, z);\n");
  fprintf(f, "  mul_residueLNS m1(s, b, sbdbz);\n");
  fprintf(f, "endmodule\n");
  fprintf(f, "\n");

  fclose(f);
}


void print_naive_add_1hot_residueLNS()
{
  int i, onehotbits;
  FILE * f;



  f = fopen(("nad1rlns.v"),"w");
  fprintf(f,"// moduli [");
  for (i=0; i<MAXMODULI; i++)
    fprintf(f, "% d", moduli.m[i]);
  fprintf(f, "] onehot\n");
  print_check_moduli_define_onehot(f,0);

  onehotbits = 0;
  for (i=0; i<MAXMODULI; i++)
    onehotbits += moduli.m[i];

  fprintf(f, "//naive algorithm\n");
  fprintf(f, "module add_1hot_residueLNS(s,a,b);\n");
  fprintf(f, "  input [%d:0] a,b;\n", onehotbits);
  fprintf(f, "  output [%d:0] s;\n", onehotbits);
  fprintf(f, "  wire [%d:0] z, sbdbz;\n", onehotbits);
  fprintf(f, "\n");
  fprintf(f, "  div_1hot_residueLNS d1(z, a, b);\n");
  fprintf(f, "  sbdb_logic_1hot     t1(sbdbz, z);\n");
  fprintf(f, "  mul_1hot_residueLNS m1(s, b, sbdbz);\n");
  fprintf(f, "endmodule\n");
  fprintf(f, "\n");

  fclose(f);
}





void print_mul_1hot_residueLNS()
{
  int i, onehotbits, bitpos, s, a, b, and_count, or_count;
  FILE * f;

  f = fopen(("mul1hlns.v"),"w");
  fprintf(f,"// onehot moduli [");
  for (i=0; i<MAXMODULI; i++)
    fprintf(f, "% d", moduli.m[i]);
  fprintf(f, "]\n");
  fprintf(f, "module mul_1hot_residueLNS(s,a,b);\n");

  onehotbits = 0;
  for (i=0; i<MAXMODULI; i++)
    onehotbits += moduli.m[i];


  fprintf(f, "  input [%d:0] a,b;\n", onehotbits);
  fprintf(f, "  output [%d:0] s;\n", onehotbits);
  fprintf(f, "  reg [%d:0] s;\n", onehotbits);
  fprintf(f, "\n");
  fprintf(f, "  always @(a or b)\n");
  fprintf(f, "    begin\n");

  and_count = 0;
  or_count = 0;

  bitpos = 0;
  i=0;
  for (i=MAXMODULI-1; i >=0; i--)
  {
    for (s=0; s<moduli.m[i]; s++)
    {
      fprintf(f, "      s[%d] =", s+bitpos);
      for (a=0; a<moduli.m[i]; a++)
      {
	and_count++;
	b = (moduli.m[i] + s - a) % moduli.m[i];
	fprintf(f, " a[%d]&b[%d]", a+bitpos, b+bitpos);
	if (a != moduli.m[i]-1)
	{
	  fprintf(f, " |");
	  or_count++;
	}
      }
      fprintf(f, "; // %d mod %d\n", s, moduli.m[i]);
    }
    bitpos += moduli.m[i];
  }
  fprintf(f, "      s[%d] = a[%d] ^ b[%d]; //sign bit\n",onehotbits, onehotbits, onehotbits);
  fprintf(f, "    end\n");

  fprintf(f, "    //%d ands, %d ors\n",and_count, or_count);
  fprintf(f, "endmodule\n");
  fprintf(f, "\n");

  fclose(f);
}



void print_div_1hot_residueLNS()
{
  int i, onehotbits, bitpos, s, a, b, and_count, or_count;
  FILE * f;

  f = fopen(("div1hlns.v"),"w");
  fprintf(f,"// onehot moduli [");
  for (i=0; i<MAXMODULI; i++)
    fprintf(f, "% d", moduli.m[i]);
  fprintf(f, "]\n");
  fprintf(f, "module div_1hot_residueLNS(s,a,b);\n");

  onehotbits = 0;
  for (i=0; i<MAXMODULI; i++)
    onehotbits += moduli.m[i];


  fprintf(f, "  input [%d:0] a,b;\n", onehotbits);
  fprintf(f, "  output [%d:0] s;\n", onehotbits);
  fprintf(f, "  reg [%d:0] s;\n", onehotbits);
  fprintf(f, "\n");
  fprintf(f, "  always @(a or b)\n");
  fprintf(f, "    begin\n");

  and_count = 0;
  or_count = 0;

  bitpos = 0;
  i=0;
  for (i=MAXMODULI-1; i >=0; i--)
  {
    for (s=0; s<moduli.m[i]; s++)
    {
      fprintf(f, "      s[%d] =", s+bitpos);
      for (b=0; b<moduli.m[i]; b++)
      {
	and_count++;
	a = (s + b) % moduli.m[i];
	fprintf(f, " a[%d]&b[%d]", a+bitpos, b+bitpos);
	if (b != moduli.m[i]-1)
	{
	  fprintf(f, " |");
	  or_count++;
	}
      }
      fprintf(f, "; // %d mod %d\n", s, moduli.m[i]);
    }
    bitpos += moduli.m[i];
  }
  fprintf(f, "      s[%d] = a[%d] ^ b[%d]; //sign bit\n",onehotbits, onehotbits, onehotbits);
  fprintf(f, "    end\n");

  fprintf(f, "    //%d ands, %d ors\n",and_count, or_count);
  fprintf(f, "endmodule\n");
  fprintf(f, "\n");

  fclose(f);
}



void print_print_residueLNS()
{
  int i,crt;
  FILE * f;



  f = fopen(("pri2rlns.v"),"w");
  fprintf(f,"// moduli [");
  for (i=0; i<MAXMODULI; i++)
    fprintf(f, "% d", moduli.m[i]);
  fprintf(f, "]\n");
  print_check_moduli_define(f,1);

  fprintf(f, "task print_residueLNS;\n");
  fprintf(f, "  input [%d:0] x;\n", MAXPACKBITS);
  for (i=0; i<MAXMODULI; i++)
    fprintf(f, "  integer dig%d;\n",i);
  fprintf(f, "  integer crt,i;\n");
  fprintf(f, "  real eps,poweps;\n");
  fprintf(f, "\n");
  fprintf(f, "    begin\n");

  fprintf(f, "      $write(%c %cb,(%c,x[%d]);\n", '"', '%', '"', MAXPACKBITS);
  for (i=0; i<MAXMODULI; i++)
  {
    fprintf(f, "      $write(%c %cd %c,x[%d:%d]);\n", '"', '%', '"',
		 left_bitpos.m[i], right_bitpos.m[i]);
    fprintf(f, "      dig%d = x[%d:%d];\n", i,
		 left_bitpos.m[i], right_bitpos.m[i]);
  }
  fprintf(f, "      $write(%c )%c);\n", '"', '"');

  fprintf(f, "      crt = (");
  for (i=0; i<MAXMODULI; i++)
  {
    crt = (M/moduli.m[i])*mul_inverse(M/moduli.m[i],moduli.m[i]);
    fprintf(f, "  %d * dig%d", crt, i);
    if (i != MAXMODULI-1)
      fprintf(f, "+");
  }
  fprintf(f, ") %c %d;\n",'%', M);

  fprintf(f, "      $write(%c %cd %c,(crt>%d)?crt-%d:crt);\n",'"','%','"',MAX,M);
  fprintf(f, "      if (crt>%d)\n",MAX);
  fprintf(f, "        begin\n");
  fprintf(f, "          crt = %d-crt;\n",M);
  fprintf(f, "          eps = %lf;\n",1.0/pow(BASE, 1.0/ZSCALE));
  fprintf(f, "        end\n");
  fprintf(f, "      else\n");
  fprintf(f, "        eps = %lf;\n",pow(BASE, 1.0/ZSCALE));
  fprintf(f, "      poweps = 1.0;\n");
  fprintf(f, "      for (i=1; i<=crt; i=i+1)\n");
  fprintf(f, "        poweps = poweps * eps;\n");
  fprintf(f, "      if (x[%d])\n",MAXPACKBITS);
  fprintf(f, "        poweps = -1.0 * poweps;\n");
  fprintf(f, "      $write(%c=%c,poweps,%c %c);\n",'"', '"', '"', '"');

  fprintf(f, "    end\n");
  fprintf(f, "endtask\n");
  fprintf(f, "\n");


  fclose(f);
}


void print_convert_residueLNS()
{
  int i,crt;
  FILE * f;



  f = fopen(("cvt2rlns.v"),"w");
  fprintf(f,"// moduli [");
  for (i=0; i<MAXMODULI; i++)
    fprintf(f, "% d", moduli.m[i]);
  fprintf(f, "]\n");
  print_check_moduli_define(f,1);



  fprintf(f, "function real residueLNS2real;\n");
  fprintf(f, "  input [%d:0] x;\n", MAXPACKBITS);
  for (i=0; i<MAXMODULI; i++)
    fprintf(f, "  integer dig%d;\n",i);
  fprintf(f, "  integer crt,i;\n");
  fprintf(f, "  real eps,poweps;\n");
  fprintf(f, "\n");
  fprintf(f, "    begin\n");

  for (i=0; i<MAXMODULI; i++)
  {
    fprintf(f, "      dig%d = x[%d:%d];\n", i,
		 left_bitpos.m[i], right_bitpos.m[i]);
  }

  fprintf(f, "      crt = (");
  for (i=0; i<MAXMODULI; i++)
  {
    crt = (M/moduli.m[i])*mul_inverse(M/moduli.m[i],moduli.m[i]);
    fprintf(f, "  %d * dig%d", crt, i);
    if (i != MAXMODULI-1)
      fprintf(f, "+");
  }
  fprintf(f, ") %c %d;\n",'%', M);

  fprintf(f, "      if (crt>%d)\n",MAX);
  fprintf(f, "        begin\n");
  fprintf(f, "          crt = %d-crt;\n",M);
  fprintf(f, "          eps = %lf;\n",1.0/pow(BASE, 1.0/ZSCALE));
  fprintf(f, "        end\n");
  fprintf(f, "      else\n");
  fprintf(f, "        eps = %lf;\n",pow(BASE, 1.0/ZSCALE));
  fprintf(f, "      poweps = 1.0;\n");
  fprintf(f, "      for (i=1; i<=crt; i=i+1)\n");
  fprintf(f, "        poweps = poweps * eps;\n");
  fprintf(f, "      if (x[%d])\n",MAXPACKBITS);
  fprintf(f, "        poweps = -1.0 * poweps;\n");
  fprintf(f, "      residueLNS2real = poweps;\n");

  fprintf(f, "    end\n");
  fprintf(f, "endfunction\n");
  fprintf(f, "\n");


  fclose(f);
}

void print_print_residueLNS_1hot()
{
  int i, onehotbits, bitpos, crt;
  FILE * f;



  f = fopen(("pri1rlns.v"),"w");
  fprintf(f,"// onehot moduli [");
  for (i=0; i<MAXMODULI; i++)
    fprintf(f, "% d", moduli.m[i]);
  fprintf(f, "]\n");
  print_check_moduli_define_onehot(f,1);
  fprintf(f, "task print_1hot_residueLNS;\n");


  onehotbits = 0;
  for (i=0; i<MAXMODULI; i++)
    onehotbits += moduli.m[i];


  fprintf(f, "  input [%d:0] x;\n", onehotbits);
  fprintf(f, "  integer digit, count;\n");
  for (i=0; i<MAXMODULI; i++)
    fprintf(f, "  integer dig%d;\n",i);
  fprintf(f, "  integer crt,i;\n");
  fprintf(f, "  real eps,poweps;\n");
  fprintf(f, "\n");
  fprintf(f, "    begin\n");

  fprintf(f, "      $write(%c %cb,(%c,x[%d]);\n", '"', '%', '"', onehotbits);
  bitpos = onehotbits;
  for (i=0; i<MAXMODULI; i++)
  {

    bitpos -= moduli.m[i];
    fprintf(f, "      count = 0;\n");
    fprintf(f, "      for (digit=0; digit<%d; digit=digit+1)\n",moduli.m[i]);
    fprintf(f, "        begin\n");
    fprintf(f, "          if (x[%d+digit])\n", bitpos);
    fprintf(f, "            begin\n");
    fprintf(f, "              count = count + 1;\n");
    fprintf(f, "              dig%d = digit;\n",i);
    fprintf(f, "              $write(%c %cd%c,digit);\n", '"', '%', '"');
    fprintf(f, "            end\n");
    fprintf(f, "        end\n");
    fprintf(f, "      if (count ==  1)\n");
    fprintf(f, "        $write(%c %c);\n", '"', '"');
    fprintf(f, "      else\n");
    fprintf(f, "        $write(%c? %c);\n", '"', '"');
  }


  fprintf(f, "      $write(%c )%c);\n", '"', '"');


  fprintf(f, "      crt = (");
  for (i=0; i<MAXMODULI; i++)
  {
    crt = (M/moduli.m[i])*mul_inverse(M/moduli.m[i],moduli.m[i]);
    fprintf(f, "  %d * dig%d", crt, i);
    if (i != MAXMODULI-1)
      fprintf(f, "+");
  }
  fprintf(f, ") %c %d;\n",'%', M);

  fprintf(f, "      $write(%c %cd %c,(crt>%d)?crt-%d:crt);\n",'"','%','"',MAX,M);
  fprintf(f, "      if (crt>%d)\n",MAX);
  fprintf(f, "        begin\n");
  fprintf(f, "          crt = %d-crt;\n",M);
  fprintf(f, "          eps = %lf;\n",1.0/pow(BASE, 1.0/ZSCALE));
  fprintf(f, "        end\n");
  fprintf(f, "      else\n");
  fprintf(f, "        eps = %lf;\n",pow(BASE, 1.0/ZSCALE));
  fprintf(f, "      poweps = 1.0;\n");
  fprintf(f, "      for (i=1; i<=crt; i=i+1)\n");
  fprintf(f, "        poweps = poweps * eps;\n");
  fprintf(f, "      if (x[%d])\n", onehotbits);
  fprintf(f, "        poweps = -1.0 * poweps;\n");
  fprintf(f, "      $write(%c=%c,poweps,%c %c);\n",'"', '"', '"', '"');

  fprintf(f, "    end\n");
  fprintf(f, "endtask\n");
  fprintf(f, "\n");

  fclose(f);
}


void print_convert_residueLNS_1hot()
{
  int i, onehotbits, bitpos, crt;
  FILE * f;



  f = fopen(("cvt1rlns.v"),"w");
  fprintf(f,"// onehot moduli [");
  for (i=0; i<MAXMODULI; i++)
    fprintf(f, "% d", moduli.m[i]);
  fprintf(f, "]\n");
  print_check_moduli_define_onehot(f,1);




  fprintf(f, "function real residueLNS1hot2real;\n");

  onehotbits = 0;
  for (i=0; i<MAXMODULI; i++)
    onehotbits += moduli.m[i];


  fprintf(f, "  input [%d:0] x;\n", onehotbits);
  fprintf(f, "  integer digit, count;\n");
  for (i=0; i<MAXMODULI; i++)
    fprintf(f, "  integer dig%d;\n",i);
  fprintf(f, "  integer crt,i;\n");
  fprintf(f, "  real eps,poweps;\n");
  fprintf(f, "\n");
  fprintf(f, "    begin\n");

  bitpos = onehotbits;
  for (i=0; i<MAXMODULI; i++)
  {

    bitpos -= moduli.m[i];
    fprintf(f, "      count = 0;\n");
    fprintf(f, "      for (digit=0; digit<%d; digit=digit+1)\n",moduli.m[i]);
    fprintf(f, "        begin\n");
    fprintf(f, "          if (x[%d+digit])\n", bitpos);
    fprintf(f, "            begin\n");
    fprintf(f, "              count = count + 1;\n");
    fprintf(f, "              dig%d = digit;\n",i);
    fprintf(f, "            end\n");
    fprintf(f, "        end\n");
  }


  fprintf(f, "      crt = (");
  for (i=0; i<MAXMODULI; i++)
  {
    crt = (M/moduli.m[i])*mul_inverse(M/moduli.m[i],moduli.m[i]);
    fprintf(f, "  %d * dig%d", crt, i);
    if (i != MAXMODULI-1)
      fprintf(f, "+");
  }
  fprintf(f, ") %c %d;\n",'%', M);

  fprintf(f, "      if (crt>%d)\n",MAX);
  fprintf(f, "        begin\n");
  fprintf(f, "          crt = %d-crt;\n",M);
  fprintf(f, "          eps = %lf;\n",1.0/pow(BASE, 1.0/ZSCALE));
  fprintf(f, "        end\n");
  fprintf(f, "      else\n");
  fprintf(f, "        eps = %lf;\n",pow(BASE, 1.0/ZSCALE));
  fprintf(f, "      poweps = 1.0;\n");
  fprintf(f, "      for (i=1; i<=crt; i=i+1)\n");
  fprintf(f, "        poweps = poweps * eps;\n");
  fprintf(f, "      if (x[%d])\n", onehotbits);
  fprintf(f, "        poweps = -1.0 * poweps;\n");
  fprintf(f, "      if (count != 1)\n");
  fprintf(f, "        poweps = 'bx;\n");
  fprintf(f, "      residueLNS1hot2real = poweps;\n");
  fprintf(f, "    end\n");
  fprintf(f, "endfunction\n");
  fprintf(f, "\n");


  fclose(f);
}




void print_naive_sbdb()
{
  int i, z;

  FILE * f;

  //f = fopen(("resLNS.v"),"w");
  f = fopen(("nsd2rlns.v"),"w");
  fprintf(f,"// moduli [");
  for (i=0; i<MAXMODULI; i++)
    fprintf(f, "% d", moduli.m[i]);
  fprintf(f, "]\n");
  print_check_moduli_define(f,0);

  fprintf(f,"// naive algorithm\n");
  fprintf(f, "module sbdb_logic(sbdb, z);\n");
  fprintf(f, "  input [%d:0] z;\n", MAXPACKBITS);
  fprintf(f, "  output [%d:0] sbdb;\n", MAXPACKBITS);
  fprintf(f, "  reg [%d:0] sbdb;\n", MAXPACKBITS);
  fprintf(f, "\n");
  fprintf(f, "  always @(z)\n");
  fprintf(f, "    begin\n");
  fprintf(f, "      case(z)\n");
  fprintf(f, "\n");

  for (z = 0; z < /*  2* only sb for jack*/ 2* MAXPACK; z++)
  {
    //printf("[%3d] = {%d,4'h%04x}; //", z,sbdb[z].xsign,res2pack(sbdb[z].x));
    fprintf(f, "          %3d : sbdb = {%d,%d'h%x}; //", z,sbdb[z].xsign,MAXPACKBITS,res2pack(sbdb[z].x));
    fprintf(f, "\t %d,",z >> (MAXPACKBITS));
    fprint_residue(f, BIN2RES(res2bin[((1 << MAXPACKBITS)-1) & z]));
    fprintf(f, " --> %d,",sbdb[z].xsign);
    fprint_residue(f, sbdb[z].x);
    fprintf(f, "\n");
  }
  fprintf(f, "      endcase\n");
  fprintf(f, "    end\n");
  fprintf(f, "endmodule\n");
  fprintf(f, "\n");

  fclose(f);
}



  int sbdblow_pos[2*MAXPACKLOW3];
  int sbdblow_neg[2*MAXPACKLOW3];
  int base_extend[2*MAXPACKLOW3];
  int sbdblowpack_pos[MAXPACKLOW3];
  int sbdblowpack_neg[MAXPACKLOW3];
  int myM = MODULUS0*MODULUS1*MODULUS2;

void init_sbdb_lowbits()
{
  int result;
  int i, z, zi, mask, zsign, zaddr,high_mrc;
  residuetype zr,zrn;
  int zn, zaddr2n, zaddr2, testpos,testneg;

  printf("%d",2*MAXPACKLOW3);

  for (z = 0; z< 2*MAXPACKLOW3; z++)
    sbdblow_pos[z] = -1;
  for (z = 0; z< 2*MAXPACKLOW3; z++)
    sbdblow_neg[z] = -1;



  for (z = 0; z< 2*MAXPACK; z++)
  {
    zi = res2bin[((1 << MAXPACKBITS)-1) & z];
    if ((zi < myM) && (( (z& ((1<<MAXPACKBITS)-1)) ==0)|(zi > 0)) )
      sbdblow_pos[z>>bitwidth.m[3]] = res2pack(sbdb[z].x);

  }

  for (z = 0; z< 2*MAXPACK; z++)
  {
    zi = res2bin[((1 << MAXPACKBITS)-1) & z];
    if ((zi >= -myM)&&(zi < 0))
      sbdblow_neg[z>>bitwidth.m[3]] = res2pack(BIN2RES(RES2BIN(sbdb[z].x)-zi));
  }

  for (z = 0; z< MAXPACKLOW3; z++)
  {
      sbdblowpack_pos[z] = sbdblow_pos[z<<1];
      sbdblowpack_neg[z] = sbdblow_neg[z<<1];
  }

  if (0)
  for (zi = -myM; zi < myM; zi++)
  {
      zr = BIN2RES(zi);
      z = res2pack(zr);
      zrn = BIN2RES(-zi);
      zn = res2pack(zrn);
      //if (zr.m[2]!=0)
      {
	zaddr = z >> bitwidth.m[3];
	zaddr2 = zaddr >> 1;
	zaddr2n = zn >> (bitwidth.m[3]+1);
	if (zaddr & 1)
	{
	  testneg = sbdblowpack_pos[zaddr2n];
	  testpos = sbdblowpack_neg[zaddr2n];
	}
	else
	{
	  testpos = sbdblowpack_pos[zaddr2];
	  testneg = sbdblowpack_neg[zaddr2];
	}
	/*
	if ((testpos == sbdblow_pos[zaddr])&&
	    (testneg == sbdblow_neg[zaddr]))
	  printf(" ok %d",zi);
	else
	{
	  printf("bad %d test p=%d n=%d  true p=%d n=%d\n",
		  zi,testpos,testneg,sbdblow_pos[zaddr],sbdblow_neg[zaddr]);
	  print_residue(zr);
	  print_residue(zrn);
	}
	getchar();
	*/
      }
  }

  for (z = 0; z< 2*MAXPACKLOW3; z++)
    base_extend[z] = -1;

  for (z = 0; z< 2*MAXPACK; z++)
  {
    zi = res2bin[((1 << MAXPACKBITS)-1) & z];
    if ((zi < myM) && (( (z& ((1<<MAXPACKBITS)-1)) ==0)|(zi > 0)) )
       base_extend [z>>bitwidth.m[3]] = zi%moduli.m[3];
  }

}


void print_mulinv_function(FILE * f, int myM, int modu)
{
  int modu_bits, i;
  modu_bits = (int) ceil(log(modu)/log(2.0));
  fprintf(f,"  function [%d:0] mulinv%d_%d;\n",modu_bits-1,myM,modu);
  fprintf(f,"    input [%d:0] x;\n",modu_bits-1);
  fprintf(f,"    begin\n");
  fprintf(f,"    case(x)\n");
  for (i=0; i<modu; i++)
  {
  fprintf(f,"      %d: mulinv%d_%d = %d;\n",i,myM,modu,
		      (mul_inverse(myM, modu) * i)%modu);
  }
  fprintf(f,"      default: mulinv%d_%d = 'bx;\n",myM,modu);
  fprintf(f,"    endcase\n");
  fprintf(f,"    end\n");
  fprintf(f,"  endfunction\n");
  fprintf(f,"\n");
}


void print_mulinvsel_function(FILE * f, int myM, int modu)
{
  int modu_bits, i, high_mrc, select;
  modu_bits = (int) ceil(log(modu)/log(2.0));
  fprintf(f,"  function [1:0] mulinv%d_%dsel;\n",myM,modu);
  fprintf(f,"    input [%d:0] x;\n",modu_bits-1);
  fprintf(f,"    begin\n");
  fprintf(f,"    case(x)\n");
  for (i=0; i<modu; i++)
  {
  high_mrc = (mul_inverse(myM, modu) * i)%modu;
  if (high_mrc == 0)
    select = 0;
  else if (high_mrc == moduli.m[3]-1 )
    select = 1;
  else if ((high_mrc > 0) & (high_mrc <=moduli.m[3]/2) )
    select = 2;
  else
    select = 3;
  fprintf(f,"      %d: mulinv%d_%dsel = %d; //%d\n",i,myM,modu,select,
		      (mul_inverse(myM, modu) * i)%modu);
  }
  fprintf(f,"      default: mulinv%d_%dsel = 'bx;\n",myM,modu);
  fprintf(f,"    endcase\n");
  fprintf(f,"    end\n");
  fprintf(f,"  endfunction\n");
  fprintf(f,"\n");
}


/*****/
void print_1hot2pack_function(FILE * f)
{
  int i, j, k, onehotbits;
  int bitpos, firsttime;
  int resbitpos;

  onehotbits = 0;
  for (i=0; i<MAXMODULI; i++)
    onehotbits += moduli.m[i];

  fprintf(f,"  function [%d:0] rlns1h2pack;\n", MAXPACKBITS);
  fprintf(f,"    input [%d:0] x;\n",onehotbits);
  fprintf(f,"    reg [%d:0] res;\n", MAXPACKBITS);
  fprintf(f,"    begin\n");



  resbitpos = 0;
  bitpos = 0;
  for (i=MAXMODULI-1; i >=0; i--)
  {
    //leftbitpos = bitpos + moduli.m[i] - 1;

    //numpackbits = (int) ceil(log(moduli.m[i])/log(2.0));
    for (k=1; k<moduli.m[i]; k=k<<1)
    {
      //fprintf(f, "res%d[%d] =",i,k);
      fprintf(f, "      res[%d] = ", resbitpos);
      firsttime = 1;
      for (j=0; j<moduli.m[i]; j++)
      {
	if (j&k)
	{
	  if (firsttime)
	    firsttime = 0;
	  else
	    fprintf(f," | ");
	  fprintf(f,"x[%d]",j+bitpos);
	}
      }
      fprintf(f,";\n");
      resbitpos++;
    }
    bitpos += moduli.m[i];
  }

  fprintf(f,"      res[%d] = x[%d];\n", MAXPACKBITS,onehotbits);
  fprintf(f,"      rlns1h2pack = res;\n");
  fprintf(f,"    end\n");
  fprintf(f,"  endfunction\n");
  fprintf(f,"\n");
}


void print_pack21hot_function(FILE * f)
{
  int i, j, k, onehotbits;
  int bitpos, firsttime;
  int resbitpos;

  onehotbits = 0;
  for (i=0; i<MAXMODULI; i++)
    onehotbits += moduli.m[i];

  fprintf(f,"  function [%d:0] rlnspack21h;\n", onehotbits);
  fprintf(f,"    input [%d:0] x;\n", MAXPACKBITS);
  fprintf(f,"    reg [%d:0] res;\n", onehotbits);
  fprintf(f,"    begin\n");



  resbitpos = 0;
  for (i=MAXMODULI-1; i >=0; i--)
  {
    for (j=0; j<moduli.m[i]; j++)
    {
      //fprintf(f, "res%d[%d] =",i,k);
      fprintf(f, "      res[%d] = ", resbitpos);
      firsttime = 1;
      bitpos = right_bitpos.m[i];
      for (k=1; k<moduli.m[i]; k=k<<1)
      {
	if (firsttime)
	  firsttime = 0;
	else
	  fprintf(f," & ");

	if (j&k)
	  fprintf(f,"x[%d]",bitpos);
	else
	  fprintf(f,"(~x[%d])",bitpos);
	bitpos++;
      }
      fprintf(f,";\n");
      resbitpos++;
    }
  }

  fprintf(f,"      res[%d] = x[%d];\n",onehotbits, MAXPACKBITS);
  fprintf(f,"      rlnspack21h = res;\n");
  fprintf(f,"    end\n");
  fprintf(f,"  endfunction\n");
  fprintf(f,"\n");
}


/*****/



void print_mulinvsel_1hot_function(FILE * f, int myM, int modu)
{
  int modu_bits, i, high_mrc, select;
  modu_bits = (int) ceil(log(modu)/log(2.0));
  fprintf(f,"  function [1:0] mulinv1h%d_%dsel;\n",myM,modu);
  fprintf(f,"    input [%d:0] x;\n",modu-1);
  fprintf(f,"    begin\n");
  fprintf(f,"    casex(x)\n");
  for (i=0; i<modu; i++)
  {
  high_mrc = (mul_inverse(myM, modu) * i)%modu;
  if (high_mrc == 0)
    select = 0;
  else if (high_mrc == moduli.m[3]-1 )
    select = 1;
  else if ((high_mrc > 0) & (high_mrc <=moduli.m[3]/2) )
    select = 2;
  else
    select = 3;
  fprintf(f,"      %d: mulinv1h%d_%dsel = %d; //%d\n",1<<i,myM,modu,select,
		      (mul_inverse(myM, modu) * i)%modu);
  }
//  fprintf(f,"      default: mulinv%d_%dsel = 'bx;\n",myM,modu);
  fprintf(f,"    endcase\n");
  fprintf(f,"    end\n");
  fprintf(f,"  endfunction\n");
  fprintf(f,"\n");
}


void print_baseextend_function(FILE * f)
{
  int modu_bits,i,in_bits;
  in_bits = (int) ceil(log(MAXPACKLOW3)/log(2.0));
  modu_bits = (int) ceil(log(moduli.m[3])/log(2.0));
  fprintf(f,"  function [%d:0] baseextend;\n", modu_bits-1);
  fprintf(f,"    input [%d:0] x;\n",in_bits);
  fprintf(f,"    begin\n");
  fprintf(f,"    casex(x[%d:0])\n",in_bits-1);
  for (i=0; i<MAXPACKLOW3; i++)
  {
  if (base_extend[i] != -1)
  fprintf(f,"      %d'h%02x: baseextend = %d;\n",in_bits,i,base_extend[i]);
  //else
  //fprintf(f,"      %d'h%02x: baseextend = %d'bx;\n",in_bits,i,in_bits);
  }
  fprintf(f,"      default: baseextend = 'bx;\n");
  fprintf(f,"    endcase\n");
  fprintf(f,"    end\n");
  fprintf(f,"  endfunction\n");
  fprintf(f,"\n");
}



void print_baseextend_1hot_function(FILE * f)
{
  int modu_bits,i,in_bits;
  unsigned long d0,d1,d2;
  unsigned long mask0,mask1,mask2;
  unsigned long onehot0, onehot1, onehot2, onehot;

  in_bits = moduli.m[0]+moduli.m[1]+moduli.m[2];
  modu_bits = moduli.m[3];

  mask2 = (1 << bitwidth.m[2]) - 1;
  mask1 = (1 << bitwidth.m[1]) - 1;
  mask0 = (1 << bitwidth.m[0]) - 1;


  fprintf(f,"  function [%d:0] baseextend1h;\n", modu_bits-1);
  fprintf(f,"    input [%d:0] x;\n",in_bits);
  fprintf(f,"    begin\n");
  fprintf(f,"    casex(x[%d:0])\n",in_bits-1);
  for (i=0; i<MAXPACKLOW3; i++)
  {
    d2 = i & mask2;
    d1 = (i >> bitwidth.m[2]) & mask1;
    d0 = (i >> (bitwidth.m[2]+bitwidth.m[1])) & mask0;

    onehot2 = 1 << d2;
    onehot1 = 1 << (d1 + moduli.m[2]);
    onehot0 = 1 << (d0 + moduli.m[1] + moduli.m[2]);
    onehot = onehot0 | onehot1 | onehot2;

    if (base_extend[i] != -1)
      fprintf(f,"      %d'h%lx: baseextend1h = 'h%x; //%d\n",in_bits,onehot,1<<base_extend[i],base_extend[i]);
  }
  fprintf(f,"      default: baseextend1h = 'bx;\n");
  fprintf(f,"    endcase\n");
  fprintf(f,"    end\n");
  fprintf(f,"  endfunction\n");
  fprintf(f,"\n");
}



//    high_mrc = (moduli.m[3] + zr.m[3] - base_extend[zaddr])%moduli.m[3];
//    high_mrc = (mul_inverse(myM, moduli.m[3]) * high_mrc)%moduli.m[3];


void print_sbdb_posneg_function(FILE * f)
{
  int modu_bits,i;
  modu_bits = (int) ceil(log(moduli.m[3])/log(2.0));
  fprintf(f,"  function [%d:0] sbdb_posneg;\n", 2*MAXPACKBITS+1);
  fprintf(f,"    input [%d:0] z;\n",MAXPACKBITS);
  fprintf(f,"    input odd;\n");
  fprintf(f,"    reg [%d:0] sbdbpos,sbdbneg;\n",MAXPACKBITS);
  fprintf(f,"    begin\n");
  fprintf(f,"    case(z>>%d)\n",modu_bits+1);
  for (i=0; i<MAXPACKLOW3; i++)
  {
  if (sbdblowpack_pos[i] != -1)
  fprintf(f,"      'h%04x: sbdbpos = 'h%04x;\n",i,sbdblowpack_pos[i]);
  }
  fprintf(f,"      default: sbdbpos = 'bx;\n");
  fprintf(f,"    endcase\n");
  fprintf(f,"   \n");

  fprintf(f,"    case(z>>%d)\n",modu_bits+1);
  for (i=0; i<MAXPACKLOW3; i++)
  {
  if (sbdblowpack_neg[i] != -1)
  fprintf(f,"      'h%04x: sbdbneg = 'h%04x;\n",i,sbdblowpack_neg[i]);
  }
  fprintf(f,"      default: sbdbneg = 'bx;\n");
  fprintf(f,"    endcase\n");

//  fprintf(f,"    if(z[%d])\n",modu_bits);
  fprintf(f,"    if(odd)\n");
  fprintf(f,"      sbdb_posneg = {sbdbneg,sbdbpos};\n");
  fprintf(f,"    else\n");
  fprintf(f,"      sbdb_posneg = {sbdbpos,sbdbneg};\n");

  fprintf(f,"    end\n");
  fprintf(f,"  endfunction\n");
  fprintf(f,"\n");
}




void print_check_relerr_1hot(FILE * f, char op)
{
  int i,onehotbits;

  onehotbits = 0;
  for (i=0; i<MAXMODULI; i++)
    onehotbits += moduli.m[i];

  fprintf(f,"\n");
  fprintf(f, "  `ifdef CHECK_RELERR\n");
  fprintf(f, "\n");
  fprintf(f, "   `include %ccvt1rlns.v%c\n",'"','"');
  fprintf(f, "   real reala, realb, reals, relerr;\n");
  fprintf(f, "   reg [%d:0] olds;\n",onehotbits);
  fprintf(f, "\n");
  fprintf(f, "   always @(s)\n");
  fprintf(f, "     begin\n");
  fprintf(f, "       olds = 'bz;\n");
  fprintf(f, "       while (olds !== s)\n");
  fprintf(f, "        begin\n");
  fprintf(f, " 	       #0; olds = s;\n");
  fprintf(f, " 	     end\n");
  fprintf(f, "       reala = residueLNS1hot2real(a);\n");
  fprintf(f, "       realb = residueLNS1hot2real(b);\n");
  fprintf(f, "       reals = residueLNS1hot2real(s);\n");
  fprintf(f, "       relerr = (reala %c realb-reals)/(reala %c realb);\n",op,op);
  fprintf(f, "       if (relerr < 0.0)\n");
  fprintf(f, "         relerr = -1.0 * relerr;\n");
  fprintf(f, "       if (relerr > `CHECK_RELERR)\n");
  fprintf(f, "        begin\n");
  fprintf(f, "          $write(%ca=%c,reala);\n",'"','"');
  fprintf(f, "          $write(%c (%ch) %c,a);\n",'"','%','"');
  fprintf(f, "          $write(%c b=%c,realb);\n",'"','"');
  fprintf(f, "          $write(%c (%ch) %c,b);\n",'"','%','"');
  fprintf(f, "          $write(%c a%cb=%c,reala %c realb);\n",'"',op,'"',op);
  fprintf(f, "          $write(%c vs. s=%c,reals);\n",'"','"');
  fprintf(f, "          $write(%c (%ch) %c,s);\n",'"','%','"');
  fprintf(f, "          $write(%c relerr=%c,relerr);\n",'"','"');
  fprintf(f, "          $display(%c time= %cd%c,$time);\n",'"','%','"');
  fprintf(f, "        end\n");
  fprintf(f, "     end\n");
  fprintf(f, "  `endif\n");
  fprintf(f, "\n");

}





void print_mul_1hotter_residueLNS()
{
  int i, onehotbits, bitpos, leftbitpos;
  FILE * f;

  f = fopen(("mul1rlns.v"),"w");
  fprintf(f,"// onehot moduli [");
  for (i=0; i<MAXMODULI; i++)
    fprintf(f, "% d", moduli.m[i]);
  fprintf(f, "]\n");
  print_check_moduli_define_onehot(f,0);


  fprintf(f, "module mul_1hot_residueLNS(s,a,b);\n");

  onehotbits = 0;
  for (i=0; i<MAXMODULI; i++)
    onehotbits += moduli.m[i];


  fprintf(f, "  input [%d:0] a,b;\n", onehotbits);
  fprintf(f, "  output [%d:0] s;\n", onehotbits);
  fprintf(f, "  reg [%d:0] s;\n", onehotbits);
  fprintf(f, "\n");

  print_check_relerr_1hot(f, '*');


  for (i=0; i<MAXMODULI; i++)
    print_add_1hot_function(f,moduli.m[i]);
  fprintf(f, "\n");
  fprintf(f, "  always @(a or b)\n");
  fprintf(f, "    begin\n");


  bitpos = 0;
  for (i=MAXMODULI-1; i >=0; i--)
  {
    leftbitpos = bitpos + moduli.m[i] - 1;
    fprintf(f, "      s[%d:%d] = add1hmod%d(a[%d:%d], b[%d:%d]);\n",
		       leftbitpos, bitpos,
		       moduli.m[i],
		       leftbitpos, bitpos,
		       leftbitpos, bitpos);

    bitpos += moduli.m[i];
  }
  fprintf(f, "      s[%d] = a[%d] ^ b[%d]; //sign bit\n",onehotbits, onehotbits, onehotbits);
  fprintf(f, "    end\n");
  fprintf(f, "endmodule\n");
  fprintf(f, "\n");

  fclose(f);
}


void print_div_1hotter_residueLNS()
{
  int i, onehotbits, bitpos, leftbitpos;
  FILE * f;

  f = fopen(("div1rlns.v"),"w");
  fprintf(f,"// onehot moduli [");
  for (i=0; i<MAXMODULI; i++)
    fprintf(f, "% d", moduli.m[i]);
  fprintf(f, "]\n");
  print_check_moduli_define_onehot(f,0);



  fprintf(f, "module div_1hot_residueLNS(s,a,b);\n");

  onehotbits = 0;
  for (i=0; i<MAXMODULI; i++)
    onehotbits += moduli.m[i];


  fprintf(f, "  input [%d:0] a,b;\n", onehotbits);
  fprintf(f, "  output [%d:0] s;\n", onehotbits);
  fprintf(f, "  reg [%d:0] s;\n", onehotbits);
  fprintf(f, "\n");

  print_check_relerr_1hot(f, '/');


  for (i=0; i<MAXMODULI; i++)
    print_sub_1hot_function(f,moduli.m[i]);
  fprintf(f, "\n");
  fprintf(f, "  always @(a or b)\n");
  fprintf(f, "    begin\n");


  bitpos = 0;
  for (i=MAXMODULI-1; i >=0; i--)
  {
    leftbitpos = bitpos + moduli.m[i] - 1;
    fprintf(f, "      s[%d:%d] = sub1hmod%d(a[%d:%d], b[%d:%d]);\n",
		       leftbitpos, bitpos,
		       moduli.m[i],
		       leftbitpos, bitpos,
		       leftbitpos, bitpos);

    bitpos += moduli.m[i];
  }
  fprintf(f, "      s[%d] = a[%d] ^ b[%d]; //sign bit\n",onehotbits, onehotbits, onehotbits);
  fprintf(f, "    end\n");
  fprintf(f, "endmodule\n");
  fprintf(f, "\n");

  fclose(f);
}


void print_check_relerr(FILE * f, char op)
{
  fprintf(f,"\n");
  fprintf(f, "  `ifdef CHECK_RELERR\n");
  fprintf(f, "\n");
  fprintf(f, "   `include %ccvt2rlns.v%c\n",'"','"');
  fprintf(f, "   real reala, realb, reals, relerr;\n");
  fprintf(f, "   reg [%d:0] olds;\n",MAXPACKBITS);
  fprintf(f, "\n");
  fprintf(f, "   always @(s)\n");
  fprintf(f, "     begin\n");
  fprintf(f, "       olds = 'bz;\n");
  fprintf(f, "       while (olds !== s)\n");
  fprintf(f, "        begin\n");
  fprintf(f, " 	       #0; olds = s;\n");
  fprintf(f, " 	     end\n");
  fprintf(f, "       reala = residueLNS2real(a);\n");
  fprintf(f, "       realb = residueLNS2real(b);\n");
  fprintf(f, "       reals = residueLNS2real(s);\n");
  fprintf(f, "       relerr = (reala %c realb-reals)/(reala %c realb);\n",op,op);
  fprintf(f, "       if (relerr < 0.0)\n");
  fprintf(f, "         relerr = -1.0 * relerr;\n");
  fprintf(f, "       if (relerr > `CHECK_RELERR)\n");
  fprintf(f, "        begin\n");
  fprintf(f, "          $write(%ca=%c,reala);\n",'"','"');
  fprintf(f, "          $write(%c (%ch) %c,a);\n",'"','%','"');
  fprintf(f, "          $write(%c b=%c,realb);\n",'"','"');
  fprintf(f, "          $write(%c (%ch) %c,b);\n",'"','%','"');
  fprintf(f, "          $write(%c a%cb=%c,reala %c realb);\n",'"',op,'"',op);
  fprintf(f, "          $write(%c vs. s=%c,reals);\n",'"','"');
  fprintf(f, "          $write(%c (%ch) %c,s);\n",'"','%','"');
  fprintf(f, "          $write(%c relerr=%c,relerr);\n",'"','"');
  fprintf(f, "          $display(%c time= %cd%c,$time);\n",'"','%','"');
  fprintf(f, "        end\n");
  fprintf(f, "     end\n");
  fprintf(f, "  `endif\n");
  fprintf(f, "\n");

}



void print_mul_residueLNS()
{
  int i;
  FILE * f;



  f = fopen(("mul2rlns.v"),"w");
  fprintf(f,"// moduli [");
  for (i=0; i<MAXMODULI; i++)
    fprintf(f, "% d", moduli.m[i]);
  fprintf(f, "]\n");
  print_check_moduli_define(f,0);

  fprintf(f, "module mul_residueLNS(s,a,b);\n");
  fprintf(f, "  input [%d:0] a,b;\n", MAXPACKBITS);
  fprintf(f, "  output [%d:0] s;\n", MAXPACKBITS);
  fprintf(f, "  reg [%d:0] s;\n", MAXPACKBITS);
  fprintf(f, "\n");

  print_check_relerr(f, '*');


  for (i=0; i<MAXMODULI; i++)
    print_add_function(f,moduli.m[i]);
  fprintf(f, "\n");
  fprintf(f, "  always @(a or b)\n");
  fprintf(f, "    begin\n");
  fprintf(f, "      s[%d] = a[%d] ^ b[%d];\n",MAXPACKBITS, MAXPACKBITS,MAXPACKBITS);

  for (i=0; i<MAXMODULI; i++)
    fprintf(f, "      s[%d:%d] = addmod%d(a[%d:%d], b[%d:%d]);\n",
		 left_bitpos.m[i], right_bitpos.m[i],
		 moduli.m[i],
		 left_bitpos.m[i], right_bitpos.m[i],
		 left_bitpos.m[i], right_bitpos.m[i]);
  fprintf(f, "    end\n");
  fprintf(f, "endmodule\n");
  fprintf(f, "\n");

  fclose(f);
}


void print_div_residueLNS()
{
  int i;
  FILE * f;



  f = fopen(("div2rlns.v"),"w");
  fprintf(f,"// moduli [");
  for (i=0; i<MAXMODULI; i++)
    fprintf(f, "% d", moduli.m[i]);
  fprintf(f, "]\n");
  print_check_moduli_define(f,0);

  fprintf(f, "module div_residueLNS(s,a,b);\n");
  fprintf(f, "  input [%d:0] a,b;\n", MAXPACKBITS);
  fprintf(f, "  output [%d:0] s;\n", MAXPACKBITS);
  fprintf(f, "  reg [%d:0] s;\n", MAXPACKBITS);
  fprintf(f, "\n");

  print_check_relerr(f, '/');


  for (i=0; i<MAXMODULI; i++)
    print_sub_function(f,moduli.m[i]);
  fprintf(f, "\n");
  fprintf(f, "\n");
  fprintf(f, "  always @(a or b)\n");
  fprintf(f, "    begin\n");
  fprintf(f, "      s[%d] = a[%d] ^ b[%d];\n",MAXPACKBITS, MAXPACKBITS,MAXPACKBITS);

  for (i=0; i<MAXMODULI; i++)
    fprintf(f, "      s[%d:%d] = submod%d(a[%d:%d], b[%d:%d]);\n",
		 left_bitpos.m[i], right_bitpos.m[i],
		 moduli.m[i],
		 left_bitpos.m[i], right_bitpos.m[i],
		 left_bitpos.m[i], right_bitpos.m[i]);
  fprintf(f, "    end\n");
  fprintf(f, "endmodule\n");
  fprintf(f, "\n");

  fclose(f);
}




void print_add_residueLNS_lowpack()
{
  int i, z;

  FILE * f;

  int modu_bits,in_bits;
  in_bits = (int) ceil(log(2*MAXPACKLOW3)/log(2.0));
  modu_bits = (int) ceil(log(moduli.m[3])/log(2.0));

  f = fopen(("add2rlns.v"),"w");
  fprintf(f,"// moduli [");
  for (i=0; i<MAXMODULI; i++)
    fprintf(f, "% d", moduli.m[i]);
  fprintf(f, "]\n");
  print_check_moduli_define(f,0);

  fprintf(f, "module add_residueLNS(s,a,b);\n");
  fprintf(f, "  input [%d:0] a,b;\n", MAXPACKBITS);
  fprintf(f, "  output [%d:0] s;\n", MAXPACKBITS);
  fprintf(f, "  reg [%d:0] s;\n", MAXPACKBITS);
  fprintf(f, "  reg [%d:0] sbdb;\n", MAXPACKBITS);
  fprintf(f, "  wire [%d:0] znew,z,negz;\n",MAXPACKBITS);
  fprintf(f, "\n");


  print_check_relerr(f, '+');

  for (i=0; i<MAXMODULI; i++)
    print_add_function(f,moduli.m[i]);
  for (i=0; i<MAXMODULI; i++)
    print_sub_function(f,moduli.m[i]);

  print_mulinvsel_function(f ,myM,moduli.m[3]);
  print_baseextend_function(f);
  print_sbdb_posneg_function(f);

  fprintf(f, "  assign z[%d] = a[%d]^b[%d];\n",MAXPACKBITS,MAXPACKBITS,MAXPACKBITS);

  for (i=0; i<MAXMODULI; i++)
  {
    fprintf(f, "  assign z[%d:%d] = submod%d(a[%d:%d], b[%d:%d]);\n",
		 left_bitpos.m[i], right_bitpos.m[i],
		 moduli.m[i],
		 left_bitpos.m[i], right_bitpos.m[i],
		 left_bitpos.m[i], right_bitpos.m[i]);
  }






  fprintf(f, "  assign negz[%d] = z[%d];\n",MAXPACKBITS,MAXPACKBITS);  //fix db bug 10/04/04

  for (i=0; i<MAXMODULI-1; i++)
  {
    fprintf(f, "  assign negz[%d:%d] = submod%d(b[%d:%d], a[%d:%d]);\n",
		 left_bitpos.m[i], right_bitpos.m[i],
		 moduli.m[i],
		 left_bitpos.m[i], right_bitpos.m[i],
		 left_bitpos.m[i], right_bitpos.m[i]);
  }


  fprintf(f, "  assign negz[%d:%d] = 0;\n",
		 left_bitpos.m[MAXMODULI-1], right_bitpos.m[MAXMODULI-1]);
  fprintf(f, "  assign znew = (z[%d]) ? negz : z;\n",
		 left_bitpos.m[MAXMODULI-1]+1);




  fprintf(f, "  wire [%d:0] sbdb_pos;\n", MAXPACKBITS);
  fprintf(f, "  wire [%d:0] sbdb_neg;\n", MAXPACKBITS);
  fprintf(f, "  wire [%d:0] sbdb_pos_plus_b;\n", MAXPACKBITS);
  fprintf(f, "  wire [%d:0] sbdb_neg_plus_a;\n", MAXPACKBITS);
  fprintf(f, "  wire [%d:0] baseext = baseextend(z>>%d);\n",modu_bits-1,modu_bits);
  fprintf(f, "  wire [%d:0] diff = submod%d(z[%d:0],baseext);\n",modu_bits-1,moduli.m[3],modu_bits-1);
  fprintf(f, "  wire [1:0] select = mulinv%d_%dsel(diff);\n",myM,moduli.m[3]);
  fprintf(f, "  assign {sbdb_pos,sbdb_neg} = sbdb_posneg(znew,z[%d]);\n",
		 left_bitpos.m[MAXMODULI-1]+1);

//?????????????????????


  fprintf(f, "  assign sbdb_neg_plus_a[%d] = sbdb_neg[%d]^a[%d];\n",MAXPACKBITS,MAXPACKBITS,MAXPACKBITS);

  for (i=0; i<MAXMODULI; i++)
  {
    fprintf(f, "  assign sbdb_neg_plus_a[%d:%d] = addmod%d(sbdb_neg[%d:%d], a[%d:%d]);\n",
		 left_bitpos.m[i], right_bitpos.m[i],
		 moduli.m[i],
		 left_bitpos.m[i], right_bitpos.m[i],
		 left_bitpos.m[i], right_bitpos.m[i]);
  }

  fprintf(f, "  assign sbdb_pos_plus_b[%d] = sbdb_pos[%d]^b[%d];\n",MAXPACKBITS,MAXPACKBITS,MAXPACKBITS);

  for (i=0; i<MAXMODULI; i++)
  {
    fprintf(f, "  assign sbdb_pos_plus_b[%d:%d] = addmod%d(sbdb_pos[%d:%d], b[%d:%d]);\n",
		 left_bitpos.m[i], right_bitpos.m[i],
		 moduli.m[i],
		 left_bitpos.m[i], right_bitpos.m[i],
		 left_bitpos.m[i], right_bitpos.m[i]);
  }


  fprintf(f, "  always @(a or b or select or sbdb_pos_plus_b or sbdb_neg_plus_a or z)\n");
  fprintf(f, "    begin\n");
  fprintf(f, "      case(select)\n");
  fprintf(f, "        0: s = (z[%d]<<%d)^sbdb_pos_plus_b;//bugfix\n",MAXPACKBITS,MAXPACKBITS);
  fprintf(f, "        1: s = (z[%d]<<%d)^sbdb_neg_plus_a;\n",MAXPACKBITS,MAXPACKBITS);
  fprintf(f, "        2: s = a; //z;\n");
  fprintf(f, "        3: s = b; //0;\n");
  fprintf(f, "      endcase\n");
  fprintf(f, "    end\n");



  fprintf(f, "endmodule\n");
  fprintf(f, "\n");



  fclose(f);
}




/****/

void print_add_1hotter_residueLNS_lowpack()
{
  int i, onehotbits, bitpos, leftbitpos;
  FILE * f;

  f = fopen(("add1rlns.v"),"w");
  fprintf(f,"// onehot moduli [");
  for (i=0; i<MAXMODULI; i++)
    fprintf(f, "% d", moduli.m[i]);
  fprintf(f, "]\n");
  print_check_moduli_define_onehot(f,0);


  fprintf(f, "module add_1hot_residueLNS(s,a,b);\n");

  onehotbits = 0;
  for (i=0; i<MAXMODULI; i++)
    onehotbits += moduli.m[i];


  fprintf(f, "  input [%d:0] a,b;\n", onehotbits);
  fprintf(f, "  output [%d:0] s;\n", onehotbits);
  fprintf(f, "  reg [%d:0] s;\n", onehotbits);
  fprintf(f, "  reg [%d:0] sbdb;\n", onehotbits);
  fprintf(f, "  wire [%d:0] znew,z,negz;\n",onehotbits);
  fprintf(f, "\n");



  print_check_relerr_1hot(f, '+');


  for (i=0; i<MAXMODULI; i++)
    print_add_1hot_function(f,moduli.m[i]);

  for (i=0; i<MAXMODULI; i++)
    print_sub_1hot_function(f,moduli.m[i]);
  fprintf(f, "\n");

  print_mulinvsel_1hot_function(f ,myM,moduli.m[3]);
  print_baseextend_1hot_function(f);
  print_sbdb_posneg_function(f);
  print_1hot2pack_function(f);
  print_pack21hot_function(f);

  fprintf(f, "  assign z[%d] = a[%d]^b[%d];\n", onehotbits, onehotbits, onehotbits);

  bitpos = 0;
  for (i=MAXMODULI-1; i >=0; i--)
  {
    leftbitpos = bitpos + moduli.m[i] - 1;
    fprintf(f, "  assign z[%d:%d] = sub1hmod%d(a[%d:%d], b[%d:%d]);\n",
		       leftbitpos, bitpos,
		       moduli.m[i],
		       leftbitpos, bitpos,
		       leftbitpos, bitpos);

    bitpos += moduli.m[i];
  }



  fprintf(f, "  assign negz[%d] = z[%d];\n",onehotbits,onehotbits);



  bitpos = 0;
  for (i=MAXMODULI-1; i >=0; i--)
  {
    leftbitpos = bitpos + moduli.m[i] - 1;
    fprintf(f, "  assign negz[%d:%d] = sub1hmod%d(b[%d:%d], a[%d:%d]);\n",
		       leftbitpos, bitpos,
		       moduli.m[i],
		       leftbitpos, bitpos,
		       leftbitpos, bitpos);

    bitpos += moduli.m[i];
  }


  fprintf(f, "  wire oddz = 0");
  for (i=1; i<moduli.m[MAXMODULI-2]; i += 2)
    fprintf(f,"|z[%d]", i + moduli.m[MAXMODULI-1]);
  fprintf(f, ";\n");
  fprintf(f, "  assign znew = (oddz) ? negz : z;\n");

  fprintf(f, "  wire [%d:0] sbdb_pos_pack;\n", MAXPACKBITS);
  fprintf(f, "  wire [%d:0] sbdb_neg_pack;\n", MAXPACKBITS);
  fprintf(f, "  wire [%d:0] sbdb_pos;\n", onehotbits);
  fprintf(f, "  wire [%d:0] sbdb_neg;\n", onehotbits);
  fprintf(f, "  wire [%d:0] sbdb_pos_plus_b;\n", onehotbits);
  fprintf(f, "  wire [%d:0] sbdb_neg_plus_a;\n", onehotbits);

  fprintf(f, "  wire [%d:0] baseext = baseextend1h(z>>%d);\n",moduli.m[3]-1,moduli.m[3]);
  fprintf(f, "  wire [%d:0] diff = sub1hmod%d(z[%d:0],baseext);\n",moduli.m[3]-1,moduli.m[3],moduli.m[3]-1);
  fprintf(f, "  wire [1:0] select = mulinv1h%d_%dsel(diff);\n",myM,moduli.m[3]);
  fprintf(f, "  assign {sbdb_pos_pack,sbdb_neg_pack} = sbdb_posneg(rlns1h2pack(znew),oddz);\n");
  fprintf(f, "  assign sbdb_pos = rlnspack21h(sbdb_pos_pack);\n");
  fprintf(f, "  assign sbdb_neg = rlnspack21h(sbdb_neg_pack);\n");

//?????????????????????


  fprintf(f, "  assign sbdb_neg_plus_a[%d] = sbdb_neg[%d]^a[%d];\n",onehotbits,onehotbits,onehotbits);


  bitpos = 0;
  for (i=MAXMODULI-1; i >=0; i--)
  {
    leftbitpos = bitpos + moduli.m[i] - 1;
    fprintf(f, "  assign sbdb_neg_plus_a[%d:%d] = add1hmod%d(sbdb_neg[%d:%d], a[%d:%d]);\n",
		       leftbitpos, bitpos,
		       moduli.m[i],
		       leftbitpos, bitpos,
		       leftbitpos, bitpos);

    bitpos += moduli.m[i];
  }


  fprintf(f, "  assign sbdb_pos_plus_b[%d] = sbdb_pos[%d]^b[%d];\n",onehotbits,onehotbits,onehotbits);


  bitpos = 0;
  for (i=MAXMODULI-1; i >=0; i--)
  {
    leftbitpos = bitpos + moduli.m[i] - 1;
    fprintf(f, "  assign sbdb_pos_plus_b[%d:%d] = add1hmod%d(sbdb_pos[%d:%d], b[%d:%d]);\n",
		       leftbitpos, bitpos,
		       moduli.m[i],
		       leftbitpos, bitpos,
		       leftbitpos, bitpos);

    bitpos += moduli.m[i];
  }




  fprintf(f, "  always @(a or b or select or sbdb_pos_plus_b or sbdb_neg_plus_a or z)\n");
  fprintf(f, "    begin\n");
  fprintf(f, "      case(select)\n");
  fprintf(f, "        0: s = (z[%d]<<%d)^sbdb_pos_plus_b;//bugfix\n",onehotbits,onehotbits);
  fprintf(f, "        1: s = (z[%d]<<%d)^sbdb_neg_plus_a;\n",onehotbits,onehotbits);
  fprintf(f, "        2: s = a; //z;\n");
  fprintf(f, "        3: s = b; //0;\n");
  fprintf(f, "      endcase\n");
  fprintf(f, "    end\n");



  fprintf(f, "endmodule\n");
  fprintf(f, "\n");



  fclose(f);
}




unsigned long residue21hot(residuetype x)
{
    unsigned long s;
    int i;

    s = 0;
    for (i  = 0; i< MAXMODULI; i++)
    {
      s = s << moduli.m[i];
      s |= 1 << x.m[i];

      //printf(" %d %d %08lx\n",x.m[i], moduli.m[i], s);
      //getchar();
    }
    return(s);
}

void print_naive_sbdb_1hot()
{
  int i, z, zi, mask, onehotbits;
  unsigned long zsign, zsignonehot;

  FILE * f;

  f = fopen(("nsd1rlns.v"),"w");
  fprintf(f,"// moduli [");
  for (i=0; i<MAXMODULI; i++)
    fprintf(f, "% d", moduli.m[i]);
  fprintf(f, "] onehot\n");
  print_check_moduli_define_onehot(f,0);

  fprintf(f, "module sbdb_logic_1hot(sbdb, z);\n");
  onehotbits = 0;
  for (i=0; i<MAXMODULI; i++)
    onehotbits += moduli.m[i];

  fprintf(f, "  input [%d:0] z;\n", onehotbits);
  fprintf(f, "  output [%d:0] sbdb;\n", onehotbits);
  fprintf(f, "  reg [%d:0] sbdb;\n", onehotbits);
  fprintf(f, "\n");
  fprintf(f, "  always @(z)\n");
  fprintf(f, "    begin\n");
  fprintf(f, "      case(z)\n");
  fprintf(f, "\n");

  mask = MAXPACK-1;
  for (zsign=0; zsign<=1; zsign++)
  for (zi = MIN; zi <= MAX; zi++)
  {
    z = zsign*MAXPACK + res2pack(BIN2RES(zi));
    zsignonehot = zsign << onehotbits;
    fprintf(f, "         %d'h%lx : sbdb = {%d,%d'h%lx}; //",
	onehotbits+1,zsignonehot|residue21hot(BIN2RES(zi)),
	sbdb[z].xsign,onehotbits,residue21hot(sbdb[z].x));
    //fprint_residue(f, sbdb[z].x);
    fprintf(f, "\n");
    //if (z%16==0)
    //  getchar();
  }
  fprintf(f, "        default: sbdb = {0,%d'h%lx};\n",onehotbits,residue21hot(BIN2RES(0)));
  fprintf(f, "      endcase\n");
  fprintf(f, "    end\n");
  fprintf(f, "endmodule\n");
  fprintf(f, "\n");

  fclose(f);
}

void print_naive_sbdb_1hotalt()
{
  int i, z, zi, mask, onehotbits;
  unsigned long zsign, zsignonehot;
  int j,bitpos;
  residuetype zres;
  FILE * f;

  f = fopen(("alt1rlns.v"),"w");
  fprintf(f,"// moduli [");
  for (i=0; i<MAXMODULI; i++)
    fprintf(f, "% d", moduli.m[i]);
  fprintf(f, "] onehot\n");
  print_check_moduli_define_onehot(f,0);

  fprintf(f, "naive algorithm\n");
  fprintf(f, "module sbdb_logic_1hot(sbdb, z);\n");
  onehotbits = 0;
  for (i=0; i<MAXMODULI; i++)
    onehotbits += moduli.m[i];

  fprintf(f, "  input [%d:0] z;\n", onehotbits);
  fprintf(f, "  output [%d:0] sbdb;\n", onehotbits);
  fprintf(f, "  reg [%d:0] sbdb;\n", onehotbits);
  fprintf(f, "\n");
  fprintf(f, "  always @(z)\n");
  fprintf(f, "    begin\n");
  fprintf(f, "      sbdb = 0;\n");
  fprintf(f, "\n");

  mask = MAXPACK-1;
  for (zsign=0; zsign<=1; zsign++)
  for (zi = MIN; zi <= MAX; zi++)
  {
    z = zsign*MAXPACK + res2pack(BIN2RES(zi));
    zres = BIN2RES(zi);
    zsignonehot = zsign << onehotbits;

    bitpos = 0;
    for (i=MAXMODULI-1; i >=0; i--)
    {
      fprintf(f,"sbdb[%d] = sbdb[%d] | ",
	  sbdb[z].x.m[i]+bitpos,
	  sbdb[z].x.m[i]+bitpos);
      for (j=0; j<MAXMODULI; j++)
	fprintf(f, "z[%d]&", zres.m[j]+bitpos);  //needs offset
      if (zsign == 0) fprintf(f,"!");
      fprintf(f, "z[signbit];// mod %d\n",moduli.m[i]);
      fprintf(f, "//%d'h%lx : sbdb = {%d,%d'h%lx}; //",
	onehotbits+1,zsignonehot|residue21hot(BIN2RES(zi)),
	sbdb[z].xsign,onehotbits,residue21hot(sbdb[z].x));

      fprintf(f, "\n");
      bitpos += moduli.m[i];
    }
  }
  fprintf(f, "        default: sbdb = {0,%d'h%lx};\n",onehotbits,residue21hot(BIN2RES(0)));
  fprintf(f, "      endcase\n");
  fprintf(f, "    end\n");
  fprintf(f, "endmodule\n");
  fprintf(f, "\n");

  fclose(f);
}


void print_sbdb_binary()
{
  int i, z, zi, zsign, mask;

  FILE * f;

  f = fopen(("binLNS.v"),"w");
  fprintf(f,"// equiv to moduli [");
  for (i=0; i<MAXMODULI; i++)
    fprintf(f, "% d", moduli.m[i]);
  fprintf(f, "]\n");
  fprintf(f, "module sbdb_logic_bin(sbdb, z);\n");
  fprintf(f, "  input [%d:0] z;\n", MAXPACKBITS);
  fprintf(f, "  output [%d:0] sbdb;\n", MAXPACKBITS);
  fprintf(f, "  reg [%d:0] sbdb;\n", MAXPACKBITS);
  fprintf(f, "\n");
  fprintf(f, "  always @(z)\n");
  fprintf(f, "    begin\n");
  fprintf(f, "      case(z)\n");
  fprintf(f, "\n");

  mask = MAXPACK-1;
  for (zsign=0; zsign<=1; zsign++)
  for (zi = MIN; zi <= MAX; zi++)
  {
    z = zsign*MAXPACK + res2pack(BIN2RES(zi));
    fprintf(f, "          %3d : sbdb = {%d,%d'h%x}; //", zsign*MAXPACK+(mask&zi),sbdb[z].xsign,MAXPACKBITS,mask&RES2BIN(sbdb[z].x));
    fprintf(f, "\n");
  }
  fprintf(f, "        default: sbdb = 0;\n");
  fprintf(f, "      endcase\n");
  fprintf(f, "    end\n");
  fprintf(f, "endmodule\n");
  fprintf(f, "\n");

  fclose(f);
}


void print_sbdb_test()
{
  int i, z, zi, zsign, mask;
  int myMin, myMax, myM;

  FILE * f;

  f = fopen(("test1LNS.v"),"w");
  fprintf(f,"// equiv to moduli [");
  for (i=0; i<MAXMODULI; i++)
    fprintf(f, "% d", moduli.m[i]);
  fprintf(f, "]\n");
  fprintf(f, "module sbdb_logic_test(sbdb, z);\n");
  fprintf(f, "  input [%d:0] z;\n", MAXPACKBITS);
  fprintf(f, "  output [%d:0] sbdb;\n", MAXPACKBITS);
  fprintf(f, "  reg [%d:0] sbdb;\n", MAXPACKBITS);
  fprintf(f, "\n");
  fprintf(f, "  always @(z)\n");
  fprintf(f, "    begin\n");
  fprintf(f, "      case(z)\n");
  fprintf(f, "\n");

  myM = moduli.m[3]*moduli.m[1]*moduli.m[2];
  myMin = -myM/2;
  myMax = myM + myMin;
  printf("%d %d %d\n",myM, myMin, myMax);

  mask = MAXPACK-1;
  for (zsign=0; zsign<=1; zsign++)
  for (zi = myMin; zi <= myMax; zi++)
  {
    z = zsign*MAXPACK + res2pack(BIN2RES(zi));
    fprintf(f, "          %3d : sbdb = {%d,%d'h%x}; //", zsign*MAXPACK+(mask&zi),sbdb[z].xsign,MAXPACKBITS,mask&RES2BIN(sbdb[z].x));
    fprintf(f, "\n");
  }
  fprintf(f, "        default: sbdb = 0;\n");
  fprintf(f, "      endcase\n");
  fprintf(f, "    end\n");
  fprintf(f, "endmodule\n");
  fprintf(f, "\n");

  fclose(f);
}



void setup_directory()
{
  char onemoduli[100];
  char s[100];
  int i;
  int success;

  sprintf(s, RLNSPATH);
  strcat(s, "/");
  for (i=0; i<MAXMODULI; i++)
  {
    sprintf(onemoduli, "%02x", moduli.m[i]);
    strcat(s, onemoduli);
  }
  printf("%s\n",s);
  success = chdir(s);
  if (success != 0)
  {
    printf("having to mkdir\n");
    mkdir(s);
    success = chdir(s);
    if (success != 0)
    {
      printf("mkdir failure\n");
      exit(1);
    }
  }
  else
  {
    printf("reusing existing directory\n");
  }
  //getchar();

}

void abnormal_exit()
{
  FILE * f;
  int success;

  success = chdir(RLNSPATH);
  //f=fopen("test.txt","w");
  //fprintf(f,"hi %d\n",success);
  //fclose(f);
  exit(1);
}


void main()
{


  init_residue_constants();
  setup_directory();

  print_print_residueLNS();
  print_print_residueLNS_1hot();
  print_convert_residueLNS();
  print_convert_residueLNS_1hot();

  init_residue_tables();
  init_residue_lookup();
  init_sbdb_lowbits();

  print_add_residueLNS_lowpack();
  print_mul_residueLNS();
  print_div_residueLNS();
  print_mul_1hotter_residueLNS();
  print_div_1hotter_residueLNS();
  print_add_1hotter_residueLNS_lowpack();

  print_naive_add_residueLNS();
  print_naive_sbdb();
  print_naive_add_1hot_residueLNS();
  print_naive_sbdb_1hot();
  print_naive_sbdb_1hotalt();

  abnormal_exit();

/*
  init_residue_bipart();
  test_residue_bipart1();
  print_bipart_residue();
*/
}


/* mn-00a.c */

/* Soluci'on num'erica de 
 *
 * y' = -2ty**2,  0<=t<=2,  y(0)=1
 *
 * usando el m'etodo de Euler y el de los 3 t'erminos de la serie 
 * de Taylor.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double f_1(double t,double y);
double  euler_1(double h,double y1,double tk);
double taylor_1(double h,double y2,double tk);

int main(int argc,char *argv[])
{
  char s_file[80];                         //En esta variable se indica cuantos caracteres puede tener el nombre del archivo, en este caso './dat/mn-00a.dat'
  double h,tk,t0,t1,y0,y1,y2,y3;
  FILE *f;

  // valores iniciales
  h = 0.01;                                //h es la delta de tiempo
  t0 = 0.;
  t1 = 2.;  
  y0 = y1 = y2 = y3 = 1.;
  //archivo
  sprintf(s_file,"./dat/mn-00a.dat");    //Para que esto funcione necesitamos que ya exista la carpeta dat, ya que no lo crea en automatico
  printf("%s\n",s_file);                 //El archivo mn-00a.dat si lo crea solito o lo sobreescribe
  f = fopen(s_file,"w");
  for (tk=t0;tk<=t1;tk+=h) {                          //tk es lo que va a valer el tiempo en cada iteracion, por eso se va sumando, en el ciclo for, se indica el inicio, la condicion final y el paso que va a dar, en este caso inicia en tk=t0, acaba en tk > a t1, y en cada iteracion va aumentando h a tk
    fprintf(f,"%f\t%f\t%f\t%f\n",tk,y0,y1,y2);        //Aca se imprimen los resultados en el formato que indica, %f es un numero flotante y \t es un tabulador para separar
    y0 = 1./(1.+tk*tk); //soluci'on anal'itica
    y1 = euler_1(h,y1,tk);
    y2 = taylor_1(h,y2,tk);
    //    y3 += frk_1(h,tk,y3);
  }
  fclose(f);
  return 0;
}
                                           //A partir de aca se ponen las funciones que se utilizan arriba en el for
double f_1(double t,double y)
{
  return -2.*t*pow(y,2.);
}

double  euler_1(double h,double y1,double tk)    //Es un aproximacion de primer orden, el error es de delta t, este metodo consiste en encontrar la integral a partir de la suma de las areas bajo la curva, con un delta t que tiende a cero
{
  return y1+f_1(tk,y1)*h;
}

double taylor_1(double h,double y2,double tk)     //Esta aproximacion es de segundo orden, entonces el error es de dt**2
{
  double y;

  y  = y2+h*f_1(tk,y2)+h*h*y2*y2*(4.*tk*tk-1.);   //Esta aproximacion es de cuarto orden, error = dt**4
  return y;
}

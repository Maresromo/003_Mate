/* mn-1.c */
/* 28-08-12 */
/* 04-09-12 */
/* 23-09-12 */

/* M'etodos num'ericos para ecuaciones y' = f(t,y)
 *
 * y0 soluci'on anal'itica
 * y1 M'etodo de Euler
 * y2 M'etodo de los 3 t'erminos de la serie de Taylor
 * y3 M'etodo de Runge-Kutta de 4o orden
 *
 * t0 = tiempo inicial
 * tmax = tiempo final
 * ya = y(t0)
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mn1a.c"

int main()
{
  char s_file[80];
  int prob; 
  double area,h,H,k,sigma,tk,t0,t1,tmax,y0,y1,y2,y3,ya;
  FILE *f_out;

  printf(" 1. y' = -2*t*y**2, y(0) = 1\n");
  printf(" 2. y' = k*y*(1-y)\n");
  printf(" 3. y' = k*y*(1-y)-H\n");
  printf(" 4. y' = y(1 + exp(-y)) + exp(t), y(0) = 0\n");
  printf(" 5. y' = y*sin(2t), y(0) = 1, 0<=t<=2pi\n");
  printf(" 6. y' = y^2 cos(2t), y(0) = 1, 0<=t<=1/4\n");
  printf(" 7. y' = (y^0.5) cos(2t), y(0) = 1\n");
  printf(" 8. y' = (1/y)cos(2t), y(0) = 2\n");
  printf(" 9. y' = y-1, y(0) = 1\n");
  printf("10. y' = (3-y)(y+1), y(0)=4\n");
  printf("11. y' = (3-y)(y+1), y(0)=0\n");
  printf("12. y' = -y^2, y(0)=1/2\n");
  printf("13. y'=t/(y-t^2y), y(0)=4, 0<=t<0.8\n");
  printf("14. y'=(1-y^2)/y, y(0)=-2, 0<=t<=3\n");
  printf("15. y'=3y^2-4ty^2,y(0)=1, 0<=t<.45\n");
  printf("16. y'=t^2y-2-2y+t^2, y(0)=1, 0<t<2.5\n");                        
  printf("17. y'=1/(ty+t+y+1), y(0)=1, 0<t<5\n");
  printf("18. y'=1+t-y, t0=0, y(t0)=0, t1=1 [y=t]\n");
  printf("19. y'=1+2t+y^2/(1+t^2)^2, t0=0, y(t0)=1, t1=1, [y=1+t^2]\n");
  printf("20. y'=2ty, t0=0, y(t0)=2, t1=1, [y=2e^{t^2}]\n");
  printf("21. w'=(3-w)*(w+1), t0=0 y(t0)=1, t1=5\n");
  printf("22. y'=y^2-y^3, t0=0, y(t0)=0.2, t1=10\n");                       
  printf("23. y'=2*y^3+t^2,  y(t0)=0, t0=0, t_max =1.5\n");
  printf("24. y'=cos(y)\n");
  printf("25. y'=(y/t^2)+4*cos(t), t0=1, y(t0)=0, t1=20\n");
  printf("26. y'=(y^2-4)*(sen^2y^3+cos y-2)/2, y(0)=1/2\n");
  printf("27. y'=(1+t^2)/(1+y^2)\n");
  printf("28. y'=(1/sqrt(M_PI))*exp(x*x/2/sigma), [-sigma:sigma]\n");
  printf("29. y'=y^2-4ty+yt^2-4y+8t-3\n");
  printf("Problema = ");
  scanf("%d",&prob);
  switch (prob) {
  case 1:
    ya = 1.;  
    h = 0.01;
    y3 = y2 = y1 = y0 = ya; 
    tmax = 2.;
    sprintf(s_file,"./dat/mn1-%d.dat",prob);
      printf("%s\n",s_file);
    f_out = fopen(s_file,"w");
    fprintf(f_out,"%f\t%f\t%f\t%f\t%f\n",tk,y0,y1,y2,y3);
    for (tk=h;tk<=tmax;tk+=h) {
      y0 = 1./(pow(tk,2.)+1);
      y1 = euler_1(h,y1,tk);
      y2 = taylor_1(h,y2,tk);
      y3 = rk_1(h,y3,tk);
      fprintf(f_out,"%f\t%f\t%f\t%f\t%f\n",tk,y0,y1,y2,y3);
      //      printf("%f\t%f\t%f\t%f\t%f\n",tk,y0,y1,y2,y3);
    }
    fclose(f_out);
    break;
  case 2:
    printf("k = ");
    scanf("%lf",&k);
    printf("0<y0<1, y0 = ");
    scanf("%lf",&ya);
    h = 0.01;
    y3 = y2 = y1 = y0 = ya;
    t0 = 0.;
    tmax = 6; 
    sprintf(s_file,"./dat/mn1-%d-y%d-t%d.dat",prob,(int) (ya*1.e5),(int) (t0)); 
    printf("%s\n",s_file);
    f_out = fopen(s_file,"w");
    for (tk=t0;tk<=tmax;tk+=h) {
      fprintf(f_out,"%f\t%f\t%f\t%f\n",tk,y1,y2,y3);
      y1 = euler_2(h,k,y1);
      y2 = taylor_2(h,k,y2);
      y3 = rk_2(h,k,y3); 
    }
    fclose(f_out);
   break;
  }
  return 0;
}


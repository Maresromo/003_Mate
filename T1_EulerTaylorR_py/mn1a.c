/* funciones de mmn-01.c */

double f_1(double t,double y)
{
  return -2.*t*pow(y,2.);
}

double  euler_1(double h,double y1,double tk)
{
  return y1+f_1(tk,y1)*h;
}

double taylor_1(double h,double y2,double tk)
{
  double y;

  y  = y2+h*f_1(tk,y2)+h*h*y2*y2*(4.*tk*tk*y-1.);
  return y;
}

double rk_1(double h,double y3,double tk) 
{
  double rk1,rk2,rk3,rk4;

  rk1 = f_1(tk,y3);
  rk2 = f_1(tk+h/2,y3+rk1*h/2.);
  rk3 = f_1(tk+h/2,y3+rk2*h/2.);
  rk4 = f_1(tk+h,y3+h*rk3);
  return y3+h*(rk1+2.*rk2+2.*rk3+rk4)/6.; 
}

///////////////////////////////////////////////////////////////////

double f_2(double k,double y)
{
  return k*y*(1-y);
}

double euler_2(double h,double k,double y1)
{
  return y1+f_2(k,y1)*h;
}

double taylor_2(double h,double k,double y2)
{
  double y;

  y  = y2+h*f_2(k,y2)+(h*h/2.)*f_2(k,y2)*k*(1.-2.*y2);
  return y;
}

double rk_2(double h,double k,double y3)
{
  double rk1,rk2,rk3,rk4;

  rk1 = f_2(k,y3);
  rk2 = f_2(k,y3+rk1*h/2.);
  rk3 = f_2(k,y3+rk2*h/2.);
  rk4 = f_2(k,y3+h*rk3);
  return y3+h*(rk1+2.*rk2+2.*rk3+rk4)/6.; 

}

///////////////////////////////////////////////////////////////////


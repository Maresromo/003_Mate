#include <stdio.h>
#include <math.h>

int main()
{

int N, t;
float xa, z, g, jo;
float det, temp, fac, r, va;
int n, pv, i, j, jc, jr, k, kc, nv;

        printf("\n\t\t Interpolaciion por Splines cubicos\n\n" );
        printf("Escribe el numero de datos que posees \n");
        scanf("%d", &N);

float x[N-1];
float y[N-1];

for(t=0; t<N; t++)
{
        printf("Escribe el valor de x %d \n", (t+1));
        scanf("%f", &x[t]);
        printf("Escribe el valor de y %d \n", (t+1));
        scanf("%f", &y[t]);
}// Para las n-2 ecuaciones, con los coeficientes almacenados en columnas
        printf("\n\n");
    printf("La matriz resultante es\n");

    //    printf("k/\t\t\tl\t\t\tm\t\t\tres\n");

float a_int[N][N];
//static float a_int[N][N];
n=N-2;
for (i=1; i<=n; i++){
        for(j=1; j<=(n+1); j++){
		a_int[i][j]=0;
//               printf("%f\t",  a_int[i][j]);
        }

}

        i=1;
        for(i=1; i<(N-1); i++)
        {

		a_int[i][i-1]=(x[i]-x[i-1]);
                a_int[i][i]=2*(x[i+1]-x[i-1]);
                a_int[i][i+1]=(x[i+1]-x[i]);
                a_int[i][N-1]=((((6)/(x[i+1]-x[i]))*(y[i+1]-y[i]))+((((6)/(x[i]-x[i-1]))*(y[i-1]-y[i]))));

                a_int[1][0]=0;
               // a_int[N-2][N-1]=0;

        }


n=N-2;

for (i=1; i<=n; i++){
        for(j=1; j<=(n+1); j++){

               printf("%f\t",  a_int[i][j]);
        }
        printf("\n\n");
}

double a[N][N];

n=N-2;
det=1;

printf("Eliminacion Gaussiana\n\n");

for (i=1; i<=n; i++){
        for(j=1; j<=(n+1); j++){
                a[i][j]=a_int[i][j];
               printf("%f\t", a[i][j]);
        }
        printf("\n\n");
}
printf("las soluciones al sistema de ecuaciones son\n");
for(i=1; i<=(n-1); i++){        /*Inicio de FOR_1*/
        pv=i;

        for(j=i+1; j<=n; j++){  /*Este for determina si es necesario un pivoteo*/
                if(fabs(a[pv][i])<fabs(a[j][i]))
                {
                        pv=j;
                 }
        }

        if(pv!=i){ /*este if pivotea la matiz*/
                for(jc=1; jc<=(n+1); jc++){
                        temp=a[i][jc];
                        a[i][jc]=a[pv][jc];
                        a[pv][jc]=temp;

                        det=(-1)*det;
                }
        } /*Aqui ya esta pivoteada la matriz*/

        if(a[i][i]==0){
                printf("La matriz es singular, no tiene solucion\n");
                break;
        }

        for(jr=i+1; jr<=n; jr++){ /*Comienza la eliminacion de los elementos debajode la diagonal */
                if(a[jr][i]!=0){
                        r=a[jr][i]/a[i][i];
                        for(kc=i+1; kc<=n+1; kc++){
                                temp=a[jr][kc];
                                a[jr][kc]=a[jr][kc]-(r*a[i][kc]);
                        }
                }
        }
}               /*Fin de For_1*/

for(i=1; i<=n; i++){
        det=det*a[i][i];        /*Se calcula el determinante de la matriz*/
}

if(det==0){
        printf("Matriz singular\n");    /*Se indica si la matriz es singular*/
}

a[n][n+1]=a[n][n+1]/a[n][n];

for(nv=(n-1); nv>=1; nv--){
        va=a[nv][n+1];
        for(k=nv+1; k<=n; k++){
                va=va-(a[nv][k]*a[k][n+1]);
        }
        a[nv][n+1]=va/a[nv][nv];
}

for(i=1; i<=n; i++)
{
       // printf("\t%d\t%16.5e\n",i, a[i][n+1]);
        printf("\t%d\t%f\n",i, a[i][n+1]);
// printf("El determinante es: %f\n", det);
printf("-------------------------------------\n");
}

float xz, V, W, res;
int B,L;

printf("Cual es el valor que deseas interpolar?\n");
printf("Tiene que estar entre\t%f", x[0]);
printf("\t y \t %f", x[N-1]);
printf("\n");
scanf("%f", &xz);
printf("\n");

B=0;

if (xz<x[0])
{
        printf("X tiene que estar entre\t%f", x[0]);
}
if (xz>x[N-1])
{
        printf("X tiene que estar entre\t%f", x[0]);
}
else
{
        B=0;
        while (xz>x[B])
        {
        B=B+1;
        }
L=B-1;
//printf("L es\t%i", L);
printf("\n");

if (L<1)
{
	V=0;
     	W=a[1][n+1];
}
if (L>n-1)
{
	W=0;
     	V=a[n][n+1];
}
else 
{
	V=a[L][n+1];
	W=a[L+1][n+1];
}
//printf("V\t%f",V);
//printf("W\t%f",W);
//printf("Xi-1\t%f",x[L]);
//printf("xi\t%f",x[L+1]);
//printf("yi-1\t%f",y[L]);
//printf("yi\t%f",y[L+1]);


res=(((V*pow((x[L+1]-xz),3)/(6*(x[L+1]-x[L])))+(W*pow((xz-x[L+1]),3)/(6*(x[L+1]-x[L])))+(((y[L]/(x[L+1]-x[L]))-((V*(x[L+1]-x[L])/6)))*(x[L+1]-xz))+((((y[L+1]/(x[L+1]-x[L]))-((W*(x[L+1]-x[L])/6)))*(xz-x[L])))));
printf("El resultado  es\t%f\n\t", res);
printf("\n");
}
}

#include "gauss_col.h"
#include "dopri8.h"

double norm(double *q,int q_len);
void input_start_params(double *x,double *p,double *t0,double *T,double *h);
void input_x0(double *x,double *p,double *t0,double *T,double *h);
void vfn(double *x,double *q,double *p);
void f(double *x,double *p,double *q,double *t0,double *T,double *h);
int method_Newton(double *A,double *x,double *p,double *q,int n,int p_len,int q_len,double *t0,double *T,double *h);
int f_end(double *q,int q_len);
int smesh(double *A,double *x,double *p,double *q,double *dp,int n,double *t0,double *T,double *h);
int gradf(double *x,double *p,double *q,double *A,int n,double *t0,double *T,double *h);

int main()
{
    //Инициализация всего необходимого
    int n=9,p_len=9,q_len=9,success,i;
    double *p,*q,*x,*A,t0[1],T[1],h[1];
    FILE *fw;
    
    if(!(fw=fopen(OUTFILE,"w")))
        printf("Невозможно открыть файл \'%s\' для вывода!\n",OUTFILE);
    else
    {
        fprintf(fw,"Начальные значения параметров пристрелки:\n");
        fprintf(fw,"    T=%.0lf секунд; Px=%.0lf; Py=%.0lf\n\n",Tstart,Px_start,Py_start);
//        fprintf(fw,"           t     %sr %sфи%su %sv %sm %sPr   %sPu  %sPv\n",space,space,space,space,space,space,space,space);
        fclose(fw);
    }
    printf("\n\n");
    
    if(!(x=new double[n]))
    {
        printf("Недостаточно памяти!\n");
        return -1;
    }
    if(!(p=new double[p_len*4]))
    {
        printf("Недостаточно памяти!\n");
        delete[] x;
        x=0;
        return -1;
    }
    if(!(q=new double[q_len*3]))
    {
        printf("Недостаточно памяти!\n");
        delete[] x;
        x=0;
        delete[] p;
        p=0;
        return -1;
    }
    if(!(A=new double[2*p_len*p_len+p_len]))
    {
        printf("Недостаточно памяти!\n");
        delete[] x;
        x=0;
        delete[] p;
        p=0;
        delete[] q;
        q=0;
        return -1;
    }
    printf("Летим с круговой орбиты Земли %.2lf км. на круговую орбиту Земли %.2lf км.\n\n",h0,hT);
    //Ввод начальных приближений
    input_start_params(x,p,t0,T,h);
    //Метод Ньютона
    success=method_Newton(A,x,p,q,n,p_len,q_len,t0,T,h);
    //Еще один пробег
    if(success==8)
        printf("  УСПЕХ \n\n");
    else
        printf("  НЕУДАЧА \n\n");
    
    input_x0(x,p,t0,T,h);
    dopri8_print(OUTFILE,*t0,x,*T,dopriACC,step_max,*h);
    vfn(x,q,p);

    printf("Время перелета: %.2lf секунд\n",*T);
    printf("Остаток массы:  %.2lf%%\n",x[4]*100);
    printf("Невязки:\n  ");
    printf("r[T]-R_Moon-hT = %e км;  ",sqrt(x[0]*x[0]+x[1]*x[1])-R_Moon-hT);
    printf("по u[T] = %e км/с;  ",q[1]);
    printf("по v[T] = %e км/с;  ",q[2]);
    printf("по сопр. перем. = %e\n",q[2]);
    printf("Полученная траектория лежит в файле \'%s\'\n",OUTFILE);
    printf("\n\n");
    printf("Teta(0)=%lf\n",p[1]);
    
    for(i=0;i<9;i++)
	printf("q[%d]=%e\n",i,q[i]);
    if(fw=fopen(OUTFILE,"a"))
    {
//       fprintf(fw,"\nr[T]-R_Earth-hT = %e;  u[T] = %e;  v[T]-v_cir_hT = %e\n",\
//                x[0]-R_Earth-hT,x[2],x[3]-v_circ_hT);
//        fprintf(fw,"H = %.18lf\n",x[5]*x[2]+sqrt(x[6]*x[6]+x[7]*x[7])*Pmax/x[4]+(x[6]/x[0])*(x[3]*x[3]-mu_Earth/x[0])-\
//                                  x[7]*x[2]*x[3]/x[0]);
        fclose(fw);
    }
    delete[] x;
    x=0;
    delete[] p;
    p=0;
    delete[] q;
    q=0;
    delete[] A;
    A=0;
    
    return 0;
}

void input_start_params(double *x,double *p,double *t0,double *T,double *h)
{
    p[0]=0;
    p[1]=1;
    p[2]=0.00001;
    p[3]=0.00001;
    p[4]=1;
    p[5]=1;
    p[6]=1;
    p[7]=0.5;
    p[8]=0.5;
    input_x0(x,p,t0,T,h);
    return;
}

//       x[0] x[1] x[2] x[3] x[4] x[5] x[6] x[7] x[8]
//   t    x    y    u    v    m    Px   Py   Pu   Pv  Chi  Teta

void input_x0(double *x,double *p,double *t0,double *T,double *h)
{
    x[0]=R_Earth+h0;
    x[1]=0.;
    x[6]=p[0];
    x[7]=p[1];
    x[2]=p[0];
    x[3]=p[1]+v_circ_h0;
    x[4]=p[2];
    x[5]=p[3];
    *T=86400*2.5;
    *t0=0.;
    *h=step;
    
    return;
}

void f(double *x,double *p,double *q,double *t0,double *T,double *h)
{
    input_x0(x,p,t0,T,h);

    dopri8(*t0,x,*T,dopriACC,step_max,*h);
	
    vfn(x,q,p);
    
    return;
}

void vfn(double *x,double *q,double *p)
{
    q[0]=x[0]*x[0]+x[1]*x[1]-(R_Moon+hT)*(R_Moon+hT);
    q[1]=x[2]+p[7]+v_circ_hT*x[1]/(R_Moon+hT);
    q[2]=x[3]+p[8]-v_circ_hT*x[0]/(R_Moon+hT);    
    q[3]=-x[5]-2*p[4]*x[0]-p[5]*(x[2]+p[7]);
    q[4]=-x[6]-2*p[4]*x[1]-p[5]*(x[3]+p[8]);
    q[5]=-x[7]-2*p[6]*(x[2]+p[7])-p[5]*x[0];
    q[6]=-x[8]-2*p[6]*(x[3]+p[8])-p[5]*x[1];
    q[7]=2*p[6]*(x[2]+p[7])+p[5]*x[0]+p[7];
    q[8]=2*p[6]*(x[3]+p[8])+p[5]*x[1]+p[8];

    return;
}

int method_Newton(double *A,double *x,double *p,double *q,int n,int p_len,int q_len,double *t0,double *T,double *h)
{
    int i,j,it=0,allpoints,flag;
    double gamma,qnorm,qnorm1,*dp=p+(p_len<<1),*p1=p+p_len*3,*q1=q+q_len*2;
    
    f(x,p,q,t0,T,h);
    if((flag=f_end(q,q_len))==0)
    {
	printf("Число интераций в методе Ньютона: %d\n",it);
	return 8;
    }
	
    for(it=1,allpoints=1,gamma=1.;it<MAXIT;it++)
    {
	gamma*=2.;
	gamma=(gamma<1.?gamma:1.);
	
	if((flag=smesh(A,x,p,q,dp,p_len,t0,T,h))<0)
        {
            printf("\n\nError smesh!\n\n\nCan't calculate derivative on %d iteration in modified Newton's method\n",it);
	    printf("Всего рассмотренных точек: %d\n",allpoints);
	    return flag;
	}
	
	qnorm=norm(q,q_len);
	
	for(j=0,allpoints++;j<30;j++,allpoints++)
	{
	    for(i=0;i<p_len;i++)
		p1[i]=p[i]+gamma*dp[i];
	    f(x,p1,q1,t0,T,h);
	    qnorm1=norm(q1,q_len);
	    if(qnorm<qnorm1)
    	        gamma*=0.5;
	    else
	    {
	        for(i=0;i<n;i++)
	        {
	            p[i]=p1[i];
		    q[i]=q1[i];
		}
		break; 
	    }
	}
printf("Метод Ньютона %d-я итерация, ошибка: %e\n",it,qnorm1);
	if(j==30)
        {
            printf("\n\nError!\n\n\n-7 Число интераций в методе Ньютона: %d, всего рассмотренных точек: %d\n",it,allpoints);
	    return -7;
	}
	if((flag=f_end(q,q_len))==0)
        {
            printf("Число интераций в методе Ньютона: %d, всего рассмотренных точек: %d\n",it,allpoints);
	    return 8;
	}
    }

    printf("\n\nError!\n\n\nreached maxit=%d limit in modified Newton's method\n",MAXIT);
    printf("Всего рассмотренных точек: %d\n",allpoints);
    if(it==MAXIT)
	return -1;
	
    return 0;
}

double norm(double *q,int q_len)
{
    double *a,x,s=0.;
    int i;
    
    for(i=0,a=q;i<q_len;i++)
    {
	x=*(a++);
	s+=x*x;
    }
    
    return sqrt(s);
}

int f_end(double *q,int n)
{
    int i;
    double *a=q;
    
    for(i=0;i<n;i++,a++)
    {
	if(*a>ACC||*a<-ACC)
	    return -1;
    }
    
    return 0;
}

int smesh(double *A,double *x,double *p,double *q,double *dp,int n,double *t0,double *T,double *h)
{
    int flag;
    
    if((flag=gradf(x,p,q,A,n,t0,T,h))<0)
        return flag;
    
    if((flag=gauss(A,dp,q,n)))
        return flag;
        
    for(flag=0;flag<n;flag++)
        dp[flag]=-dp[flag];
        
    return 0;
}

int gradf(double *x,double *p,double *q,double *A,int n,double *t0,double *T,double *h)
{
    double *p1=p+n,*q1=q+n+1,h_dp;
    int i,j,flag;
    for(j=0;j<n;j++)
    {
        for(i=0;i<n;i++)
            p1[i]=p[i];

        if(p1[j]>-1. && p1[j]<1.)
            h_dp=h_dp_eps;
        else
            h_dp=p1[j]*h_dp_eps*(p1[j]>0.?1.:-1.);
        
        p1[j]+=h_dp;
        
        f(x,p1,q1,t0,T,h);

        for(i=0;i<n;i++)
            A[j+i*n]=(q1[i]-q[i])/h_dp;
    }
    
    return 0;
}

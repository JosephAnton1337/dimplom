#include "gauss_col.h"
// n - размер всей матрицы; a - матрица n*n(заполненная); 
int gauss_col(double *a,double *b,int n)
{
    int i,j,q,nmax,nmax_n,i_n,j_n,n_i;
    double *c,*e,*f,*g,x,tmp,tmp8;
    
/////////////////////////////////////////////////////////////////////////////////////////

    //Прямой ход для первого столбца - записываем правую часть в первый столбец матрица A
    for(j=0,tmp=0.,nmax=0,j_n=0;j<n;j++,j_n+=n)
    {	//Ищем номер элемента, максимального по модулю в первом столбце
	tmp8=a[j_n];
	tmp8=(tmp8>0?tmp8:-tmp);
	if(tmp>tmp8-gauss_EPS)
	    continue;
	else
	{
	    tmp=tmp8;
	    nmax=j;
	}
    }
    
    if((tmp=a[nmax_n=nmax*n])<-gauss_EPS || tmp>gauss_EPS)
	tmp=1./tmp;
    else
	return 1;
	
    if(nmax!=0)
    {	//Переставляем местами строку i и nmax-овую
	tmp8=*a;
        a[nmax_n]=b[0]-(a[0]=b[nmax]*tmp)*tmp8;
	for(j=1,c=a+nmax_n;j<n;j++)
	{
	    x=a[j];
	    c[j]=x-(a[j]=c[j]*tmp)*tmp8;
	}
	x=*a;
    }
    else
    {
	x=*a=*b*tmp;
	for(j=1;j<n;j++)
	    a[j]*=tmp;
    }
    for(j=1,e=a+n,g=b+1;j<nmax;j++,e+=n,g++)
    {
	*e=*g-x*(tmp8=*e);
	for(q=1;q<n;q++)
	    e[q]-=a[q]*tmp8;
    }
    for(j=nmax+1,e=a+nmax_n+n,g=b+nmax+1;j<n;j++,e+=n,g++)
    {	//первый шаг прямого хода для nmax+1,..,n-1 строки
	*e=*g-x*(tmp8=*e);
	for(q=1;q<n;q++)
	    e[q]-=a[q]*tmp8;
    }
/////////////////////////////////////////////////////////////////////////////////////////
    //Прямой ход
    
    for(i=1,c=a+n+1,n_i=n-1,f=a+n;i<n;i++,c+=n+1,n_i--,f+=n)
    {
	for(j=0,tmp=0.,nmax=0,j_n=0;j<n_i;j++,j_n+=n)
	{
	    tmp8=c[j_n];
	    tmp8=(tmp8>0?tmp8:-tmp);
    	    if(tmp>tmp8-gauss_EPS)
		continue;
    	    else
	    {
		tmp=tmp8;
		nmax=j;
	    }
	}
	if((tmp=c[nmax_n=nmax*n])<-gauss_EPS || tmp>gauss_EPS)
	    tmp=1./tmp;
	else
	    return 1;
	if(nmax!=0)
	{
	    tmp8=*c;
	    x=*f;
	    f[nmax_n]=x-(*f=f[nmax_n]*tmp)*tmp8;
	    for(j=1,e=c+nmax_n;j<n_i;j++)
	    {
	        x=c[j];
		e[j]=x-(c[j]=e[j]*tmp)*tmp8;
	    }
	    x=*f;
	}
	else
	{
	    x=*f*=tmp;
	    for(j=1;j<n_i;j++)
	        c[j]*=tmp;
	}
	for(j=1,e=c+n,g=f+n;j<nmax;j++,e+=n,g+=n)
	{
	    *g-=x*(tmp8=*e);
	    for(q=1;q<n_i;q++)
	        e[q]-=c[q]*tmp8;
	}
	nmax_n+=n;
	for(j=nmax+1,e=c+nmax_n,g=f+nmax_n;j<n_i;j++,e+=n,g+=n)
	{
    	    *g-=x*(tmp8=*e);
	    for(q=1;q<n_i;q++)
	        e[q]-=c[q]*tmp8;
	}
    }
    
/////////////////////////////////////////////////////////////////////////////////////////    
    //Обратный ход
    for(i_n=(i=n-1)*n;i>0;i--,i_n-=n)
    {
	for(x=*(a+i_n),n_i=i_n-n,j=i;j>0;j--,n_i-=n)
	{
	    g=a+n_i;
	    *g-=*(e=g+i)*x;
	}
    }
    
    return 0;        
}

int gauss(double *A,double *x,double *b,int n)
{
    double *C=A+n*n,*e,*g,*acc,s,norm,norm1;
    int i,j;
    acc=C+n*n;
    
    for(i=0;i<n*n;i++)
        C[i]=A[i];
    
    if((i=gauss_col(A,b,n)))
        return i;
        
    for(i=0,e=x;i<n;i++)
        *(e++)=A[i*n];
    
    for(i=0,e=C,g=acc,norm=0.;i<n;i++)
    {
        for(j=0,s=0.;j<n;j++)
            s+=*(e++)*x[j];
            
        *(g++)=s-b[i];
        norm+=(s-b[i]>0.?s-b[i]:b[i]-s);
    }
    

    for(i=0;i<n*n;i++)
        A[i]=C[i];    

    gauss_col(C,acc,n);
    
    for(i=0;i<n;i++)
        acc[i]=x[i]-C[i*n];
    
    for(i=0,e=A,norm1=0.;i<n;i++)
    {
        for(j=0,s=0.;j<n;j++)
            s+=*(e++)*acc[j];
            
        norm1+=(s-b[i]>0.?s-b[i]:b[i]-s);
    }
    
    if(norm1<norm)
    {
        for(i=0;i<n;i++)
            x[i]=acc[i];
    }
    
    return 0;
}

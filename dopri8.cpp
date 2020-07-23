#include "dopri8.h"
#define NMAX 23000
//максимальное число шагов dopri
#define max(a,b) ((a)>(b) ? (a) : (b))
#define min(a,b) ((a)<(b) ? (a) : (b))
#define N 8
//размер массива 'x' для dopri


int len(double x);
int output(const char *in,const char *out,int *lengthes,int n,char *doprisays);
void fcn(double t,double *y,double *f);
double sgn(double x);


//p09:   x[0] x[1] x[2] x[3] x[4] x[5] x[6] x[7] x[8] x[9]
//   t    x    y    u    v    m    Px   Py   Pu   Pv   Pm  Chi  Teta


void fcn(double t,double *x,double *dx)
{
    double r=sqrt(x[0]*x[0] + x[1]*x[1]);
    double r3 =r * r * r;
    double r5 =r3 * r * r;
    double rho = sqrt(x[7]*x[7]+x[6]*x[6]);
    dx[0] = x[2];
    dx[1] = x[3];
    dx[2] = - mu_Moon*x[0]/r3;
    dx[3] = - mu_Moon*x[1]/r3;
    dx[4] = x[6] * (mu_Moon/r3 - 3.*mu_Moon*x[0]*x[0]/r5) - x[7] * 3.*mu_Moon*x[0]*x[1]/r5;
    dx[5] = x[7] * (mu_Moon/r3 - 3.*mu_Moon*x[1]*x[1]/r5) - x[6] * 3.*mu_Moon*x[0]*x[1]/r5;
    dx[6] = -x[4];
    dx[7] = -x[5];
    return;
    t=t;
}

static double   c2=1./18.,
                c3=1./12.,
                c4=0.125,
                c5=0.3125,
                c6=0.375,
                c7=0.1475,
                c8=0.465,
                c9=5490023248./9719169821.,
                c10=0.65,
                c11=1201146811./1299019798.,
                a21=1./18.,
                a31=1./48.,
                a32=0.0625, 
                a41=0.03125,
                a43=0.09375,
                a51=0.3125,
                a53=-1.171875,
                a54=1.171875,
                a61=0.0375,
                a64=0.1875,
                a65=0.15,
                a71=29443841./614563906.,
                a74=77736538./692538347.,
                a75=-28693883./1125000000.,
                a76=23124283./1800000000.,
                a81=16016141./946692911.,
                a84=61564180./158732637.,
                a85=22789713./633445777.,
                a86=545815736./2771057229.,
                a87=-180193667./1043307555.,
                a91=39632708./573591083.,
                a94=-433636366./683701615.,
                a95=-421739975./2616292301.,
                a96=100302831./723423059.,
                a97=790204164./839813087.,
                a98=800635310./3783071287.,
                a101=246121993./1340847787.,
                a104=-37695042795./15268766246.,
                a105=-309121744./1061227803.,
                a106=-12992083./490766935.,
                a107=6005943493./2108947869.,
                a108=393006217./1396673457.,
                a109=123872331./1001029789.,
                a111=-1028468189./846180014.,
                a114=8478235783./508512852.,
                a115=1311729495./1432422823.,
                a116=-10304129995./1701304382.,
                a117=-48777925059./3047939560.,
                a118=15336726248./1032824649.,
                a119=-45442868181./3398467696.,
                a1110=3065993473./597172653.,
                a121=185892177./718116043.,
                a124=-3185094517./667107341.,
                a125=-477755414./1098053517.,
                a126=-703635378./230739211.,
                a127=5731566787./1027545527.,
                a128=5232866602./850066563.,
                a129=-4093664535./808688257.,
                a1210=3962137247./1805957418.,
                a1211=65686358./487910083.,
                a131=403863854./491063109.,
                a134=-5068492393./434740067.,
                a135=-411421997./543043805.,
                a136=652783627./914296604.,
                a137=11173962825./925320556.,
                a138=-13158990841./6184727034.,
                a139=3936647629./1978049680.,
                a1310=-160528059./685178525.,
                a1311=248638103./1413531060.,
                b1=14005451./335480064.,
                b6=-59238493./1068277825.,
                b7=181606767./758867731.,
                b8=561292985./797845732.,
                b9=-1041891430./1371343529.,
                b10=760417239./1151165299.,
                b11=118820643./751138087.,
                b12=-528747749./2220607170.,
                b13=0.25,
                bh1=13451932./455176623.,
                bh6=-808719846./976000145.,
                bh7=1757004468./5645159321.,
                bh8=656045339./265891186.,
                bh9=-3867574721./1518517206.,
                bh10=465885868./322736535.,
                bh11=53011238./667516719.,
                bh12=2./45.;


double sgn(double x)
{
    if(x>0.) return 1.;
    if(x<0.) return -1.;
    return 0.;
}  


void dopri8(double x,double *y,double xend,double eps,double hmax,double h)
{
    int i,nfcn=0,nstep=0,naccpt=0,nrejct=0,reject=0;
    double posneg,y11s,y12s,xph,err,denom,fac,hnew,tmp;
    static double k1[N],k2[N],k3[N],k4[N],k5[N],k6[N],k7[N],y1[N];
    posneg=sgn(xend-x);
    hmax=fabs(hmax);
    h=min(max(1.e-10,fabs(h)),hmax);
    h=fabs(h)*posneg;
    eps=max(eps,13.*EPS);
    for(;;)
    {
        if(nstep>NMAX)
        {
            printf("Dopri: Слишком много шагов!\n");
            return; 
        }
        if(!(x+0.03*h>x) && !(x+0.03*h<x))
        {
            printf("Dopri:  Слишком маленький шаг! t=%.18lf\n",x);
            printf("        Обращений к fcn: %d; Общее число шагов: %d; Принятых шагов: %d; Смен шага: %d\n",nfcn,nstep,naccpt,nrejct);
            return;
        }
        if((x-xend)*posneg+EPS>0.) 
            return;
        if((x+h-xend)*posneg>0.)
            h=xend-x;
        fcn(x,y,k1);
        for(;;)
        {
            nstep++;
            for(i=0;i<N;i++)
                y1[i]=y[i]+h*a21*k1[i];
            fcn(x+c2*h,y1,k2);
            for(i=0;i<N;i++)
                y1[i]=y[i]+h*(a31*k1[i]+a32*k2[i]);
            fcn(x+c3*h,y1,k3);
            for(i=0;i<N;i++)
                y1[i]=y[i]+h*(a41*k1[i]+a43*k3[i]);
            fcn(x+c4*h,y1,k4);
            for(i=0;i<N;i++)
                y1[i]=y[i]+h*(a51*k1[i]+a53*k3[i]+a54*k4[i]);
            fcn(x+c5*h,y1,k5);
            for(i=0;i<N;i++)
                y1[i]=y[i]+h*(a61*k1[i]+a64*k4[i]+a65*k5[i]);
            fcn(x+c6*h,y1,k6);
            for(i=0;i<N;i++)
                y1[i]=y[i]+h*(a71*k1[i]+a74*k4[i]+a75*k5[i]+a76*k6[i]);
            fcn(x+c7*h,y1,k7);
            for(i=0;i<N;i++)
                y1[i]=y[i]+h*(a81*k1[i]+a84*k4[i]+a85*k5[i]+a86*k6[i]+a87*k7[i]);
            fcn(x+c8*h,y1,k2);
            for(i=0;i<N;i++)
                y1[i]=y[i]+h*(a91*k1[i]+a94*k4[i]+a95*k5[i]+a96*k6[i]+a97*k7[i]+a98*k2[i]);
            fcn(x+c9*h,y1,k3);
            for(i=0;i<N;i++)
                y1[i]=y[i]+h*(a101*k1[i]+a104*k4[i]+a105*k5[i]+a106*k6[i]+a107*k7[i]+a108*k2[i]+a109*k3[i]);
            for(i=0;i<N;i++) 
            {
                y11s=a111*k1[i]+a114*k4[i]+a115*k5[i]+a116*k6[i]+a117*k7[i]+a118*k2[i]+a119*k3[i];
                y12s=a121*k1[i]+a124*k4[i]+a125*k5[i]+a126*k6[i]+a127*k7[i]+a128*k2[i]+a129*k3[i];
                k4[i]=a131*k1[i]+a134*k4[i]+a135*k5[i]+a136*k6[i]+a137*k7[i]+a138*k2[i]+a139*k3[i];
                k5[i]=b1*k1[i]+b6*k6[i]+b7*k7[i]+b8*k2[i]+b9*k3[i];
                k6[i]=bh1*k1[i]+bh6*k6[i]+bh7*k7[i]+bh8*k2[i]+bh9*k3[i];
                k2[i]=y11s;
                k3[i]=y12s;
            }
            fcn(x+c10*h,y1,k7);
            for(i=0;i<N;i++)
                y1[i]=y[i]+h*(k2[i]+a1110*k7[i]);
            fcn(x+c11*h,y1,k2);
            xph=x+h;
            for(i=0;i<N;i++)
                y1[i]=y[i]+h*(k3[i]+a1210*k7[i]+a1211*k2[i]);
            fcn(xph,y1,k3);
            for(i=0;i<N;i++)
                y1[i]=y[i]+h*(k4[i]+a1310*k7[i]+a1311*k2[i]);
            fcn(xph,y1,k4);
            nfcn+=13;
            for(i=0;i<N;i++) 
            {
                k5[i]=y[i]+h*(k5[i]+b10*k7[i]+b11*k2[i]+b12*k3[i]+b13*k4[i]);
                k6[i]=y[i]+h*(k6[i]+bh10*k7[i]+bh11*k2[i]+bh12*k3[i]);
            }
            err=0.;
            for(i=0;i<N;i++) 
            {
                denom=max(max(1.e-6,fabs(k5[i])),max(fabs(y[i]),2.*EPS/eps));
                err+=((k5[i]-k6[i])/denom)*((k5[i]-k6[i])/denom);
            }
            err=sqrt(err/(double)(N));
            fac=max((1./6.),min(3.,pow(err/eps,(1./8.))/0.9));
            hnew=h/fac;
            if(err>eps)
            {
                reject=1;
                h=hnew;
                if(naccpt>1)
                    nrejct++;
                nfcn--;
                continue;
            }
            naccpt++;
            for(i=0;i<N;i++)
                y[i]=k5[i];
            x=xph;
            if(fabs(hnew)>hmax)
                hnew=posneg*hmax;
            if(reject==1)
                hnew=posneg*min(fabs(hnew),fabs(h));
            reject=0;
            h=hnew;
            break;
        }
    }
    
    return;
}


void dopri8_print(const char *filename,double x,double *y,double xend,double eps,double hmax,double h)
{
    int i,nfcn=0,nstep=0,naccpt=0,nrejct=0,reject=0,lengthes[N+1];
    double posneg,y11s,y12s,xph,err,denom,fac,hnew,tmp;
    static double k1[N],k2[N],k3[N],k4[N],k5[N],k6[N],k7[N],y1[N];
    char doprisays[180];
    posneg=sgn(xend-x);
    hmax=fabs(hmax);
    h=min(max(1.e-10,fabs(h)),hmax);
    h=fabs(h)*posneg;
    eps=max(eps,13.*EPS);
    FILE *fw=0;

    if(!(fw=fopen(DOPRITMP,"w")))
    {
        printf("Dopri: \"Невозможно вывести траектрорию в файл!!!\"\n");
        return;
    }
    fprintf(fw,"В этом файле лежит неформатированный вывод из Dopri\n\n");
    fprintf(fw,"%.18lf",x);
    lengthes[0]=len(x);
    for(i=0;i<N;i++)
    {
        fprintf(fw," %.18lf",y[i]);
        lengthes[i+1]=len(y[i]);
    }
    fprintf(fw,"\n");

    for(;;)
    {
        if(nstep>NMAX)
        {
            printf("Dopri: Слишком много шагов!\n");
            fclose(fw);
            return; 
        }
        if(!(x+0.03*h>x) && !(x+0.03*h<x))
        {
            printf("Dopri:  Слишком маленький шаг! t=%.18lf\n",x);
            printf("        Обращений к fcn: %d; Общее число шагов: %d; Принятых шагов: %d; Смен шага: %d\n",nfcn,nstep,naccpt,nrejct);
            fclose(fw);
            return;
        }
        if((x-xend)*posneg+EPS>0.) 
        {
            fclose(fw);
            sprintf(doprisays,"Dopri:  \"Обращений к fcn: %d; Общее число шагов: %d; Принятых шагов: %d; Смен шага: %d\"\n",nfcn,nstep,naccpt,nrejct);
            if((i=output(DOPRITMP,filename,lengthes,N+1,doprisays)))
                printf("Dopri: \"Ошибка вывода результата в файл!\"\n");                
            return;
        }
        if((x+h-xend)*posneg>0.)
            h=xend-x;
        fcn(x,y,k1);
        for(;;)
        {
            nstep++;
            for(i=0;i<N;i++)
                y1[i]=y[i]+h*a21*k1[i];
            fcn(x+c2*h,y1,k2);
            for(i=0;i<N;i++)
                y1[i]=y[i]+h*(a31*k1[i]+a32*k2[i]);
            fcn(x+c3*h,y1,k3);
            for(i=0;i<N;i++)
                y1[i]=y[i]+h*(a41*k1[i]+a43*k3[i]);
            fcn(x+c4*h,y1,k4);
            for(i=0;i<N;i++)
                y1[i]=y[i]+h*(a51*k1[i]+a53*k3[i]+a54*k4[i]);
            fcn(x+c5*h,y1,k5);
            for(i=0;i<N;i++)
                y1[i]=y[i]+h*(a61*k1[i]+a64*k4[i]+a65*k5[i]);
            fcn(x+c6*h,y1,k6);
            for(i=0;i<N;i++)
                y1[i]=y[i]+h*(a71*k1[i]+a74*k4[i]+a75*k5[i]+a76*k6[i]);
            fcn(x+c7*h,y1,k7);
            for(i=0;i<N;i++)
                y1[i]=y[i]+h*(a81*k1[i]+a84*k4[i]+a85*k5[i]+a86*k6[i]+a87*k7[i]);
            fcn(x+c8*h,y1,k2);
            for(i=0;i<N;i++)
                y1[i]=y[i]+h*(a91*k1[i]+a94*k4[i]+a95*k5[i]+a96*k6[i]+a97*k7[i]+a98*k2[i]);
            fcn(x+c9*h,y1,k3);
            for(i=0;i<N;i++)
                y1[i]=y[i]+h*(a101*k1[i]+a104*k4[i]+a105*k5[i]+a106*k6[i]+a107*k7[i]+a108*k2[i]+a109*k3[i]);
            for(i=0;i<N;i++) 
            {
                y11s=a111*k1[i]+a114*k4[i]+a115*k5[i]+a116*k6[i]+a117*k7[i]+a118*k2[i]+a119*k3[i];
                y12s=a121*k1[i]+a124*k4[i]+a125*k5[i]+a126*k6[i]+a127*k7[i]+a128*k2[i]+a129*k3[i];
                k4[i]=a131*k1[i]+a134*k4[i]+a135*k5[i]+a136*k6[i]+a137*k7[i]+a138*k2[i]+a139*k3[i];
                k5[i]=b1*k1[i]+b6*k6[i]+b7*k7[i]+b8*k2[i]+b9*k3[i];
                k6[i]=bh1*k1[i]+bh6*k6[i]+bh7*k7[i]+bh8*k2[i]+bh9*k3[i];
                k2[i]=y11s;
                k3[i]=y12s;
            }
            fcn(x+c10*h,y1,k7);
            for(i=0;i<N;i++)
                y1[i]=y[i]+h*(k2[i]+a1110*k7[i]);
            fcn(x+c11*h,y1,k2);
            xph=x+h;
            for(i=0;i<N;i++)
                y1[i]=y[i]+h*(k3[i]+a1210*k7[i]+a1211*k2[i]);
            fcn(xph,y1,k3);
            for(i=0;i<N;i++)
                y1[i]=y[i]+h*(k4[i]+a1310*k7[i]+a1311*k2[i]);
            fcn(xph,y1,k4);
            nfcn+=13;
            for(i=0;i<N;i++) 
            {
                k5[i]=y[i]+h*(k5[i]+b10*k7[i]+b11*k2[i]+b12*k3[i]+b13*k4[i]);
                k6[i]=y[i]+h*(k6[i]+bh10*k7[i]+bh11*k2[i]+bh12*k3[i]);
            }
            err=0.;
            for(i=0;i<N;i++) 
            {
                denom=max(max(1.e-6,fabs(k5[i])),max(fabs(y[i]),2.*EPS/eps));
                err+=((k5[i]-k6[i])/denom)*((k5[i]-k6[i])/denom);
            }
            err=sqrt(err/(double)(N));
            fac=max((1./6.),min(3.,pow(err/eps,(1./8.))/0.9));
            hnew=h/fac;
            if(err>eps)
            {
                reject=1;
                h=hnew;
                if(naccpt>1)
                    nrejct++;
                nfcn--;
                continue;
            }
            naccpt++;
            for(i=0;i<N;i++)
                y[i]=k5[i];
            x=xph;
            fprintf(fw,"%.18lf",x);
            lengthes[0]=max(lengthes[0],len(x));
            for(i=0;i<N;i++) 
            {
                fprintf(fw," %.18lf",y[i]);
                lengthes[i+1]=max(lengthes[i+1],len(y[i]));
            }
            fprintf(fw,"\n");   
            if(fabs(hnew)>hmax)
                hnew=posneg*hmax;
            if(reject==1)
                hnew=posneg*min(fabs(hnew),fabs(h));
            reject=0;
            h=hnew;
            break;
        }
    }
    
    return;
}


int len(double x)
{
    int length;
    double y=x,tmp;
    
    for(length=1,tmp=(x>0.?x:-x);tmp>10.-1000*EPS;length++,tmp*=0.1);
    
    return length+(y<0.?1:0);
}


int output(const char *in,const char *out,int *lengthes,int n,char *doprisays)
{
    int i,j,colnameslen=3;
    double x;
    FILE *fr,*fw;
    char *s,buf[2000],colnames[N+1][4]={" t "," x "," y "," u "," v ","Px ","Py ","Pu ","Pv "};
    
    if(!(fr=fopen(in,"r")))
        return -1;
    
    if(!(fw=fopen(out,"a")))
    {
        fclose(fr);
        return -2;
    }
    
    for(i=0;i<N+1;i++)
    {
        for(j=0;j<(lengthes[i]+19)/2-colnameslen/2;j++)
            fprintf(fw," ");
        fprintf(fw,"%s",colnames[i]);
        for(j=(lengthes[i]+19)/2-colnameslen/2+colnameslen;j<lengthes[i]+20;j++)
            fprintf(fw," ");
    
    }
    
    fprintf(fw,"\n");
    for(i=1;i<N+1;i++)
        lengthes[i]+=1;
    fgets(buf,2000,fr);
    fgets(buf,2000,fr);
    for(;;)
    {
        if(fgets(buf,2000,fr))
        {
            for(i=0,s=buf;i<N+1;i++)
            {
                sscanf(s,"%lf",&x);
                for(j=0;j<lengthes[i]-len(x);j++)
                    fprintf(fw," ");
                fprintf(fw,"%.18lf",x);
                s+=len(x)+20;
            }
            fprintf(fw,"\n");
        }
        else
            break;
    }

    for(i=1;i<N+1;i++)
        lengthes[i]-=1;
    for(i=0;i<N+1;i++)
    {
        for(j=0;j<(lengthes[i]+19)/2-colnameslen/2;j++)
            fprintf(fw," ");
        fprintf(fw,"%s",colnames[i]);
        for(j=(lengthes[i]+19)/2-colnameslen/2+colnameslen;j<lengthes[i]+20;j++)
            fprintf(fw," ");
    }
    fprintf(fw,"\n\n");

    fprintf(fw,"%s\n",doprisays);
    fclose(fr);
    fclose(fw);
    
    return 0;
}

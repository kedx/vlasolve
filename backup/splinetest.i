d=read_ascii("test_spline1D_X1.txt");
e=read_ascii("test_spline1D_X2.txt");
f=read_ascii("test_spline1DALL_X1.txt");
g=read_ascii("test_spline1DALL_X2.txt");

WS,1,dpi=150;
plg,g(3,),g(2,),width=4,color=__grey;
plg2,f(3,),f(2,),msize=10;

plg,e(3,),e(2,),color=__red;
plg,e(3+2,),e(2+2,),color=__blue;
plg2,d(3,),d(2,),msize=10;
plg2,d(3+2,),d(2+2,),msize=10;

w=[2.21431e-06,-1.77145e-05,7.97152e-05,-0.000301146,0.0011138,-0.00414519,0.0154647,-0.0577138,0.21539,-0.803848];

xb=d(2+2,);
yb=d(3+2,);
dfb=(w(::-1)*yb(2:12-1))(sum)/(xb(dif)(avg));

xa=d(2,);
ya=d(3,);
dfa=(w*ya(-10:-1))(sum)/(xa(dif)(avg));

dfa;dfb;
/*
xx=g(2,);
yy=g(3,);
yp=yy(dif)/xx(dif);
xp=xx(zcen);
yp=smooth(yp,20);
yp2=yp(dif)/xp(dif);
xp2=xp(zcen);
yp2=smooth(yp2,20);
yp3=yp2(dif)/xp2(dif);
xp3=xp2(zcen);
yp3=smooth(yp3,20);

WS,3,dpi=150;
plg,yp,xp,color=__red;
plg,yp2,xp2,color=__blue;
plg,yp3,xp3,color=__black;
*/
// derivative ...


RES=20;
FAC=88;
d=read_ascii("test_spline2D_X1.txt");
e=read_ascii("test_spline2D_X2.txt");
f=read_ascii("test_spline2DALL_X1.txt");
g=read_ascii("test_spline2DALL_X2.txt");

x=d(2,);y=d(3,);z=d(4,);
x2=e(2,);y2=e(3,);z2=e(4,);
dim=[2,RES+1,RES+1];
x=reform(x,dim);y=reform(y,dim);z=reform(z,dim);
dim=[2,RES*FAC+1,RES*FAC+1];
x2=reform(x2,dim);y2=reform(y2,dim);z2=reform(z2,dim);

xx=d(2+3,);yy=d(3+3,);zz=d(4+3,);
xx2=e(2+3,);yy2=e(3+3,);zz2=e(4+3,);
dim=[2,RES+1,RES+1];
xx=reform(xx,dim);yy=reform(yy,dim);zz=reform(zz,dim);
dim=[2,RES*FAC+1,RES*FAC+1];
xx2=reform(xx2,dim);yy2=reform(yy2,dim);zz2=reform(zz2,dim);

xa=f(2,);ya=f(3,);za=f(4,);
xa2=g(2,);ya2=g(3,);za2=g(4,);

dim=[2,2*RES+1,RES+1];//dim=[2,RES+1,2*RES+1];
xa=reform(xa,dim);ya=reform(ya,dim);za=reform(za,dim);
//xa=transpose(xa);ya=transpose(ya);za=transpose(za);

dim=[2,2*RES*FAC+1,RES*FAC+1];//dim=[2,RES*FAC+1,2*RES*FAC+1];
xa2=reform(xa2,dim);ya2=reform(ya2,dim);za2=reform(za2,dim);
//xa2=transpose(xa2);ya2=transpose(ya2);za2=transpose(za2);

WS,6,dpi=150;

for (i=1;i<=20;i++)
  {
    plg,za2(,1+FAC*(i-1)),xa2(,1+FAC*(i-1)),width=4,color=__grey;
    plg2,za(,i),xa(,i),msize=10;
    
    plg,z2(,1+FAC*(i-1)),x2(,1+FAC*(i-1)),color=__red;
    plg,zz2(,1+FAC*(i-1)),xx2(,1+FAC*(i-1)),color=__blue;
    plg2,z(,i),x(,i),color=__red,msize=10;
    plg2,zz(,i),xx(,i),color=__blue,msize=10;
  }
/*
WS,7,dpi=150;

for (i=1;i<=20;i++)
  {
    plg,za2(1+FAC*(i-1),),ya2(1+FAC*(i-1),),width=4,color=__grey;
    plg2,za(i,),ya(i,),msize=10;
    
    plg,z2(1+FAC*(i-1),),y2(1+FAC*(i-1),),color=__red;
    plg,zz2(1+FAC*(i-1),),yy2(1+FAC*(i-1),),color=__blue;
    plg2,z(i,),y(i,),color=__red,msize=10;
    plg2,zz(i,),yy(i,),color=__blue,msize=10;
  }
*/
z2(1)=zz2(1)=min(za2);
z2(2)=zz2(2)=max(za2);

WS,2,dpi=150;
pli,z,min(x),min(y),max(x),max(y);c_c;
pli,zz,min(xx),min(yy),max(xx),max(yy);c_c;

WS,3,dpi=150;
pli,za,min(xa),min(ya),max(xa),max(ya);c_c;

WS,4,dpi=150;
pli,z2,min(x2),min(y2),max(x2),max(y2);c_c;
pli,zz2,min(xx2),min(yy2),max(xx2),max(yy2);c_c;

WS,5,dpi=150;
pli,za2,min(xa2),min(ya2),max(xa2),max(ya2);c_c;



WS,10,dpi=150;
u=span(0.,2*2*pi,800);
v=span(0,2*pi,800);
uu=u(,-::799);
vv=v(-::799,);
uv=[uu,vv];
w=cos(uv(,,1))+sin(uv(,,1)*0.84+uv(,,2)*0.46-0.13)*cos(uv(,,1)*0.56-uv(,,2));
//w=cos(uv(,,1))+cos(uv(,,1)-uv(,,2));
pli,w,min(u),min(v),max(u),max(v);
c_c;


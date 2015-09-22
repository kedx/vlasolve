fname_stat="UDF_200_200_224_0.005_stats/snapshot";
fname_ref="UDF_reflective/snapshot";
fstats="UDF_200_200_224_0.005_stats/statistics.txt";
fstats_half="UDF_200_200_224_0.0025_stats/statistics.txt";
step=50.;
count=60;

quickmovie,2.5,80,cc=1,fname=fname_stat;
quickmovie,2.5,80,dim=2,j=100,cc=1,fname=fname_stat;
// Problem reflective kernel
d=NDfield_read(fname_stat+"_20.0000.ND",hd);
WS,9,dpi=150;pliND,d,hd,1;c_c;logxy,1,0;
quickmovie,1,40,cc=1,dim=2,j=100,fname=fname_ref,pt=200;
// Convergence rho
quickmovie,step,count,cc=1,plrho=1,fname=fname_stat;
//kernel problem
quickmovie,step,count,cc=1,dim=2,j=100,maxval=0.01,fname=fname_stat;
quickmovie,step,count,cc=1,dim=2,j=100,maxval=0.001,fname=fname_stat;
// stats
plotStats,[fstats,fstats_half];
//kernel problem
d=NDfield_read(fname_stat+"_2950.0000.ND",hd);
WS,10,dpi=150;pliND,d,hd,18;c_c;logxy,1,0;
d=NDfield_read(fname_stat+"_30.0000.ND",hd);
WS,11,dpi=150;pliND,d,hd,18;c_c;logxy,1,0;

func plotStats(fname,plType=,nofma=,dpi=)
{
  if (is_void(fname)) fname="statistics.txt";
  if (is_void(plType)) plType=1;
  if (is_void(nofma)) nofma=0;
  if (is_void(dpi)) dpi=150;
  
  if (numberof(fname)>1)
    {
      for (i=1;i<=numberof(fname);i++) plotStats,fname(i),plType=plType+i-1,nofma=(i-1)+nofma;
      return;
    }
  
  d=read_ascii(fname);
  win=0;
  alpha = 1.-(plType-1)*0.3;
  
  if (!nofma) WS,win++,dpi=dpi;
  else window,win++;
  
  plg,d(-1,),d(1,),width=2,color=char(__red*alpha),type=plType;
  plg,d(3,),d(1,),width=2,color=char(__blue*alpha),type=plType;

  if (!nofma) WS,win++,dpi=dpi;
  else window,win++;
 
  plg,-d(-2,),d(1,),width=2,color=char(__red*alpha),type=plType;

  if (!nofma) WS,win++,dpi=dpi;
  else window,win++;

  plg,-d(-3,),d(1,),width=2,color=char(__purple*alpha),type=plType;
  plg,d(-4,),d(1,),width=2,color=char(__cyan*alpha),type=plType;
  plg,-d(-5,),d(1,),width=2,color=char(__red*alpha),type=plType;
}


/*
WS,0,dpi=150;
fma;
a=read_ascii("statistics.txt");
b=read_ascii("statistics_0.0025.txt");
d=a;
i=0;

plg,d(i--,),d(1,),color=__orange,width=3;
plg,d(i--,),d(1,),color=__green,width=3;
plg,-d(i--,),d(1,),color=__purple,width=3;
plg,d(i--,),d(1,),color=__red,width=3;
plg,d(i--,),d(1,),color=__blue,width=3;
plg,d(i--,),d(1,),color=__black,width=3;
plg,d(i--,),d(1,),color=__black,width=3,type=2;
d=b;i=0;
plg,d(i--,),d(1,)
plg,d(i--,),d(1,)
plg,-d(i--,),d(1,)
plg,d(i--,),d(1,)
plg,d(i--,),d(1,)
plg,d(i--,),d(1,)
plg,d(i--,),d(1,),type=2;
*/

func genScale(x0,xmax,N,sclt,valt,delta=)
{
  if (sclt==0) r=span(x0,xmax,N+(valt!=0));
  else if (sclt==1) r=exp(span(log(x0),log(xmax),N+(valt!=0)));
  else r=pow(span(pow(x0,1./sclt),pow(xmax,1./sclt),N+(valt!=0)),sclt);

  if (!is_void(delta)) return r(dif);
  
  if (valt) r=r(zcen);
 
  return r;
}


func getScales(dims,header,x0=,xmax=)
{
  //if (h.comment=="no comment")
  //if (is_void(header))
    {
      scl=[1,0,2];
      val=[0,0,1];
    }
  
    if (is_void(x0)) x0=*header.x0;
    if (is_void(xmax))
      {
        if (!is_void(delta)) xmax=x0+*header.delta;
        else xmax=x0+*header.delta;
      }
    
    result=[];
    for (i=1;i<=numberof(scl);i++) grow,result,&genScale(x0(i),xmax(i),dims(1+i),scl(i),val(i));
    for (i=1;i<=numberof(scl);i++) grow,result,&genScale(x0(i),xmax(i),dims(1+i),scl(i),val(i),delta=1);
    
    return result;
}


func density(d,h=,scl=)
{
  if (is_void(scl))
    {
      if (is_void(h)) error,"density: I need the header or the scales !!!";
      scl=getScales(dimsof(d),h)
    }
  
  dim=dimsof(d)(2:);
  d1=dim(1)-1;
  d2=dim(2)-1-1;
  d3=dim(3)-1;
  
  r=*scl(1);
  du=*scl(3+2);
  j=*scl(3);
  dj=*scl(3+3);

  du=du(-::d1,,-::d3);
  j=j(-::d1,-::d2,);
  dj=dj(-::d1,-::d2,);
  
  rho = ((2*pi*j*dj)*du*d(,zcen,))(,sum,sum);
  
  return [(rho/(r*r))(zcen),r(zcen)];
}

func pliNDmulti(d,h,N,dim=,minval=,maxval=,noyscale=,cmin=,cmax=,sliceAt=)
{
  if (is_void(dim)) dim=3;
  if (dim==1) {e=d(N,,..);}
  if (dim==2) {e=d(,N,..);}
  if (dim==3) {e=d(,,N,..);}
  
  if (is_void(minval)) minval=min(e(*));
  if (is_void(maxval)) maxval=max(e(*));
  for (i=1;i<=dimsof(d)(0);i++) pliND,d(..,i),h(..,i),N,minval=minval,maxval=maxval,noyscale=noyscale,cmin=cmin,cmax=cmax,sliceAt=sliceAt;
}

func getSliceAt(d,h,sliceAt)
{
  N=1;
  scl=*getScales(dimsof(d),h)(dim);
  w=where(scl<=sliceAt)(max);
  if (!is_void(w))
    {
      write,format=" Slice at U[%d](%d) = %e\n",dim,w(1),scl(w(1));
      N=w(1);
    }
  return N;
}

func pliND(d,h,N,dim=,minval=,maxval=,noyscale=,cmin=,cmax=,printVal=,sliceAt=,useLog=)
{
  x0=*h.x0;
  x1=x0+*h.delta;
  
  dm= str2long(strtok(h.comment," ",6));
  if (is_void(dim)) dim=3;
  if (is_void(N)) N=-1;
  /*
  if (!is_void(sliceAt))   
    {
      scl=*getScales(dimsof(d),h)(dim);
      w=where(scl<=sliceAt)(max);
      if (!is_void(w))
        {
          write,format=" Slice at U[%d](%d) = %e\n",dim,w(1),scl(w(1));
          N=w(1);
        }
    }
  */
  if (!is_void(sliceAt))
    {
      N = getSliceAt(d,h,sliceAt);
      if (is_void(printVal)) printVal=1;
    }
  
  if ((!is_void(printVal))&&(N>=0))
    {
      scl=*getScales(dimsof(d),h)(dim);
      write,format=" Slice at U[%d](%d) = %e\n",dim,N,scl(N);
    }
  
  if (N<0)
    {
      if (is_void(h)) error,"density: I need the header or the scales !!!";
      scl=getScales(dimsof(d),h);
      dms=dimsof(d)(2:);
      d1=dms(1)-1;
      d2=dms(2)-1;
      d3=dms(3)-1;
  
      r=*scl(1);
      //du=*scl(3+2);
      
      dj=*scl(3+3);

      //du=du(-::d1,,-::d3);
      j=j(-::d1,-::d2,);
      dj=dj(-::d1,-::d2,);
      if (dim==3) {
        e=(2*pi*j*dj*d)(,,sum);
        e/=(r*r);
        x0=x0([1,2]);x1=x1([1,2]);dm=dm([1,2,3,4]);
      }
    }
  else
    {
      if (dim==1) {e=d(N,,);x0=x0([2,3]);x1=x1([2,3]);dm=dm([3,4,5,6]);}
      if (dim==2) {e=d(,N,);x0=x0([1,3]);x1=x1([1,3]);dm=dm([1,2,5,6]);}
      if (dim==3) {e=d(,,N);x0=x0([1,2]);x1=x1([1,2]);dm=dm([1,2,3,4]);}
    }
  dm;
  for (i=1;i<=2;i++)
    {
      n=(i-1)*2 +1;
      if (dm(n+1)==0) continue;
      if (dm(n)==0)
        {
          delta= 0.5 * (x1(i)-x0(i))/(dimsof(e)(2:)(i)-1);
          x0(i)-=delta;x1(i)+=delta;
        }
      else if (dm(n)==1)
        {
          delta=0.5*(log(x1(i))-log(x0(i)))/(dimsof(e)(2:)(i)-1);
          //x0(i);x1(i);
          x0(i)=exp(log(x0(i))-delta);x1(i)=exp(log(x1(i))+delta);
          //x0(i);x1(i);
        }
      else
        {
          delta= 0.5 * (pow(x1(i),1./dm(n))-pow(x0(i),1./dm(n)))/(dimsof(e)(2:)(i)-1);
          x0(i)=pow( pow(x0(i),1./dm(n))-delta ,dm(n));
          x1(i)=pow( pow(x1(i),1./dm(n))-delta ,dm(n));
        }
    }
     
  if (!is_void(minval)) {w=where(e<minval);if (numberof(w)) e(w)=minval; else e(*)(1)=minval;}
  if (!is_void(maxval)) {w=where(e>maxval);if (numberof(w)) e(w)=maxval; else e(*)(2)=maxval;}
  
  if (is_void(noyscale))
    {
       if (is_void(useLog))
         pli,e,x0(1),x0(2),x1(1),x1(2);
       else
         pli,log(1.0+useLog*e),x0(1),x0(2),x1(1),x1(2);
    }
  else
    {
      if (is_void(useLog))
        pli,e,x0(1),0,x1(1),dimsof(e)(3),cmin=cmin,cmax=cmax;
      else
        pli,log(1.0+useLog*e),x0(1),0,x1(1),dimsof(e)(3),cmin=cmin,cmax=cmax;
    }

  return N;
  /*
  info,e;info,x0;
  plg,e(,1),*scl(1);
  smwrite("profile.dat",[e(1,),*scl(2)]);
  */
}

func pluv(muv,j=,color=,w=)
{
  if (is_void(j)) j=1;
  if (is_void(w)) pldj,muv(1,j,),muv(2,j,),muv(3,j,),muv(4,j,),color=color;
  else  pldj,muv(1,j,w),muv(2,j,w),muv(3,j,w),muv(4,j,w),color=color;
}

// quickmovie,1.0,100,j=137,cc=1,lmts=[0.1,25,-1.5,1.5],jpg=1,minval=0,maxval=0.065
func quickmovie(dt,N,fname=,pt=,j=,dpi=,win=,plrho=,i0=,nopli=,pldif=,dim=,minval=,maxval=,logx=,logy=,cc=,noyscale=,lmts=,jpg=,sliceAt=,useLog=)
{
  if (is_void(fname)) fname ="snapshot";
  if (is_void(pt)) pt=20;
  if (is_void(j)) j=-1;
  if (is_void(dpi)) dpi=150;
  if (is_void(i0)) i0=0.;
  if (is_void(dim)) dim=3;
  if (is_void(logx)) logx=1;
  if (is_void(logy)) logy=0;  

  nfn=numberof(fname);

  col=[__black,__red,__blue,__green];
  
  //name=array(string,nfn);
  //dp=array(pointer,nfn);
  sclp=array(pointer,nfn);
  dref=[];
  for (ct=1;ct<=nfn;ct++)
    {
      name=swrite(format=fname(ct)+"_%4.4f.ND",(i0+0)*dt);
      d=NDfield_read(name,h);
      if (!is_void(sliceAt)) j=getSliceAt(d,h,sliceAt);
      if (!is_void(pldif)) grow,dref,&d;
      sclp(ct)=&getScales(dimsof(d),h);
      write,format=" Slice at U[%d](%d) = %e\n",dim,j,(*(*sclp(ct))(dim))(j);
    }
  
  
  if (is_void(win)) win=1;
  
  for (ct=1;ct<=nfn;ct++) {
    WS,win+ct-1,dpi=dpi;logxy,logx,logy;
    if (!is_void(cc)) c_c;
    if (is_void(jpg)) animate,1; 
    if (!is_void(lmts)) limits,lmts(1),lmts(2),lmts(3),lmts(4);
  }
  if (!is_void(plrho)) {
    WS,win+nfn,dpi=dpi;
    
    scl=*sclp(1);
    limits,(*scl(1))(min),(*scl(1))(max),1.E-8,10;
    logxy,1,1;animate,1;
  }

  dimsofd = dimsof(d);
  //(*(*sclp(1))(3));
  for (i=0;i<N;i++) {
    
    for (ct=1;ct<=nfn;ct++)
      {
        name=swrite(format=fname(ct)+"_%4.4f.ND",(i0+i)*dt);
        
        if ((dim==3)&&(is_void(plrho))) d=NDfield_read(name,h,first=1+(j-1)*dimsofd(2)*dimsofd(3),last=1+j*dimsofd(2)*dimsofd(3));
        else d=NDfield_read(name,h);
        
        scl=*sclp(ct);
        if (is_void(nopli)) {
          window,win+ct-1;
          if (is_void(pldif)) pliND,d,h,j,dim=dim,minval=minval,maxval=maxval,noyscale=noyscale,useLog=useLog;
          else pliND,d-(*dref(ct)),h,j,dim=dim,minval=minval,maxval=maxval,noyscale=noyscale,useLog=useLog;
          //plg,d(,int(dimsof(d)(3)/2),j)*10,spanl((*scl(1))(1),(*scl(1))(0),dimsof(d)(2)),color=__blue;
          add_txt,swrite(format="T = %3.3f",(i0+i)*dt),pos=[0.45,0.82];
          if (dim==3) add_txt,swrite(format="J = %4.3e",(*scl(3))(j)),pos=[0.45,0.8];
        }
        if (!is_void(jpg)) {
          if (!is_string(jpg)) jpg="SNAP";
          jpeg,swrite(format="%s_%5.5d.jpg",jpg,i),nodisp=1,background=1;
        }
        
        if (!is_void(plrho))
          {
            window,win+nfn;
            
            rho=density(d,scl=scl);
            plg,rho(,1),rho(,2),color=col(,ct);
            if (ct==1) x=spanl(0.5,20,500);
            if (ct==1) {plg,0.04*pow(x,-4),x,color=__red;plg,pow(x,-6),x,color=__blue;}
            if (ct==1) {add_txt,swrite(format="T = %3.3f",(i0+i)*dt),pos=[0.45,0.82];}
            if (!is_void(jpg)) {
              if (!is_string(jpg)) jpg="SNAP";
              jpeg,swrite(format="%s_density_%5.5d.jpg",jpg,i),nodisp=1,background=1;
            }
          }
    
        //window,win+ct-1;fma;
        //if (!is_void(plrho)) {if (ct==nfn) {window,win+nfn;fma;} }
      }
    for (ct=1;ct<=nfn;ct++) {window,win+ct-1;fma;}
    if (!is_void(plrho)) {window,win+nfn;fma;}
    pause,pt;
  }
  i-=1;
  for (ct=1;ct<=nfn;ct++)
    {
      WS,win+ct-1,dpi=dpi;//animate,0;fma;
      logxy,logx,logy;
      if (!is_void(cc)) c_c;
      
      name=swrite(format=fname(ct)+"_%4.4f.ND",(i0+i)*dt);
      if (dim==3) d=NDfield_read(name,h,first=1+(j-1)*dimsofd(2)*dimsofd(3),last=1+j*dimsofd(2)*dimsofd(3));
      else d=NDfield_read(name,h);
      
      if (is_void(pldif)) pliND,d,h,j,dim=dim,minval=minval,maxval=maxval,noyscale=noyscale,useLog=useLog;
      else pliND,d-(*dref(ct)),h,j,dim=dim,minval=minval,maxval=maxval,noyscale=noyscale,useLog=useLog;
   
      scl=*sclp(ct);
      if (ct==1) add_txt,swrite(format="T = %3.3f",(i+i0)*dt),pos=[0.45,0.82];
      if (ct==1) add_txt,swrite(format="J = %4.3e",(*scl(3))(j)),pos=[0.45,0.8];

      if (!is_void(plrho))
        {
          if (ct==1) WS,win+nfn,dpi=dpi;//fma;
          else window,win+nfn;
          rho=density(d,scl=scl);
          logxy,1,1;
          plg,rho(,1),rho(,2),color=col(,ct);
          if (ct==1) x=spanl(0.5,20,500);
          if (ct==1) {plg,0.04*pow(x,-4),x,color=__red;plg,pow(x,-6),x,color=__blue;}
          if (ct==1) {add_txt,swrite(format="T = %3.3f",(i+i0)*dt),pos=[0.45,0.82];}
        }
    }
  
    
}

#ifndef __MPI_REQUESTS_HXX__
#define __MPI_REQUESTS_HXX__

#include "mpiCommunication.hxx"

//typename gridSliceItT::pointer
template <typename ptrT, class mpiComT, int DIMS>
class mpiValAtCoordsRequestT
  {
  public:    
    typedef mpiValAtCoordsRequestT<ptrT,mpiComT,DIMS> myT;

    std::vector<int> procID;
    std::vector<ptrT> ptr;
    std::vector<double> coords;
    std::vector<double> val;
    int nfields;
    int nProcs;
    std::vector< std::pair<long,long> > range;    
    
    void clear()
    {
      procID.clear();
      ptr.clear();
      coords.clear();
      range.clear();
      val.clear();
      nProcs=0;
    }

    struct compareLess
    {
      const std::vector<int> &t;
      compareLess(const std::vector<int> &tab):t(tab){}
      bool operator() (const long a,const long b) { return t[a]<t[b];}
    };
   
    void setupLocal(int nfields_=1)
    {
      nfields=nfields_;
      if (procID.size()==0) 
	{
	  nProcs=0;
	  return;
	}

      std::vector<long> srtID;
      if (srtID.size()!=procID.size())
	{
	  srtID.resize(procID.size());
	  for (long j=0;j<procID.size();j++) srtID[j]=j;
	}

      compareLess cmp(procID);
      std::sort(srtID.begin(),srtID.end(),cmp);
      for (long j=0;j<procID.size();j++)
	{
	  unsigned int k=srtID[j];
	  
	  std::swap(ptr[j],ptr[k]);
	  std::swap(procID[j],procID[k]);
	  for (int i=0;i<DIMS;i++) std::swap(coords[DIMS*j+i],coords[DIMS*k+i]);

	  std::swap(srtID[j],srtID[k]);
	}

      range.clear();
      nProcs=0;
      long ref=0;
      for (long i=0;i<procID.size();i++)
	{
	  if (procID[i]!=procID[ref])
	    {
	      nProcs++;
	      range.push_back(std::make_pair(ref,i));
	      ref=i;
	    }
	} 
      nProcs++;
      
      range.push_back(std::make_pair(ref,procID.size()));
      for (long i=0;i<range.size();i++) procID[i]=procID[range[i].first];
      procID.resize(range.size());
      val.assign(ptr.size()*nfields,0);
    }
    
    
    void setupRemote(myT& myRequests,mpiComT *mpiCom)
    {
      int NMAX=2*hlp::IntPower<3,DIMS>::value+1;
      std::vector<int> buffer;
      //MPI_Status status;
      long i,j;

      nfields=myRequests.nfields;

      for (i=0;i<mpiCom->size();i++)
	{
	  buffer.assign(NMAX,-1);
	
	  if (i==mpiCom->rank()) 
	    {
	      for (j=0;j<myRequests.procID.size();j++)
		{
		  buffer[2*j]=myRequests.procID[j];
		  buffer[2*j+1]=myRequests.range[j].second-myRequests.range[j].first;
		}
	    }
	
	  //MPI_Bcast(&buffer[0],buffer.size(),MPI_INT,i,MPI_COMM_WORLD);
	  mpiCom->Bcast(buffer,i);

	  if (i!=mpiCom->rank())
	    {
	      j=0;
	      while (buffer[j]>-1)
		{
		  if (buffer[j++]==mpiCom->rank()) 
		    {
		      procID.push_back(i);
		      long start=coords.size()/DIMS;
		      coords.resize(DIMS*(start+buffer[j]));
		      //MPI_Recv(&coords[DIMS*start],DIMS*buffer[j],MPI_DOUBLE,i,0,MPI_COMM_WORLD, &status);
		      mpiCom->Recv(&coords[DIMS*start],DIMS*buffer[j],i,0);
		      range.push_back(std::make_pair(start,start+buffer[j]));
		      nProcs++;
		    }
		  j++;
		};
	    }
	  else
	    {
	      MPI_Request mpiReq[myRequests.procID.size()];
	    
	      for (j=0;j<myRequests.procID.size();j++)
		{
		  long size=DIMS*(myRequests.range[j].second-myRequests.range[j].first);
	
		  //MPI_Isend(&myRequests.coords[DIMS*myRequests.range[j].first],size,MPI_DOUBLE,myRequests.procID[j],0,MPI_COMM_WORLD,&mpiReq[j]);
		  mpiCom->Isend(&myRequests.coords[DIMS*myRequests.range[j].first],size,myRequests.procID[j],&mpiReq[j],0);
		}
	   
	      for (j=0;j<myRequests.procID.size();j++)
		{
		  mpiCom->Wait(&mpiReq[j]);
		  //MPI_Wait(&mpiReq[j], &status);
		}
	    }
	}
    
      if (nProcs) val.resize(range.back().second*nfields);   
    }
    
    // req is the neighbor request (on the remote processes)
    void synchronize(myT& req, mpiComT *mpiCom)
    {
      	MPI_Request mpiReqS[req.nProcs];
	MPI_Request mpiReqR[nProcs];
	long i;

	for (i=0;i<req.nProcs;i++)
	  {
	    long size=nfields*(req.range[i].second-req.range[i].first);
	    //printf ("Proc %d sending %ld to %d\n",bT::mpiCom->rank(),size,req.procID[i]);
	    //MPI_Isend(&req.val[nfields*req.range[i].first],size,MPI_DOUBLE,req.procID[i],0,MPI_COMM_WORLD,&mpiReqS[i]);
	    mpiCom->Isend(&req.val[nfields*req.range[i].first],size,req.procID[i],&mpiReqS[i],0);
	  }

	for (i=0;i<nProcs;i++)
	  {
	    //printf("Proc %d say hello %ld %ld %ld = %ld\n",bT::mpiCom->rank(),range[0].second,range[0].first,range[0].second-range[0].first,val.size());
	    long size=nfields*(range[i].second-range[i].first);
	    //printf ("Proc %d receiving %ld from %d\n",bT::mpiCom->rank(),size,procID[i]);
	    //MPI_Irecv(&val[nfields*range[i].first],size,MPI_DOUBLE,procID[i],0,MPI_COMM_WORLD,&mpiReqR[i]);
	    mpiCom->Irecv(&val[nfields*range[i].first],size,procID[i],&mpiReqR[i],0);
	  }
    
	MPI_Status status;

	for (i=0;i<nProcs;i++)
	  mpiCom->Wait(&mpiReqR[i]);
	//MPI_Wait(&mpiReqR[i], &status);
    
	for (i=0;i<val.size();i++)
	  (*ptr[i]) = val[i];

	for (i=0;i<req.nProcs;i++)
	  mpiCom->Wait(&mpiReqS[i]);
	//MPI_Wait(&mpiReqS[i], &status);  
	//printf("Proc %d all done !\n",bT::mpiCom->rank());
	
    }

};

#endif

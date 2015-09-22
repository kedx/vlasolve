#ifndef __PARAMS_PARSER_HXX__
#define __PARAMS_PARSER_HXX__

#include "global.h"

#include <string>
#include <sstream> 
#include <cstddef>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <list>

#include <stdlib.h>
#include <stdio.h>


class paramsParser {
public:
  static const std::string defaultCategory() {return "default";}

private:
  template <class T>
  bool from_string(T& t, const std::string& s) const
  {
    std::istringstream iss(s);
    return !(iss >> t).fail();
  }
 
  template <typename Container>
  Container& split(Container& result,
		   const typename Container::value_type& s,
		   const typename Container::value_type& delimiters) const
  {
    result.clear();
    size_t current;
    size_t next = -1;
    do
      {
	next = s.find_first_not_of( delimiters, next + 1 );
	if (next == Container::value_type::npos) break;
	next -= 1;
	current = next + 1;
	next = s.find_first_of( delimiters, current );
	result.push_back( s.substr( current, next - current ) );
      }
    while (next != Container::value_type::npos);
    return result;
  }

  struct paramValueT {
    std::vector< std::vector<std::string> > val;
    std::vector< std::string > fname;
    std::vector< int > line;
    std::vector< bool > wasRead;
    //paramValue():wasRead(false){}
  };

  typedef std::map< std::string, paramValueT > paramsT;
  typedef paramsT::iterator paramsT_it;
  typedef paramsT::const_iterator paramsT_cit;

  typedef std::list< paramsT > paramsListT;
  typedef paramsListT::iterator paramsListT_it;
  typedef paramsListT::const_iterator paramsListT_cit;

  typedef std::map<std::string,paramsListT_it> categoryT;
  typedef std::map<std::string,paramsListT_it>::iterator categoryT_it;
  typedef std::map<std::string,paramsListT_it>::const_iterator categoryT_cit;

  std::string execFileName;
  paramsListT params;
  categoryT category;
  std::string paramFileName;

  std::pair<std::string,std::string> getKey(const std::string &rawKey)
  {
    std::pair<std::string,std::string> result;
    size_t p;
    p=rawKey.find_first_of(std::string("."));
    //printf("key=%s\n",rawKey.c_str());
    if (p==std::string::npos)
      return std::make_pair(defaultCategory(),rawKey);

    if ((p!=rawKey.find_last_of(std::string(".")))||(p==0)||(p==rawKey.length()-1))
      {
	fprintf(stderr,"ERROR in paramsParseT.\n");
	fprintf(stderr,"  ill formated keyword '%s'.\n",rawKey.c_str());
	exit(-1);
      }

    return std::make_pair(rawKey.substr(0,p),rawKey.substr(p+1));
  }

  void insertItem(const std::string &rawKey,const std::vector<std::string> &value)
  {
    std::pair<std::string,std::string> key=getKey(rawKey);

    categoryT_it cit=category.find(key.first);
    paramsListT_it lit;
    
    if (cit==category.end())
      {
	//printf("paramsParser: new category '%s'\n",key.first.c_str());
	lit=params.insert(params.begin(),paramsT());
	category.insert(std::make_pair(key.first,lit));
      }
    else lit=cit->second;

    paramsT_it it = lit->find(key.second);

    if (it==lit->end())
      {
	paramValueT p;
	it = lit->insert(std::make_pair(key.second,p)).first;	
      }
    
    it->second.val.push_back(value);
    it->second.wasRead.push_back(false);
    //it->second.val.push_back(fname);
    //it->second.val.push_back(line);

    /*
    std::pair<paramsT_it,bool> lit->insert(std::make_pair(key.second,value));

    if (! lit->insert(std::make_pair(key.second,value)).second)
      {
	//printf("WARNING: multiple definitions of '%s' ignored.\n",rawKey.c_str());
      }
    */
  }

  bool parse(std::string fname)
  {
    std::ifstream t(fname.c_str());

    if (!t) return false;

    std::string line;
    std::vector<std::string> tokens;

    std::string key;
    std::vector<std::string> value;

    while (!t.eof())
      {
	tokens.clear();
	std::getline(t,line);

	//printf("line:  %s",line.c_str());
	if (line.find_first_of("#")!=std::string::npos)
	  line=line.substr(0,line.find_first_of("#"));	  
	//printf("-> %s\n",line.c_str());

	split(tokens,line,std::string(" =:{}[],\t"));

	if (tokens.size())
	  {
	    std::string key=tokens[0];
	    std::vector<std::string> value(tokens.begin()+1,tokens.end());
	    insertItem(key,value);
	  }
      }
    return true;
  }

  bool parse(int argc, char **argv)
  {
    int n=0;
    std::string key;
    std::vector<std::string> value;

    while (n<argc)
      {
	if (argv[n][0]!='-')
	  {
	    fprintf(stderr,"ERROR parsing command line.\n");
	    fprintf(stderr,"  AT: %s\n",argv[n]);
	    exit(-1);
	  }
	
	value.clear();
	key=std::string(&argv[n++][1]);
	
	while ((n<argc)&&((argv[n][0]!='-')||(isdigit(argv[n][1])))) value.push_back(argv[n++]);
	
	insertItem(key,value);
      }

    return true;
  }

public:
  paramsParser(int argc=0, char **argv=NULL,std::string fname=std::string("params.ini"))
  {
    int i0=0;

    if (argc>i0)
      execFileName=std::string(argv[i0++]);

    if (argc>i0)
      {
	if (argv[i0][0]!='-') 
	  fname=std::string(argv[i0++]);
      }
    //printf("Parsing file : %s\n",fname.c_str());
    parse(argc-i0, &argv[i0]);
    if (parse(fname)) 
      paramFileName=fname;
    else if (fname!=std::string("params.ini"))
      {
	fprintf(stderr,"ERROR in paramsParser:");
	fprintf(stderr,"could not read parameter file '%s'.\n",fname.c_str());
	exit(0);
      }
    else paramFileName=std::string();
  }

  ~paramsParser()
  {

  }

  void report()
  {
    long i,j;

    printf("The following parameters will be used:\n");
    printf(" Input file: '%s'\n",paramFileName.c_str());
    for (categoryT_it cit=category.begin();cit!=category.end();cit++)
      {
	printf(" Category: '%s'\n",cit->first.c_str());
	for (paramsT_it pit=cit->second->begin();pit!=cit->second->end();pit++)
	  {
	    printf("   %s = ",pit->first.c_str());
	    for (j=0;j<pit->second.val.size();j++)
	      {
		if (j!=0) printf("      = ");
		for (i=0;i<pit->second.val[j].size();i++)
		  {
		    printf("%s ",pit->second.val[j][i].c_str());
		  }
		if (i==0) printf("true\n");
		else printf("\n");
	      }
	  }
      }
  }
  
  void reportUnused()
  {
    bool found=false;
    long j,i;

    for (categoryT_it cit=category.begin();cit!=category.end();cit++)
      {
	for (paramsT_it pit=cit->second->begin();pit!=cit->second->end();pit++)
	  {
	    for (j=0;j<pit->second.wasRead.size();j++)
	      {
		if (!pit->second.wasRead[j])
		  {
		    if (!found)
		      {
			printf("WARNING! The following parameter definitions were ignored:\n");
			found=true;
		      }
		    printf("   %s.%s = ",cit->first.c_str(),pit->first.c_str());
		    for (i=0;i<pit->second.val[j].size();i++)
		      {
			printf("%s ",pit->second.val[j][i].c_str());
		      }
		    if (i==0) printf("true\n");
		    else printf("\n");
		  }
	      }
	  }
      }
  }

  template <typename T>
  T getOrDie(const int index, const std::string &what,const std::string &cat,const int which=-1) const
  {
    categoryT_cit cit=category.find(cat);
    if (cit==category.end())
      {
	fprintf(stderr,"ERROR in paramParserT::getOrDie\n");
	fprintf(stderr," Category '%s' not found for keyword '%s'\n",cat.c_str(),what.c_str());
	exit(-1);
      }

    paramsT &p=*(cit->second);
    paramsT_it it=p.find(what);

    if (it==p.end())
      {
	fprintf(stderr,"ERROR in paramParserT::getOrDie\n");
	fprintf(stderr," keyword '%s' not found in category '%s'.\n",what.c_str(),cat.c_str());
	exit(-1);
      }

    if ((index>=0)&&(it->second.val.size()<=index))
      {
	fprintf(stderr,"ERROR in paramParserT::getOrDie\n");
	fprintf(stderr," '%dth' definition of keyword '%s.%s' not available (%ld total).\n",index,cat.c_str(),what.c_str(),it->second.val.size());
	exit(-1);
      }

    int id=(index<0)?0:index;

    if ((which>=0)&&(it->second.val[id].size()<=which))
      {
	fprintf(stderr,"ERROR in paramParserT::getOrDie\n");
	fprintf(stderr," for keyword '%s' in category '%s'.\n",what.c_str(),cat.c_str());
	fprintf(stderr," cannot retrieve value number %d (%ld found).\n",which,it->second.val[0].size());
	exit(-1);
      }
 
    int w=(which<0)?0:which;
    T result;
    if (!it->second.val[id].size()) result=true;
    else from_string<T>(result,it->second.val[id][w]);

    it->second.wasRead[id]=true;

    if (verbose>=2)
      {
	if (which<0) std::cout<<"New value assigned to "<<cat<<"."<<what<<" : "<<result<<".\n";
	else std::cout<<"New value assigned to "<<cat<<"."<<what<<"["<<which<<"] : "<<result<<".\n";
      }
    return result;   
  }
  
  template <typename T>
  T getOrDie(const std::string &what,const std::string &cat,const int which=-1) const
  {
    return getOrDie<T>(-1,what,cat,which);
  }

  template <typename T>
  T get(int index, const std::string &what,const std::string &cat,const T &def, const int which=-1) const
  {
    categoryT_cit cit=category.find(cat);
    if (cit==category.end())
      {
	if (verbose>=1) {
	  if (which<0) std::cout<<"Default value assigned to "<<cat<<"."<<what<<" : '"<<def<<"'.\n";
	  else std::cout<<"Default value assigned to "<<cat<<"."<<what<<"["<<which<<"] : '"<<def<<"'.\n";
	}
	return def;
      }

    paramsT &p=*(cit->second);
    paramsT_it it=p.find(what);

    if (it==p.end())
      {
	if (verbose>=1) {
	  if (which<0) std::cout<<"Default value assigned to "<<cat<<"."<<what<<" : '"<<def<<"'.\n";
	  else std::cout<<"Default value assigned to "<<cat<<"."<<what<<"["<<which<<"] : '"<<def<<"'.\n";
	}
	return def;
      }

    if ((index>=0)&&(it->second.val.size()<=index))
      {
	if (verbose>=1) {
	  if (which<0) std::cout<<"Default value assigned to "<<cat<<"."<<what<<" : '"<<def<<"'.\n";
	  else std::cout<<"Default value assigned to "<<cat<<"."<<what<<"["<<which<<"] : '"<<def<<"'.\n";
	}
	return def;
      }

    int id = (index<0)?0:index;

    if ((which>=0)&&(it->second.val[id].size()<=which))
      {	
	if (verbose>=1) {
	  if (which<0) std::cout<<"Default value assigned to "<<cat<<"."<<what<<" : '"<<def<<"'.\n";
	  else std::cout<<"Default value assigned to "<<cat<<"."<<what<<"["<<which<<"] : '"<<def<<"'.\n";
	}
	return def;
      }

    int w=(which<0)?0:which;
    T result;
    if (!it->second.val[id].size()) result=true;
    else from_string<T>(result,it->second.val[id][w]);
    
    it->second.wasRead[id]=true;

    if (verbose>=2)
      {
	if (which<0) std::cout<<"New value assigned to "<<cat<<"."<<what<<" : '"<<result<<"'.\n";
	else std::cout<<"New value assigned to "<<cat<<"."<<what<<"["<<which<<"] : '"<<result<<"'.\n";
      }
    return result;
  }

  template <typename T>
  T get(const std::string &what,const std::string &cat,const T &def, const int which=-1) const
  {
    return get<T>(-1,what,cat,def,which);
  }

   long isDefined(const std::string &what,const std::string &cat) const
  {
    categoryT_cit cit=category.find(cat);
    if (cit==category.end()) return 0;
   
    paramsT &p=*(cit->second);
    paramsT_cit it=p.find(what);

    if (it==p.end()) return 0;

    return it->second.val.size();
  }


};

#endif

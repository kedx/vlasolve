#ifndef __VAL_LOCATION_TYPE_HXX__
#define __VAL_LOCATION_TYPE_HXX__

#include "typeSelect.hxx"

struct valLocationVal {
  enum type {CELL=0, VERTEX=1, UNDEFINED=-1};
};

typedef valLocationVal::type valLocationType;

struct valLocationSelect : public typeSelect<valLocationVal> {
  valLocationSelect()
  {
    insert("cell",valLocationVal::CELL);
    insert("vertex",valLocationVal::VERTEX);
  }
  std::string name() {return "val location";}
};

#endif

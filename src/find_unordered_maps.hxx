#ifndef _FIND_UNORDERED_MAPS_HEADER_
#define _FIND_UNORDERED_MAPS_HEADER_

#ifndef HAVE_TR1

#include <map>
template<typename KeyType, typename MappedType>
struct my_unordered_map
{
  typedef std::map<KeyType, MappedType> type;
};

#else //HAVE_TR1

#ifdef HAVE_TR1_HEADER_PREFIX
#include <tr1/unordered_map>
#else //HAVE_TR1_HEADER_PREFIX
#include <unordered_map>
#endif //HAVE_TR1_HEADER_PREFIX

#ifdef HAVE_TR1_NAMESPACE_PREFIX
template<typename KeyType, typename MappedType>
struct my_unordered_map
{
  typedef std::tr1::unordered_map<KeyType, MappedType> type;
};
#else //HAVE_TR1_NAMESPACE_PREFIX
template<typename KeyType, typename MappedType>
struct my_unordered_map
{
  typedef std::unordered_map<KeyType, MappedType> type;
};
#endif //HAVE_TR1_NAMESPACE_PREFIX

#endif //HAVE_TR1


#endif

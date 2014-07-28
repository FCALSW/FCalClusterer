#ifndef SortingClass_h
#define SortingClass_h 1

/* --------------------------------------------------------------------------
   class ....
   -------------------------------------------------------------------------- */
class SortingClass {

public:
  SortingClass(int idNow, double weightNow);
  ~SortingClass();

  int	Id;
  double	Weight;
};

SortingClass::SortingClass(int idNow, double weightNow):  
  Id (idNow),
  Weight(weightNow)
{
}

SortingClass::~SortingClass() {
}

/* --------------------------------------------------------------------------
   sort rule
   -------------------------------------------------------------------------- */
//in descending order (lowest value is first)
inline bool cmpRuleDesc( SortingClass const& a , SortingClass const& b ) {
  return a.Weight < b.Weight;
}

#endif // SortingClass_h

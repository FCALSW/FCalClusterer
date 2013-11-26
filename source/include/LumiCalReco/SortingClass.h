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
bool cmpRuleDesc( SortingClass * a , SortingClass * b ) {

  return a->Weight < b->Weight;
}

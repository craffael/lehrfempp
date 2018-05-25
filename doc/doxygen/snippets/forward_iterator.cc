
#include <lf/base/base.h>
#include <chrono>

using namespace lf::base;

class Ingredient {};

//! [recipe]
class Recipe {
public:

  /* The preparation time of the recipe */
  virtual std::chrono::duration<long> prepTime() = 0;

  /* Points to the first ingredient and can be incremented */
  virtual ForwardIterator<Ingredient> beginIngredients() = 0;

  /* Points to one behind the last ingredient */
  virtual ForwardIterator<Ingredient> endIngredients() = 0;

  /* virtual destructor */
  virtual ~Recipe() {}
};
//! [recipe]

void foo() {
  
  //! [usage]
  std::vector<std::string> vector{"hello", "world"};
  std::list<std::string> list{"Gugus"};

  // assign any iterator to the wrapper and dereference it
  ForwardIterator<std::string> fi = vector.begin();
  assert(*fi == "hello");

  // assign a list iterator to the wrapper
  fi = list.begin();
  assert(*fi == "Gugus");
  assert(fi->length() == 5);

  // increment the iterator:
  ++fi;

  // compare it to another iterator:
  ForwardIterator<std::string> fi2 = list.end();
  assert(fi == fi2);
  //! [usage]

  //! [equality]
  std::vector<int> numbers{0,1,2};
  std::vector<int>::iterator it = numbers.begin();
  std::vector<int>::const_iterator const_it = numbers.begin();

  assert(ForwardIterator<const int>(it) != ForwardIterator<const int>(const_it));
  //! [equality]


}
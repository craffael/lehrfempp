/**
 * This code is a walkthrough for the usage of static variables in C++
 * (1) Static Variables
 * (2) Macros
 * (3) Track-Class
 * (4) Combination of the above / The Greater Picture
 * (5) Summary
 */

# include <iostream>

/**
 * (1) Static Variables
 * These are global variables, defined outside of any function.
 */

// Ordinary variables. They will be accessible in any function defined later on.
static int a; // a = 0
static double x = 9.2;

// Classes may have static variables as well. These member variables will have
// the same value across all instances of that class.
class A {
  public:
    int individual;
    static int same;
};

// Static member variables have to be defined outside of the class,
// otherwise the compilter doesn't allocate memory for it.
int A::same = 0;

void part1(){
  std::cout << "\n-- Part 1: Static variables ----------------------------\n\n";
  A Alice, Bob;
  Alice.individual = 100;
  Bob.individual = 3;
  std::cout << "Alice.individual = " << Alice.individual << "\n"
            << "Bob.individual   = " << Bob.individual << "\n\n";

  Alice.same = 200; // For both Alice & Bob the variable `same` is 200.
  Bob.same = 6; // Overrides Alice's setting, for both now `same` is 6!
  std::cout << "Alice.same = " << Alice.same << "\n"
            << "Bob.same   = " << Bob.same << "\n\n";
}

/**
 * (2) Macros
 * Macros are preprocessor instructions, that get evaluated during compile time.
 * It can be used for declaring variables, replace parts of your code or
 * for writing functions (see (4)).
 */
# define X 1
# define Y "Hello there!"
# define Z                        \
  struct reindeer {               \
    std::string name = "Rudolph"; \
    unsigned int mph = 200;       \
  };

Z // define the class reindeer

// Or we can declare functions:
// as one-liner
# define F(a, b) ((a + b) % 5)
// or more complex
# define G(a, b) int myadd() { int c = (a + b) % 5; \
                               return c - 1;        \
                             }
// Define myadd with a = X = 1 and b = 27.
// Note, that (a) we omit the semicolon, the compiler adds it automatically
// and (b) that this cannot be inside a function, since G is a "placeholder"
// for the full function definition of myadd.
G(X, 27)

// Some syntactic sugar..
// -- `#msg` stringifies: msg -> "msg"
# define H(msg) void mymsg() { std::cout << #msg << "\n"; }
H(barbershop)
// -- `a##b` concatenates a and b
# define I(a,b) int my##a() { return b; }
I(price,4000) // defines a function myprice that returns 4000

void part2() {
  std::cout << "\n-- Part 2: Macros --------------------------------------\n\n";
  std::cout << Y << "\n";
  std::cout << "X + 2               = " << X + 2 << "\n";
  std::cout << "F(X, 27)            = " << F(X, 27) << "\n";
  std::cout << "myadd() [a=X, b=27] = " << myadd() << "\n\n";

  reindeer r;
  std::cout << r.name << " goes " << r.mph << " mph!\n\n";

  std::cout << "I went to a [mymsg()] ";
  mymsg();
  std::cout << "and that cost me [myprice()]" << myprice() << " CHF.\n";
}

/**
 * (3) Track-Class
 * Basically a linked list with a stored value per element.
 * The root element is the one that's added the latest.
 */
template <class T>
class Track {
  public:
   Track *next_;         // next variable in list
   Track *prev_;         // previous variable in list
   std::string name_;    // name of variable
   std::string comment_; // optional documentation
   T &ref_;              // stored value

   /**
    * Constructor: Places a new item of the global info list in the list
    * Usually called via a macro (DECLARE, COUNTER).
    * The second version of the constructor also permits to add a comment to the
    * information item.
    */
   Track(std::string name, T &ref, Track<T> *&root,
         std::string comment = std::string({}));
   Track() = delete;      // tells compiler to not create a default constructor
   Track(const Track &) = delete;
   Track(Track &&) = delete;
   Track &operator=(const Track &) = delete;
   Track &operator=(Track &&) = delete;
   ~Track();

   static int count;      // global counter (= #instances of Track)
};

template <class T>
int Track<T>::count = 0;

template <class T>
Track<T>::Track(std::string name, T &ref, Track<T> *&root, std::string comment)
     : prev_(nullptr),
       name_(std::move(name)),
       comment_(std::move(comment)),
       ref_(ref) {
   /**
    * Create new element in the list. The new element will replace the root
    * element. Increases the static variable `counter`.
    *
    *      old    ---- prev_ --->   new
    *      root   <--- next_ ----   root (this element)
    */
   if (root) {
     root->prev_ = this;  // previous element of old root is this
   }
   next_ = root;          // next element is the old root
   root = this;           // this element is the new root
   count++;
}

template <class T>
Track<T>::~Track() {
   if (prev_) {
     prev_->next_ = next_;
   }
   if (next_) {
     next_->prev_ = prev_;
   }
   count--;
}

// Track<int> will function as a static variable. Each variable will know
// what the next and the previous one is, have a name and a value.
using StaticVar = Track<int>;

// define a global root element (we need an initial one)
StaticVar *root = nullptr;

void part3() {
  std::cout << "\n-- Part 3: Track-Class ---------------------------------\n\n";

  int local_i = 20;
  StaticVar i("global_i", local_i, root, "first variable");
  // now `i` is the new root

  int local_j = 40,
      local_k = 80;
  StaticVar j("global_j", local_j, root, "second variable");
  StaticVar k("global_k", local_k, root, "third variable");

  std::cout << "#(StaticVars): " << i.count << "(i), "
                                 << j.count << "(j), "
                                 << k.count << "(j)\n\n";
  std::cout << "Name / Value\n"
            << i.name_ << " / " << i.ref_ << "\n"
            << j.name_ << " / " << j.ref_ << "\n"
            << k.name_ << " / " << k.ref_ << "\n\n";

  const StaticVar *iter = root;
  std::cout << "Starting at root (should be k) and iterating through \
next elements until a nullptr (the end) is reached.\n";
  std::cout << "Name / Value\n";
  while (iter != nullptr) {
    std::cout << iter->name_ << " / " << iter->ref_ << "\n";
    iter = iter->next_;
  }
}

/**
 * (4) The Greater Picture
 * Now we want to simplify the declaration of global variables.
 * We define two Macro functions, one for declaring the variables, one for
 * setting it (with an ordinary function, comes at the end).
 */

// Setting a new list of variables
StaticVar *ctrl_root = nullptr;

// goal: a new static variable called ctrlvar<value> with the name `name`
// and value 0.
// e.g. CONTROLDECLARE(f, "variable f")
//  --> ctrlvarf with ctrlvarf.ref_ == 0 and ctrlvarf.name_ == "variable f"
# define CONTROLDECLARE(value, name)                            \
  int value = 0;                                                \
  static StaticVar ctrlvar##value(name, value, ctrl_root)

// Note, that all values are set to 0!
CONTROLDECLARE(f, "variable f"); // new static variable, called ctrlvarf
CONTROLDECLARE(g, "variable g");
CONTROLDECLARE(h, "variable h"); // current root element

// What if we embed everything into a class?
class Var {
  public:
    static int u, v, w; // should be static so we set the variables globally
};

# define CLASSCONTROLDECLARE(class, value, name)                \
  int class::value = 0;                                         \
  static StaticVar class##value(name, class::value, ctrl_root)

// Does two things:
// (1) Sets Var::u to 0
// (2) Defines a new StaticVar variable called Varu with value 0 and
//     name "class variable u"
CLASSCONTROLDECLARE(Var, u, "class variable u");
CLASSCONTROLDECLARE(Var, v, "class variable v");
CLASSCONTROLDECLARE(Var, w, "class variable w");
// Error: x is not a member of Var, so Var::x makes no sense!
// CLASSCONTROLDECLARE(Var, x, "class variable x?")

// Setting a variable directly is easy, see the first lines of part4().
// But how do we set a variable, if we only know its name ("variable f")?
// Since we have a linked list of all static variables, we simply iterate
// through all elements and if we found the right one we set its value
bool SetCtrlVar(const std::string &name, int value, const StaticVar *root) {
  const StaticVar *iter = root;
  // iterate over all static variables
  while (iter != nullptr) {
    // are we at the right one?
    if (iter->name_ == name) {
      iter->ref_ = value;
      return true;  // done!
    }
    iter = iter->next_; // try the next one
  }
  return false; // couldn't find it..
}

void part4() {
  std::cout << "\n-- Part 4: The Greater Picture -------------------------\n\n";

  ctrlvarh.ref_ = 101;
  Var::w = -52;
  // Varw.ref_ = -52; // already set! Since StaticVar took Var::w via reference.

  std::cout << "Initally everyhing is 0 (except Var::w and ctrlvarh, which \
we already set manually.)\n"
            << "Variable / Name / Value\n"
            << "ctrlvarf / " << ctrlvarf.name_ << " / " << ctrlvarf.ref_ << "\n"
            << "ctrlvarg / " << ctrlvarg.name_ << " / " << ctrlvarg.ref_ << "\n"
            << "ctrlvarh / " << ctrlvarh.name_ << " / " << ctrlvarh.ref_ << "\n"
            << "Varu / " << Varu.name_ << " / " << Varu.ref_ << "\n"
            << "Varv / " << Varv.name_ << " / " << Varv.ref_ << "\n"
            << "Varw / " << Varw.name_ << " / " << Varw.ref_ << "\n\n";
  std::cout << "We could access the values of the class directly (Var::u). \
Note, however, that they do not have a name.\n"
            << "Variable / Value\n"
            << "Var::u / " << Var::u << "\n"
            << "Var::v / " << Var::v << "\n"
            << "Var::w / " << Var::w << "\n\n";

  std::cout << "Let's set some variables using their names:\n";
  SetCtrlVar("variable f", 9, ctrl_root);
  std::cout << "SetCtrlVar('variable f', 9, ctrl_root);\n"
            << "--> ctrlvarf / " << ctrlvarf.name_ << " / " << ctrlvarf.ref_
            << "\n\n";
  std::cout << "This works for class values as well, since Varu and Var::u are \
connected via reference:\n";
  SetCtrlVar("class variable u", 2147483647, ctrl_root);
  std::cout << "SetCtrlVar('class variable u', 2147483647, ctrl_root);\n"
            << "--> Varu / " << Varu.name_ << " / " << Varu.ref_ << "\n"
            << "    Var::u / " << Var::u << "\n\n";

  std::cout << "-- To summarise:\n"
            << "Declare a new static (int) variable `i` with name 'my i': "
            << "CONTROLDECLARE(i, 'my i');\n"
            << "Set the variable to the value 701: "
            << "\t SetCtrlVar('my i', 701);\n"
            << "Alternatively, directly: "
            << "\t\t ctrlvari = 701;\n\n"
            << "If there is a class Var with a static int member u, set it via:"
            << "\t CLASSCONTROLDECLARE(Var, u, 'my u');\n"
            << "Set it via: "
            << "\t SetCtrlVar('my u', 702); \n"
            << "Or: "
            << "\t\t Var::u = 702; \n"
            << "Or: "
            << "\t\t Varu = 702; \n";
}


/**
 * Run everything :)
 */
int main() {
  part1();
  part2();
  part3();
  part4();
  return 0;
}
